#!/usr/bin/env python3
"""R1-S5 자기시험 실행기: 모든 fake_r1_*.py를 check_r1.py에 넣어 적발 여부를 모은다.

정답 구현(fake_r1_correct.py)은 오탐이 0이어야 하고, 나머지 가짜는 전부 FAIL이어야 한다.
검사기는 가짜마다 별도 프로세스(python -I)로 돌린다(dyad 정밀도 전역 상태 분리).

사용법: python run_r1_selftest.py [--starts N] [--workers K] [--seeds 20261011,7,424242]
"""
import argparse
import concurrent.futures as cf
import glob
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
CHECKER = os.path.join(HERE, "check_r1.py")
OUTDIR = os.path.join(HERE, "r1_selftest_runs")


def run_one(module, tag, extra, timeout):
    os.makedirs(OUTDIR, exist_ok=True)
    out = os.path.join(OUTDIR, "R1_CHECK_%s.json" % tag)
    cmd = [sys.executable, "-I", CHECKER, module, "--out", out] + extra
    t0 = time.time()
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, encoding="utf-8", errors="replace")
        rc = p.returncode
        tail = (p.stdout or "")[-1500:] + (p.stderr or "")[-800:]
    except subprocess.TimeoutExpired:
        rc, tail = -9, "TIMEOUT"
    res = None
    if os.path.exists(out):
        with open(out, encoding="utf-8") as f:
            res = json.load(f)
    return tag, module, rc, time.time() - t0, res, tail


def summarize(res):
    if not res or "criteria" not in res:
        return None
    crit = {}
    for k, c in res["criteria"].items():
        crit[k] = {"verdict": c["verdict"], "violations": c["counts"]["VIOLATION"], "total": c["total"],
                   "refused": c["counts"]["REFUSED"], "undecided": c["counts"]["UNDECIDED"],
                   "not_run": c["counts"]["NOT_RUN"]}
    return crit


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--starts", type=int, default=50)
    ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--seeds", default="20261011,7,424242")
    ap.add_argument("--timeout", type=float, default=3600)
    ap.add_argument("--fakes-only", action="store_true")
    ap.add_argument("--correct-only", action="store_true")
    args = ap.parse_args()
    fakes = sorted(p for p in glob.glob(os.path.join(HERE, "fake_r1_*.py"))
                   if not p.endswith("fake_r1_correct.py"))
    correct = os.path.join(HERE, "fake_r1_correct.py")
    jobs = []
    base = ["--starts", str(args.starts), "--time-limit", str(args.timeout - 60)]
    if not args.fakes_only:
        for i, sd in enumerate(args.seeds.split(",")):
            extra = base + ["--seed", sd] + ([] if i == 0 else ["--no-ladder"])
            jobs.append((correct, "correct_seed%s" % sd, extra))
    if not args.correct_only:
        for fk in fakes:
            name = os.path.splitext(os.path.basename(fk))[0]
            jobs.append((fk, name, base + ["--no-ladder"]))
    results = []
    with cf.ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = [ex.submit(run_one, m, t, e, args.timeout) for (m, t, e) in jobs]
        for f in cf.as_completed(futs):
            tag, module, rc, dt, res, tail = f.result()
            sm = summarize(res)
            print("[%6.1fs] %-34s rc=%s sealed=%s" % (dt, tag, rc, (res or {}).get("verdict_sealed_criteria")), flush=True)
            results.append({"tag": tag, "module": os.path.basename(module), "returncode": rc,
                            "elapsed_sec": round(dt, 1),
                            "verdict_sealed": (res or {}).get("verdict_sealed_criteria"),
                            "verdict_all": (res or {}).get("verdict_all_including_supplementary"),
                            "criteria": sm, "min_counterexample": (res or {}).get("min_counterexample"),
                            "tail": tail if res is None else None})
    results.sort(key=lambda r: r["tag"])
    correct_rows = [r for r in results if r["tag"].startswith("correct_seed")]
    fake_rows = [r for r in results if not r["tag"].startswith("correct_seed")]
    caught_sealed = [r["tag"] for r in fake_rows if r["verdict_sealed"] == "FAIL"]
    caught_any = [r["tag"] for r in fake_rows if r["verdict_all"] == "FAIL"]
    only_supp = [t for t in caught_any if t not in caught_sealed]
    false_pos = [r["tag"] for r in correct_rows if r["verdict_all"] == "FAIL"]
    not_clean = [r["tag"] for r in correct_rows if r["verdict_sealed"] != "PASS" or r["verdict_all"] != "PASS"]
    summary = {
        "starts": args.starts, "fakes_total": len(fake_rows),
        "caught_by_sealed_criteria": len(caught_sealed), "caught_by_any_including_supplementary": len(caught_any),
        "caught_only_by_supplementary": only_supp,
        "missed": [r["tag"] for r in fake_rows if r["verdict_all"] != "FAIL"],
        "correct_runs": len(correct_rows), "false_positives(FAIL)": false_pos,
        "correct_runs_not_all_PASS": not_clean,
    }
    out = {"summary": summary, "rows": results}
    path = os.path.join(HERE, "R1_CHECK_SELFTEST_RESULT.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=True, indent=1, default=str)
    print(json.dumps(summary, ensure_ascii=True, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
