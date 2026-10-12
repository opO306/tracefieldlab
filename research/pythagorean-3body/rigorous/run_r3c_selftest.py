"""R3-S3 자기시험 실행기: 참조 구현(오탐) + 가짜 구현(적발)을 check_r3c.py 로 각각 별도 프로세스에서 검사.

사용법: python run_r3c_selftest.py [--workers 6] [--s1-steps 1000] [--only 이름,...] [--seed N]
결과: r3c_selftest_runs/R3C_CHECK_<이름>.json, R3C_CHECK_SELFTEST_RESULT.json
참조 구현을 먼저 돌려 기준값 캐시를 만들고, 가짜들은 그 캐시를 읽기 전용으로 쓴다(기준값은 구현과 무관).
"""
import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, "r3c_selftest_runs")
CHECKER = os.path.join(HERE, "check_r3c.py")
CACHE = os.path.join(RUNS, "ref_cache.pkl")

REQUIRED = ["inv_float", "point_J", "no_remainder", "R_not_updated", "no_center_offset", "no_apriori"]   # 과제가 요구한 6종
EXTRA = ["box_ignores_B", "h_used_lie", "poly_drops_last_term", "J_drops_last_term", "qt_inverse"]
# 결함이 아닌 유효한 변형(오탐 시험): 위반이 나오면 검사기가 틀린 것이다
VALID_VARIANTS = ["J_no_remainder", "PY_is_identity", "pivot_qr"]


def sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def run_one(name, seed, s1_steps, write_cache, extra):
    mod = os.path.join(HERE, "fake_r3c_%s.py" % name)
    out = os.path.join(RUNS, "R3C_CHECK_%s_seed%d.json" % (name, seed))
    log = os.path.join(RUNS, "R3C_CHECK_%s_seed%d.log" % (name, seed))
    cmd = [sys.executable, "-I", CHECKER, mod, "--out", out, "--seed", str(seed), "--s1-steps", str(s1_steps),
           "--ref-cache", CACHE] + (["--ref-cache-write"] if write_cache else []) + extra
    t0 = time.time()
    with open(log, "w", encoding="utf-8") as f:
        p = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
    return name, p.returncode, time.time() - t0, out


def summarize(name, out):
    d = json.load(open(out, encoding="utf-8"))
    crit = d.get("criteria", {})
    viol = {k: v["counts"]["VIOLATION"] for k, v in crit.items()}
    tot = {k: v["total"] for k, v in crit.items()}
    sealed = [k for k, v in crit.items() if v["sealed"]]
    return {
        "name": name,
        "verdict_sealed": d.get("verdict_sealed_criteria", d.get("verdict")),
        "verdict_all": d.get("verdict_all_including_supplementary"),
        "violations": viol, "totals": tot,
        "detected_by_sealed": any(viol[k] > 0 for k in sealed),
        "detected_any": any(viol.values()),
        "undecided": {k: v["counts"]["UNDECIDED"] for k, v in crit.items() if v["counts"]["UNDECIDED"]},
        "refused": {k: v["counts"]["REFUSED"] for k, v in crit.items() if v["counts"]["REFUSED"]},
        "not_run": {k: v["counts"]["NOT_RUN"] for k, v in crit.items() if v["counts"]["NOT_RUN"]},
        "min_counterexample": d.get("min_counterexample"),
        "elapsed_sec": d.get("elapsed_sec"),
        "checker_sha256": d.get("checker_sha256"),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--s1-steps", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=20261011)
    ap.add_argument("--only", default=None)
    ap.add_argument("--extra", default="", help="check_r3c.py 에 넘길 추가 인자(공백 구분)")
    a = ap.parse_args()
    os.makedirs(RUNS, exist_ok=True)
    extra = a.extra.split() if a.extra else []
    names_fake = REQUIRED + EXTRA + VALID_VARIANTS
    if a.only:
        sel = [x.strip() for x in a.only.split(",")]
        names_fake = [n for n in names_fake if n in sel]
    results = {}
    t0 = time.time()
    if not a.only or "correct" in a.only.split(","):
        n, rc, dt, out = run_one("correct", a.seed, a.s1_steps, True, extra)
        results[n] = summarize(n, out)
        results[n].update(returncode=rc, wall_sec=round(dt, 1))
        print("correct:", results[n]["verdict_sealed"], results[n]["violations"], "%.0fs" % dt, flush=True)
    with ThreadPoolExecutor(max_workers=a.workers) as ex:
        futs = [ex.submit(run_one, n, a.seed, a.s1_steps, False, extra) for n in names_fake]
        for fu in futs:
            n, rc, dt, out = fu.result()
            try:
                results[n] = summarize(n, out)
                results[n].update(returncode=rc, wall_sec=round(dt, 1))
                print("%-22s sealed=%s any=%s viol=%s %.0fs" % (n, results[n]["verdict_sealed"], results[n]["detected_any"],
                                                              {k: v for k, v in results[n]["violations"].items() if v}, dt), flush=True)
            except Exception as e:  # noqa: BLE001
                results[n] = {"name": n, "error": repr(e), "returncode": rc}
                print(n, "ERROR", e, flush=True)
    required = [results[n] for n in REQUIRED if n in results and "error" not in results[n]]
    extras = [results[n] for n in EXTRA if n in results and "error" not in results[n]]
    valids = [results[n] for n in VALID_VARIANTS if n in results and "error" not in results[n]]
    corr = results.get("correct")
    summary = {
        "checker": "check_r3c.py", "checker_sha256": sha(CHECKER), "seed": a.seed, "s1_steps": a.s1_steps,
        "false_positives_reference": (sum(v for v in corr["violations"].values()) if corr else None),
        "reference_undecided_refused_notrun": ({"undecided": corr["undecided"], "refused": corr["refused"],
                                                "not_run": corr["not_run"]} if corr else None),
        "required_fakes_total": len(required),
        "required_fakes_detected": sum(1 for r in required if r["detected_any"]),
        "required_fakes_detected_by_sealed_criteria": sum(1 for r in required if r["detected_by_sealed"]),
        "extra_fakes_total": len(extras),
        "extra_fakes_detected": sum(1 for r in extras if r["detected_any"]),
        "extra_fakes_detected_by_sealed_criteria": sum(1 for r in extras if r["detected_by_sealed"]),
        "undetected": [r["name"] for r in required + extras if not r["detected_any"]],
        "valid_variants_total": len(valids),
        "valid_variants_flagged_as_violation": [r["name"] for r in valids if r["detected_any"]],
        "per_fake": results,
        "wall_sec": round(time.time() - t0, 1),
    }
    with open(os.path.join(HERE, "R3C_CHECK_SELFTEST_RESULT.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=True, indent=1)
    print("required detected %d/%d (sealed-only %d); extra detected %d/%d (sealed-only %d); reference violations=%s; undetected=%s; valid variants flagged=%s"
          % (summary["required_fakes_detected"], summary["required_fakes_total"],
             summary["required_fakes_detected_by_sealed_criteria"],
             summary["extra_fakes_detected"], summary["extra_fakes_total"], summary["extra_fakes_detected_by_sealed_criteria"],
             summary["false_positives_reference"], summary["undetected"], summary["valid_variants_flagged_as_violation"]))


if __name__ == "__main__":
    main()
