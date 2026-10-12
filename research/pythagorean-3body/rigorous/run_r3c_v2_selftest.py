"""R3c 검사기 v2 자기시험 실행기.

v1 자기시험(run_r3c_selftest.py)과 같은 참조 구현·가짜 구현(fake_r3c_*.py)을 쓰되, 각 구현을
"계약 엄수 래퍼"로 감싸서 검사한다. 래퍼는 lohner_step 에 들어오는 state 의 c·B 항과 R 끝점이
이진 분수(가수 <= 현재 정밀도)가 아니면 실제 구현이 낸 것과 같은 문구의 ValueError 를 던진다
(R0 결정 A). 이렇게 해야 v1 이 만든 계약 밖 시작 집합이 자기시험에서도 드러난다
(v1 의 자기시험 때는 래퍼가 없어서 참조 구현이 분모 1000 의 중심을 그냥 받아 들였다).

래퍼의 lohner_init 은 비이진 x0(예: 9/10)를 받으면 가까운 이진 분수 중심 + 그 오차를 덮는 R 로
바꿔 준다(계약을 지키는 구현이 해야 할 일; 실제 구현이 어떻게 하는지는 보지 않았다).

사용법:
  python run_r3c_v2_selftest.py [--checker v2|v1] [--workers 2] [--only a,b,..] [--outdir DIR]
         [--extra "--only S2,X4"] [--s1-steps 1000] [--seed N] [--compare-v1 R3C_CHECK_SELFTEST_RESULT.json]
  --checker v1 은 재현 시험용: v1 검사기 × 계약 엄수 참조 구현에서 S2 8 + X4 2 가 REFUSED 인지 본다.
참조 구현('correct')을 먼저 돌려 기준값 캐시를 만든 뒤 가짜들은 그 캐시를 읽기 전용으로 쓴다.
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

REQUIRED = ["inv_float", "point_J", "no_remainder", "R_not_updated", "no_center_offset", "no_apriori"]
EXTRA = ["box_ignores_B", "h_used_lie", "poly_drops_last_term", "J_drops_last_term", "qt_inverse"]
VALID_VARIANTS = ["J_no_remainder", "PY_is_identity", "pivot_qr"]

WRAPPER = '''"""계약 엄수 래퍼: fake_r3c_{name} (자동 생성; 수정하지 말 것)."""
import os, sys
from fractions import Fraction
sys.path.insert(0, {here!r})
import dyad
from dyad import DI
import fake_r3c_{name} as _m
from fake_r3c_{name} import make_harmonic, make_kepler, make_nbody
import _fake_r3c_core as _core


def _dyadic_ok(q, prec):
    q = Fraction(q)
    d = q.denominator
    if d & (d - 1):
        return False
    n = abs(q.numerator)
    if n == 0:
        return True
    n >>= (n & -n).bit_length() - 1
    return n.bit_length() <= prec


def _check_state(state):
    prec = dyad.get_prec()
    entries = list(state.c) + [x for row in state.B for x in row] + [e for d in state.R for e in (d.lo_q, d.hi_q)]
    for q in entries:
        if not _dyadic_ok(q, prec):
            raise ValueError("state entry %s is not a dyadic with <= %d mantissa bits (R0 decision A: one precision per run)" % (Fraction(q), prec))


def lohner_init(system, x0):
    st = _m.lohner_init(system, x0)
    prec = dyad.get_prec()
    c = list(st.c)
    R = list(st.R)
    for i, q in enumerate(c):
        if not _dyadic_ok(q, prec):
            c2 = _core.round_dyadic(q, prec)
            R[i] = DI(Fraction(q) - c2)          # 오차를 R 로 바깥 반올림해 덮는다
            c[i] = c2
    st.c, st.R = c, R
    return st


def lohner_step(system, state, h, order):
    _check_state(state)
    return _m.lohner_step(system, state, h, order)


def enclosure_box(state):
    return _m.enclosure_box(state)
'''


def sha(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def make_wrapper(outdir, name):
    d = os.path.join(outdir, "strict")
    os.makedirs(d, exist_ok=True)
    p = os.path.join(d, "strict_fake_r3c_%s.py" % name)
    with open(p, "w", encoding="utf-8", newline="\n") as f:
        f.write(WRAPPER.format(name=name, here=HERE))
    return p


def run_one(checker, tag, outdir, cache, name, seed, s1_steps, write_cache, extra):
    mod = make_wrapper(outdir, name)
    out = os.path.join(outdir, "R3C_CHECK_%s_%s_seed%d.json" % (tag, name, seed))
    log = out[:-5] + ".log"
    cmd = [sys.executable, "-I", checker, mod, "--out", out, "--seed", str(seed), "--s1-steps", str(s1_steps),
           "--ref-cache", cache] + (["--ref-cache-write"] if write_cache else []) + extra
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
        "counts": {k: v["counts"] for k, v in crit.items()},
        "detected_by_sealed": any(viol[k] > 0 for k in sealed),
        "detected_any": any(viol.values()),
        "undecided": {k: v["counts"]["UNDECIDED"] for k, v in crit.items() if v["counts"]["UNDECIDED"]},
        "refused": {k: v["counts"]["REFUSED"] for k, v in crit.items() if v["counts"]["REFUSED"]},
        "not_run": {k: v["counts"]["NOT_RUN"] for k, v in crit.items() if v["counts"]["NOT_RUN"]},
        "min_counterexample": d.get("min_counterexample"),
        "elapsed_sec": d.get("elapsed_sec"),
        "checker_sha256": d.get("checker_sha256"),
        "module_sha256": d.get("module_sha256"),
        "start_set_contract_audit_v2": d.get("start_set_contract_audit_v2"),
        "error": d.get("error") if d.get("verdict") in ("CHECKER_ERROR", "LOAD_FAILED", "API_MISSING") else None,
    }


def compare_v1(results, v1_path):
    """v1 자기시험 요약(R3C_CHECK_SELFTEST_RESULT.json)과 기준별 위반·시도 수를 나란히 비교."""
    v1 = json.load(open(v1_path, encoding="utf-8"))["per_fake"]
    rows = {}
    for n, r in results.items():
        if "error" in r and "violations" not in r:
            continue
        o = v1.get(n)
        if not o or "violations" not in o:
            continue
        row = {}
        for k in sorted(set(o["violations"]) | set(r["violations"])):
            row[k] = {"v1_viol": o["violations"].get(k), "v2_viol": r["violations"].get(k),
                      "v1_total": o["totals"].get(k), "v2_total": r["totals"].get(k)}
        rows[n] = {"verdict_sealed_v1": o.get("verdict_sealed"), "verdict_sealed_v2": r.get("verdict_sealed"),
                   "same_detection_class": bool(o.get("detected_any")) == bool(r.get("detected_any")),
                   "per_criterion": row}
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--checker", default="v2", choices=["v1", "v2"])
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--s1-steps", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=20261011)
    ap.add_argument("--only", default=None, help="콤마 구분 이름(correct 포함 가능)")
    ap.add_argument("--extra", default="", help="검사기에 넘길 추가 인자(공백 구분), 예: '--only S2,X4'")
    ap.add_argument("--outdir", default=os.path.join(HERE, "r3c_v2_selftest_runs"))
    ap.add_argument("--cache", default=None, help="기준값 캐시 경로(기본: outdir/ref_cache.pkl)")
    ap.add_argument("--no-cache-write", action="store_true", help="correct 도 캐시를 쓰지 않는다(읽기만)")
    ap.add_argument("--reuse-correct", action="store_true",
                    help="outdir 에 같은 seed 의 correct 결과 JSON 이 있으면 다시 돌리지 않고 그것을 요약에 쓴다")
    ap.add_argument("--compare-v1", default=None)
    ap.add_argument("--summary-name", default=None)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    checker = os.path.join(HERE, "check_r3c_v2.py" if a.checker == "v2" else "check_r3c.py")
    tag = a.checker
    cache = a.cache or os.path.join(a.outdir, "ref_cache.pkl")
    extra = a.extra.split() if a.extra else []
    names_fake = REQUIRED + EXTRA + VALID_VARIANTS
    sel = [x.strip() for x in a.only.split(",")] if a.only else None
    if sel:
        names_fake = [n for n in names_fake if n in sel]
    results = {}
    t0 = time.time()
    cj = os.path.join(a.outdir, "R3C_CHECK_%s_correct_seed%d.json" % (tag, a.seed))
    if a.reuse_correct and os.path.exists(cj):
        results["correct"] = summarize("correct", cj)
        results["correct"].update(returncode=None, wall_sec=None, reused_existing_result=True)
        print("correct (reused %s):" % cj, results["correct"]["verdict_sealed"], results["correct"]["violations"],
              results["correct"]["refused"], flush=True)
    elif sel is None or "correct" in sel:
        n, rc, dt, out = run_one(checker, tag, a.outdir, cache, "correct", a.seed, a.s1_steps, not a.no_cache_write, extra)
        results[n] = summarize(n, out)
        results[n].update(returncode=rc, wall_sec=round(dt, 1))
        print("correct:", results[n]["verdict_sealed"], results[n]["violations"], results[n]["refused"], "%.0fs" % dt, flush=True)
    with ThreadPoolExecutor(max_workers=a.workers) as ex:
        futs = [ex.submit(run_one, checker, tag, a.outdir, cache, n, a.seed, a.s1_steps, False, extra) for n in names_fake]
        for fu in futs:
            n, rc, dt, out = fu.result()
            try:
                results[n] = summarize(n, out)
                results[n].update(returncode=rc, wall_sec=round(dt, 1))
                print("%-22s sealed=%s any=%s viol=%s refused=%s %.0fs" % (
                    n, results[n]["verdict_sealed"], results[n]["detected_any"],
                    {k: v for k, v in results[n]["violations"].items() if v}, results[n]["refused"], dt), flush=True)
            except Exception as e:  # noqa: BLE001
                results[n] = {"name": n, "error": repr(e), "returncode": rc}
                print(n, "ERROR", e, flush=True)
    ok = lambda names: [results[n] for n in names if n in results and "violations" in results[n]]  # noqa: E731
    required, extras, valids = ok(REQUIRED), ok(EXTRA), ok(VALID_VARIANTS)
    corr = results.get("correct")
    summary = {
        "checker": os.path.basename(checker), "checker_sha256": sha(checker), "seed": a.seed, "s1_steps": a.s1_steps,
        "wrapper": "strict-contract wrapper (ValueError on non-dyadic / over-precision state entries) around each fake",
        "extra_args": extra,
        "false_positives_reference": (sum(v for v in corr["violations"].values()) if corr and "violations" in corr else None),
        "reference_undecided_refused_notrun": ({"undecided": corr["undecided"], "refused": corr["refused"],
                                                "not_run": corr["not_run"]} if corr and "violations" in corr else None),
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
    if a.compare_v1:
        summary["compare_with_v1_selftest"] = compare_v1(results, a.compare_v1)
    sname = a.summary_name or ("R3C_CHECK_%s_SELFTEST_RESULT.json" % tag.upper())
    with open(os.path.join(a.outdir, sname), "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=True, indent=1)
    print("required detected %d/%d (sealed-only %d); extra detected %d/%d (sealed-only %d); reference violations=%s; "
          "undetected=%s; valid variants flagged=%s"
          % (summary["required_fakes_detected"], summary["required_fakes_total"],
             summary["required_fakes_detected_by_sealed_criteria"],
             summary["extra_fakes_detected"], summary["extra_fakes_total"], summary["extra_fakes_detected_by_sealed_criteria"],
             summary["false_positives_reference"], summary["undetected"], summary["valid_variants_flagged_as_violation"]))


if __name__ == "__main__":
    main()
