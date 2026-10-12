"""summarize_results.py - results/*.json 을 보고서용 표(마크다운)로 요약. 계산 없음.

사용: python3 summarize_results.py <results 폴더> > SUMMARY.md
"""
import json
import os
import sys


def load(d, name):
    p = os.path.join(d, name)
    if not os.path.exists(p):
        return None
    with open(p, encoding="utf-8") as f:
        return json.load(f)


def fmt_slack(s):
    return "—" if s is None else "%.1e" % s


def t2_table(name, d):
    if d is None:
        print("\n### %s: 결과 없음 (미완)\n" % name)
        return
    s = d["summary"]
    print("\n### %s — 사례 %d개, 걸음 예외 %d개, 소요 %.0f s" % (name, s["n_cases"], s["n_step_errors"], d.get("wall_s", 0)))
    print("| 검사 | PASS | UNDECIDED | VIOLATION | 최소 상대 여유 |")
    print("|---|---|---|---|---|")
    for k, v in s["by_check"].items():
        print("| %s | %d | %d | %d | %s |" % (k, v["PASS"], v["UNDECIDED"], v["VIOLATION"], fmt_slack(v["min_rel_slack"])))
    if s["step_errors"]:
        print("\n걸음 예외(구현이 상자를 내지 않고 멈춘 경우; 위반 아님):")
        for e in s["step_errors"][:10]:
            print("- " + e)
    hx = [r for r in d["results"] if "harness_exception" in r]
    if hx:
        print("\n시험 도구 예외 %d건:" % len(hx))
        for r in hx[:3]:
            print("- %s: %s" % (r["case"], r["harness_exception"].strip().splitlines()[-1]))
    red = [r for r in d["results"] if r.get("h_used") and r.get("h") != r.get("h_used")]
    print("\nh 축소(사전 상자) 사례: %d / %d" % (len(red), len(d["results"])))
    und = {}
    for r in d["results"]:
        for k, t in r.get("tallies", {}).items():
            if t["counts"]["UNDECIDED"]:
                und.setdefault(k, []).append(r["case"])
    for k, cs in und.items():
        print("- UNDECIDED가 나온 %s 검사 사례 %d개 (예: %s)" % (k, len(cs), cs[0]))


def t5_table(d, title):
    if d is None:
        print("\n### %s: 결과 없음 (미완)\n" % title)
        return
    print("\n### %s" % title)
    print("| 돌연변이 | 잡힘 | 진동자 위반(검사별) | 케플러 위반 | 삼체 위반 |")
    print("|---|---|---|---|---|")
    for name, m in d["mutants"].items():
        cells = []
        for su in ("harmonic", "kepler", "nbody"):
            ps = m["per_suite"].get(su)
            if ps is None:
                cells.append("—")
                continue
            v = {k: x["VIOLATION"] for k, x in ps["by_check"].items() if x["VIOLATION"]}
            e = ps["n_step_errors"]
            cells.append((", ".join("%s %d" % kv for kv in v.items()) or "0") + (" (예외 %d)" % e if e else ""))
        print("| %s | %s | %s |" % (name, "예" if m["caught"] else "**아니오**", " | ".join(cells)))
    if "original_control" in d:
        oc = d["original_control"]
        tot = sum(x["VIOLATION"] for s in oc.values() for x in s["by_check"].values())
        print("\n원본 대조(같은 quick 사례): 위반 합계 %d" % tot)


def t3_table(names, ds):
    print("\n### T3 사슬")
    print("| 사례 | 정지 | 걸음 | t_final | 점 수 | set PASS/UND/VIOL | box VIOL | aux PASS/UND/VIOL | set 최소 상대 여유 | 소요 s |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for n, d in zip(names, ds):
        if d is None:
            print("| %s | 결과 없음 | | | | | | | | |" % n)
            continue
        t = d["tallies"]
        s, b, a = t["set"]["counts"], t["box"]["counts"], t["aux"]["counts"]
        print("| %s | %s | %d | %.6f | %d | %d/%d/%d | %d | %d/%d/%d | %s | %.0f |" % (
            n, d["stop"], d["steps"], d["t_final"], d["n_points"], s["PASS"], s["UNDECIDED"], s["VIOLATION"],
            b["VIOLATION"], a["PASS"], a["UNDECIDED"], a["VIOLATION"], fmt_slack(t["set"]["min_rel_slack"]), d["wall_s"]))
        if d.get("curve"):
            last = d["curve"][-1]
            print("|  ↳ 마지막 기록 | k=%d t=%.4f | log10 외접폭 %.2f | log10 기준오차 %.1f | | | | | | |" % (
                last["k"], last["t"], last["wbox_log10"] or 0, last["eps_cum_log10"] or 0))


def main():
    d = sys.argv[1]
    t1 = load(d, "t1_inverse_format.json")
    print("## 자동 요약 (summarize_results.py)")
    if t1:
        print("\n### T1 역행렬·형식: %d개 중 통과 %d, 실패 %d" % (t1["summary"]["n_tests"], t1["summary"]["n_ok"], t1["summary"]["n_fail"]))
    for nm, f in (("T2 진동자", "t2_harmonic.json"), ("T2 케플러", "t2_kepler.json"), ("T2 삼체", "t2_nbody.json"),
                  ("T2-gap 삼체(검사기 미시험 부류)", "t2_nbody_gap.json")):
        t2_table(nm, load(d, f))
    t5_table(load(d, "t5_mutants.json"), "T5 돌연변이(검사력)")
    t5_table(load(d, "t5_valid_variants.json"), "T5 정당한 변형(위반 0이어야 맞음)")
    names = ["harm", "harm53", "kep", "nbody"]
    t3_table(names, [load(d, "t3_%s.json" % n) for n in names])
    srcs = set()
    for f in os.listdir(d):
        if f.endswith(".json"):
            x = load(d, f)
            if isinstance(x, dict) and "sources" in x:
                srcs.add(tuple(sorted((k, v["match_seal"]) for k, v in x["sources"].items())))
    print("\n원본 해시 확인(결과 파일들의 sources): %s" % sorted(srcs))


if __name__ == "__main__":
    main()
