"""star_dynamic_falsify.py 결과(JSONL) 요약. 실행: python summarize_dynamic.py out_*.jsonl"""
import json
import sys

CHECKS = ["C1", "C2", "C3", "C4", "C5"]
for path in sys.argv[1:]:
    rows = [json.loads(l) for l in open(path, encoding="utf-8")]
    meta = rows[0]["meta"]
    rows = rows[1:]
    cls = meta["class"]
    rel = [r for r in rows if r["reliable"]]
    print(f"== {cls}  ({path})  표본 {len(rows)} / 계획 {meta['n']}, 수치 신뢰(에너지 오차<=1e-9, rmin>=1e-3) {len(rel)}")
    kinds = {}
    for r in rows:
        kinds[r["kind"]] = kinds.get(r["kind"], 0) + 1
    print("   표본 종류:", kinds, "| 종료 상태:", {s: sum(1 for r in rows if r["status"] == s) for s in set(r["status"] for r in rows)})
    for subset, name in ((rows, "전체"), (rel, "신뢰")):
        if not subset:
            continue
        viol = {k: sum(1 for r in subset if k in r["first_violation_t"]) for k in CHECKS + ["fallback"]}
        mins = {k: min(r["min_margins"][k] for r in subset) for k in CHECKS}
        print(f"   [{name}] 위반 표본 수:", viol)
        print(f"   [{name}] 최소 여유:", {k: f"{v:.3g}" for k, v in mins.items()})
    if rows:
        print("   최대 에너지 상대오차: %.2e, 최소 rmin: %.3g, t_end 범위: %.1f-%.1f, 벽시계 합 %.0fs" % (
            max(r["maxdE_rel"] for r in rows), min(r["rmin"] for r in rows),
            min(r["t_end"] for r in rows), max(r["t_end"] for r in rows), sum(r["wall"] for r in rows)))
        print("   초기 여유 C 범위: %.3g - %.3g, r0/d1 최대 %.6f, rho0/(lam d1)-1 최소 %.3g" % (
            min(r["C"] for r in rows), max(r["C"] for r in rows), max(r["r0"] / r["d1"] for r in rows),
            min(r["rho0"] / (5 / 9 * r["d1"]) - 1 for r in rows)))
    ex = [r for r in rel if r["first_violation_t"]]
    if ex:
        r = min(ex, key=lambda r: min(r["first_violation_t"].values()))
        print("   가장 이른 위반(신뢰 표본): i=%d kind=%s r0=%.6f d1=%.6f rho0=%.6f rhod0=%.6f C=%.4g 위반시각=%s" % (
            r["i"], r["kind"], r["r0"], r["d1"], r["rho0"], r["rhod0"], r["C"],
            {k: round(v, 3) for k, v in r["first_violation_t"].items()}))
        print("      y0 =", [round(v, 12) for v in r["y0"]])
    print()
