"""fc_summary.py - Gate FC1 자체 시험 결과 JSON 들을 한 표로 모은다 (보고서용, 읽기만).

사용: python3 -I -B fc_summary.py <selftest 결과 폴더>
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


def main():
    d = sys.argv[1]
    rows = []
    for name in ("fc1.json", "fc2_a.json", "fc2_b.json", "fc3.json", "fcx.json"):
        r = load(d, name)
        if r is None:
            rows.append((name, "MISSING"))
            continue
        combos = r.get("combos", [])
        boxes = [c.get("boxes") for c in combos if c.get("boxes") is not None]
        rows.append((name, "cases=%d mismatch=%d intervals=%d combos=%d min_boxes_per_combo=%s "
                     "t_orig=%.1fs t_fast=%.1fs wall=%.1fs precs=%s orders=%s" % (
                         r["cases_total"], r["mismatch_total"], r["intervals_compared"], len(combos),
                         min(boxes) if boxes else "-", r.get("t_orig_s", 0), r.get("t_fast_s", 0),
                         r.get("wall_s", 0), r.get("precs"), r.get("orders"))))
        if r.get("first_mismatches"):
            rows.append((name + " first", repr(r["first_mismatches"][:2])[:800]))
        if name == "fcx.json":
            for k, v in r.get("by_case", {}).items():
                rows.append(("  fcx " + k, "n=%d mismatch=%d kinds=%s" % (v["n"], v["mismatch"], v["kinds"])))
    r = load(d, "fc0.json")
    if r is None:
        rows.append(("fc0.json", "MISSING"))
    else:
        rows.append(("fc0.json", "ops=%d mismatch=%d int_roundtrip=%d mismatch=%d wall=%.1fs" % (
            r["n_ops"], r["ops_mismatch"], r["int_roundtrip_n"], r["int_roundtrip_mismatch"], r["wall_s"])))
        rows.append(("fc0 op kinds", "%d kinds, min count %d" % (len(r["counts"]), min(r["counts"].values()))))
    for a, b in rows:
        print("%-28s %s" % (a, b))


if __name__ == "__main__":
    main()
