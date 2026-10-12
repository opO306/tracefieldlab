"""Exact range of total momentum and centre of mass (linear forms, exact on a box)
over the saved R4 boxes. True values are exactly 0 for the Pythagorean problem.
Run from research/pythagorean-3body: python analysis_width/an4_conservation.py
"""
import json
import math
from fractions import Fraction as F

m = [3, 4, 5]


def lg(x):
    return math.log10(float(x)) if x > 0 else float("-inf")


B = [json.loads(l) for l in open("results_r4/r4_boxes_128w.jsonl")]
for b in B[::50] + [B[-1]]:
    lo = [F(x) for x in b["lo"]]
    hi = [F(x) for x in b["hi"]]
    w = max(h - l for l, h in zip(lo, hi))
    widths = []
    for base in (6, 0):  # 6: velocities (momentum), 0: positions (centre of mass)
        for ax in (0, 1):
            c = [0] * 12
            for i in range(3):
                c[base + 2 * i + ax] = F(m[i])
            a = sum(ci * (l if ci > 0 else h) for ci, l, h in zip(c, lo, hi))
            z = sum(ci * (h if ci > 0 else l) for ci, l, h in zip(c, lo, hi))
            widths.append((z - a) / 12)
    print(f"t={float(F(b['t'])):6.2f} log10 box width={lg(w):7.2f} "
          "range width/12 (px py cx cy):", " ".join(f"{lg(x):7.2f}" for x in widths))
