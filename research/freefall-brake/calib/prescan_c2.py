"""C2 사전 조사(비엄밀): python calib/prescan_c2.py (research/freefall-brake 에서) → results/prescan_c2.json"""
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ff_prescan import prescan, refine  # noqa: E402

# Li & Liao 2019 표 A1, A2, A5 (x, y, T)
ORBITS = {
    "c2": ([1, 1, 1], 0.2053886532, 0.1952668419, 1.6896364928),
    "c3": ([0.6, 0.8, 1], 0.2061730906, 0.4463433325, 2.6925803371),
    "n1": ([1, 0.8, 0.8], 0.0009114239, 0.3019805958, 1.8286248401),
}
out = {}
for k, (m, x, y, T) in ORBITS.items():
    out[k] = {"table": {"masses": m, "x": x, "y": y, "T": T},
              "prescan": prescan(m, x, y, T), "refined": refine(m, x, y, T / 2)}
    print(k, out[k]["prescan"], out[k]["refined"])
    if out[k]["prescan"]["vmax_half"] > 1e-6:
        print("STOP: vmax_half > 1e-6 for", k)
with open(os.path.join(ROOT, "results", "prescan_c2.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, indent=1)
