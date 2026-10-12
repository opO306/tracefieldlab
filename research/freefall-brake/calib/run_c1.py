"""C1 교정 실행: 결과 JSON 저장 + 독립 검사. 실행: python calib/run_c1.py (research/freefall-brake 에서)"""
import json
import math
import os
import sys
from fractions import Fraction as F

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ff_common import dyad, ts  # noqa: E402
from ff_problem import Problem, run_problem  # noqa: E402
import ff_check  # noqa: E402

dyad.set_prec(128)
R = F(1, 2 ** 40)
RUNS = [
    ("c1a_harmonic",
     Problem("C1a_harmonic", ts.make_harmonic(), [F(1), F(0)], [], [1], [F(0)], F(1, 8)),
     [F(math.pi)], [R]),
    ("c1b_kepler",
     Problem("C1b_kepler", ts.make_kepler(1), [F(3, 2), F(0), F(0), F(0)], [3], [0, 1],
             [F(-1, 2), F(0)], F(1, 64)),
     [F(1 / math.sqrt(3)), F(math.pi)], [R, R]),
]
os.makedirs(os.path.join(ROOT, "results"), exist_ok=True)
for tag, prob, m, r in RUNS:
    rec = run_problem(prob, m, r)
    rec["independent_check"] = ff_check.check(rec) if rec["verdict"] == "PASS" else None
    with open(os.path.join(ROOT, "results", tag + ".json"), "w", encoding="utf-8") as f:
        json.dump(rec, f, indent=1)
    ic = rec["independent_check"]["verdict"] if rec["independent_check"] else "-"
    print(tag, rec["verdict"], rec.get("reason", ""), "indep:", ic, "wall_s=%.1f" % rec["wall_s"])
