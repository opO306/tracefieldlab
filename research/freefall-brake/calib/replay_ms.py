"""PASS 기록의 마지막 시도를 같은 중심·상자로 다시 계산해 gm, Dg, K가 비트 단위로 같은지 비교(교차 확인).
python3 calib/replay_ms.py <record.json> <c2|c3|n1|kepler> [--fast] [--workers N]"""
import argparse, json, os, sys
from fractions import Fraction as F
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
ap = argparse.ArgumentParser()
ap.add_argument("record"); ap.add_argument("which", choices=["c2", "c3", "n1", "kepler"])
ap.add_argument("--fast", action="store_true"); ap.add_argument("--workers", type=int, default=1)
a = ap.parse_args()
if a.fast:
    sys.path.insert(0, os.path.join(ROOT, "..", "pythagorean-3body", "fastcore"))
    import fast_system as FS
from ff_common import dyad, DI, ts  # noqa: E402
from ff_brake import brake_problem  # noqa: E402
from ff_problem import Problem  # noqa: E402
from ff_mshoot import evaluate_ms  # noqa: E402
from ff_krawczyk import krawczyk, di2s  # noqa: E402

dyad.set_prec(128)
rec = json.load(open(a.record, encoding="utf-8"))
ORB = {"c2": ("c2_f2_111", [1, 1, 1]), "c3": ("c3_f1_0608_1", [F(3, 5), F(4, 5), 1]),
       "n1": ("n1_f1_1_08_08", [1, F(4, 5), F(4, 5)])}
if a.which == "kepler":
    prob = Problem("kep", ts.make_kepler(1), [F(3, 2), 0, 0, 0], [3], [0, 1], [F(-1, 2), 0], F(1, 8), tol_control=True)
else:
    tag, masses = ORB[a.which]
    prob = brake_problem(tag, masses, system=FS.make_nbody_fast(masses) if a.fast else None)
X = [DI(F(lo), F(hi)) for lo, hi in rec["X"]]
m = [F(v) for v in rec["m"]]
times = [F(t) for t in rec["times"]]
G, Dg, P, segs = evaluate_ms(prob, times, X, m, a.workers)
kr = krawczyk(G, Dg, X, m)
same = {"gm": [di2s(g) for g in G] == rec["gm"], "Dg": [[di2s(e) for e in r] for r in Dg] == rec["Dg"],
        "K": [di2s(k) for k in kr["K"]] == rec["K"], "verdict": kr["verdict"] == rec["verdict"]}
print("IDENTICAL" if all(same.values()) else "DIFFERENT", same, "verdict", kr["verdict"], "fast", a.fast, flush=True)
