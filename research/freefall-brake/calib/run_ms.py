"""점선 사격 실행: python3 calib/run_ms.py <c2|c3|n1> [--fast] [--workers N] [--retries R]"""
import argparse, json, os, sys, time
from fractions import Fraction as F
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
ap = argparse.ArgumentParser()
ap.add_argument("which", choices=["c2", "c3", "n1"])
ap.add_argument("--fast", action="store_true")
ap.add_argument("--workers", type=int, default=1)
ap.add_argument("--retries", type=int, default=2)
ap.add_argument("--r-exp", type=int, default=30, help="box radius 2^-r_exp")
ap.add_argument("--k", type=int, default=8, help="segments (1 = single shooting)")
a = ap.parse_args()
if a.fast:
    sys.path.insert(0, os.path.join(ROOT, "..", "pythagorean-3body", "fastcore"))
    import fast_system as FS
from ff_common import dyad  # noqa: E402
from ff_brake import brake_problem  # noqa: E402
from ff_mshoot import run_ms  # noqa: E402
import ff_check  # noqa: E402
ORBITS = {"c2": ("c2_f2_111", [1, 1, 1]), "c3": ("c3_f1_0608_1", [F(3, 5), F(4, 5), 1]),
          "n1": ("n1_f1_1_08_08", [1, F(4, 5), F(4, 5)])}
dyad.set_prec(128)
tag, masses = ORBITS[a.which]
pre = json.load(open(os.path.join(ROOT, "results", "prescan_c2.json"), encoding="utf-8"))[a.which]["refined"]
times = [F(t) for t in json.load(open(os.path.join(ROOT, "results", "ms_times.json"), encoding="utf-8"))[a.which]] if a.k > 1 else []
if a.k not in (1, 8):
    raise SystemExit("only k=1 or k=8 (ms_times.json has 8-segment stations)")
prob = brake_problem(tag, masses, system=FS.make_nbody_fast(masses) if a.fast else None)
R = F(1, 2 ** a.r_exp)
rec = run_ms(prob, times, [F(pre["x"]), F(pre["y"])], F(pre["tau"]), R, R, R, retries=a.retries, workers=a.workers)
rec["independent_check"] = ff_check.check(rec) if rec["verdict"] == "PASS" else None
rec["fast"], rec["workers"] = a.fast, a.workers
rec["r_exp"], rec["k_req"] = a.r_exp, a.k
json.dump(rec, open(os.path.join(ROOT, "results", "ms_%s_k%d_r%d.json" % (tag, a.k, a.r_exp)), "w", encoding="utf-8"), indent=1)
print(tag, rec["verdict"], rec["reason"], rec.get("segments"), "attempts", rec["attempts"], "wall_s=%.0f" % rec["wall_s"])
