"""C2 실행: python calib/run_c2.py <c2|c3|n1> [--fast]  (research/freefall-brake 에서)"""
import argparse, json, os, sys, time
from fractions import Fraction as F

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
ap = argparse.ArgumentParser()
ap.add_argument("which", choices=["c2", "c3", "n1"])
ap.add_argument("--fast", action="store_true")
a = ap.parse_args()
if a.fast:  # fast_system 은 원본 모듈을 먼저 import 해야 한다(다른 경로로 먼저 import 되면 ImportError)
    sys.path.insert(0, os.path.join(ROOT, "..", "pythagorean-3body", "fastcore"))
    import fast_system as FS  # 리눅스 서버 전용; 실패하면 여기서 종료
from ff_common import dyad  # noqa: E402
from ff_brake import brake_problem  # noqa: E402
from ff_problem import run_problem  # noqa: E402
import ff_check  # noqa: E402

ORBITS = {"c2": ("c2_f2_111", [1, 1, 1]),
          "c3": ("c3_f1_0608_1", [F(3, 5), F(4, 5), 1]),
          "n1": ("n1_f1_1_08_08", [1, F(4, 5), F(4, 5)])}
dyad.set_prec(128)
tag, masses = ORBITS[a.which]
pre = json.load(open(os.path.join(ROOT, "results", "prescan_c2.json"), encoding="utf-8"))[a.which]
system = FS.make_nbody_fast(masses) if a.fast else None
prob = brake_problem(tag, masses, system=system)
m = [F(pre["refined"]["x"]), F(pre["refined"]["y"]), F(pre["refined"]["tau"])]
t0 = time.time()
rec = run_problem(prob, m, [F(1, 2 ** 30)] * 3)
rec["independent_check"] = ff_check.check(rec) if rec["verdict"] == "PASS" else None
rec["prescan"] = pre
rec["fast"] = a.fast
with open(os.path.join(ROOT, "results", tag + ".json"), "w", encoding="utf-8") as f:
    json.dump(rec, f, indent=1)
print(tag, rec["verdict"], rec.get("reason", ""), "wall_s=%.0f" % (time.time() - t0))
