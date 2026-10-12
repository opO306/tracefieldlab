"""서버 점검(비판정): 케플러 점선 사격을 workers=1, 3으로 돌려 결과가 같은지. python3 calib/ms_par_check.py"""
import math, os, sys
from fractions import Fraction as F
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ff_common import dyad, ts  # noqa: E402
from ff_problem import Problem  # noqa: E402
from ff_mshoot import run_ms  # noqa: E402
dyad.set_prec(128)
R = F(1, 2 ** 40)
p = Problem("kep_ms", ts.make_kepler(1), [F(3, 2), 0, 0, 0], [3], [0, 1], [F(-1, 2), 0], F(1, 8), tol_control=True)
a = run_ms(p, [F(1), F(2)], [F(1 / math.sqrt(3))], F(math.pi), R, R, R, workers=1)
b = run_ms(p, [F(1), F(2)], [F(1 / math.sqrt(3))], F(math.pi), R, R, R, workers=3)
keys = [k for k in a if k != "wall_s"]
print("IDENTICAL" if all(a[k] == b[k] for k in keys) else "DIFFERENT", a["verdict"], b["verdict"],
      "wall %.0f vs %.0f" % (a["wall_s"], b["wall_s"]))
