"""시간 측정(비판정): 점 초기조건을 보폭 제어로 τ까지, 최대 WALL초. python calib/time_smoke.py n1 600"""
import json, os, sys, time
from fractions import Fraction as F
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ff_common import dyad, ts, ls  # noqa: E402
from ff_flow import _controlled_step  # noqa: E402
MASSES = {"c2": [1, 1, 1], "c3": [F(3, 5), F(4, 5), 1], "n1": [1, F(4, 5), F(4, 5)]}
which, wall = sys.argv[1], float(sys.argv[2])
dyad.set_prec(128)
pre = json.load(open(os.path.join(ROOT, "results", "prescan_c2.json"), encoding="utf-8"))[which]["refined"]
s = ts.make_nbody(MASSES[which])
st = ls.lohner_init(s, [F(-1, 2), 0, F(1, 2), 0, F(pre["x"]), F(pre["y"])] + [0] * 6)
tau = F(pre["tau"]); t0 = time.time(); n = 0
while st.t < tau and time.time() - t0 < wall:
    st, _ = _controlled_step(s, st, min(F(1, 8), tau - st.t)); n += 1
    if n % 20 == 0:
        print("steps %d t=%.4f / %.4f  elapsed %.0fs" % (n, float(st.t), float(tau), time.time() - t0), flush=True)
print("DONE" if st.t == tau else "WALL", "steps", n, "t=%.5f" % float(st.t), "elapsed %.0fs" % (time.time() - t0))
