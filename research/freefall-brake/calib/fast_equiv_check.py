"""서버 점검(비판정): fastcore 경로와 원본 파이썬 경로가 삼체 보폭 제어 걸음 N개에서 비트 단위로 같은지.
python3 calib/fast_equiv_check.py <n1|c2|c3> <N>"""
import json, os, sys
from fractions import Fraction as F
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "..", "pythagorean-3body", "fastcore"))
import fast_system as FS  # noqa: E402
from ff_common import dyad, ts, ls  # noqa: E402
from ff_flow import _controlled_step  # noqa: E402
MASSES = {"c2": [1, 1, 1], "c3": [F(3, 5), F(4, 5), 1], "n1": [1, F(4, 5), F(4, 5)]}
which, N = sys.argv[1], int(sys.argv[2])
dyad.set_prec(128)
pre = json.load(open(os.path.join(ROOT, "results", "prescan_c2.json"), encoding="utf-8"))[which]["refined"]
x0 = [F(-1, 2), 0, F(1, 2), 0, F(pre["x"]), F(pre["y"])] + [0] * 6
out = []
for s in (ts.make_nbody(MASSES[which]), FS.make_nbody_fast(MASSES[which])):
    st = ls.lohner_init(s, x0)
    for _ in range(N):
        st, d = _controlled_step(s, st, F(1, 8))
    out.append((st.t, st.c, st.B, [(r.lo_q, r.hi_q) for r in st.R], [[(e.lo_q, e.hi_q) for e in row] for row in d["J"]]))
print("IDENTICAL" if out[0] == out[1] else "DIFFERENT", "t=%s" % float(out[0][0]), FS.fc_stats(s))
