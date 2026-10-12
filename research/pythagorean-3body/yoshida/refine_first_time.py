"""taylor_B 의 t=59.72 상태(고정밀)에서 DOP853 로 짧게(0.02) 적분하며 첫 성립 시각을 이분법으로 정밀화 (비엄밀)."""
import json, os, math, numpy as np
from scipy.integrate import solve_ivp
import importlib.util
here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("at", os.path.join(here, "analyze_taylor.py")); at = importlib.util.module_from_spec(spec); spec.loader.exec_module(at)
spec2 = importlib.util.spec_from_file_location("yc", os.path.join(here, "yoshida_check.py")); yc = importlib.util.module_from_spec(spec2); spec2.loader.exec_module(yc)
S = json.load(open(os.path.join(here, "taylor_B.json")))["samples"]
z = [float(q) for q in S[5972][1:]]
# taylor 순서(body별 x,y,vx,vy) -> yc 순서(x 6개, v 6개)
y0 = np.array([z[0],z[1],z[4],z[5],z[8],z[9], z[2],z[3],z[6],z[7],z[10],z[11]])
sol = solve_ivp(yc.rhs, (59.72, 59.74), y0, method="DOP853", rtol=1e-13, atol=1e-15, dense_output=True).sol
def jac(t):
    y = sol(t); rho, rhodot, r, _ = yc.jacobi(y); return rho, rhodot, r
dgrid = [0.002*i for i in range(1, 2501)]
def cstar(t):
    c = at.T1(*jac(t), at.d1); return c is not None and c[0] > 0
def t1opt(t):
    rho, rhodot, r = jac(t); return max((at.T1(rho, rhodot, r, d) or (-1e9,))[0] for d in dgrid) > 0
for name, f in [("C* (d=d1)", cstar), ("Theorem1 d-opt", t1opt)]:
    lo, hi = 59.72, 59.73
    assert not f(lo) and f(hi)
    for _ in range(30):
        mid = (lo+hi)/2
        if f(mid): hi = mid
        else: lo = mid
    print(f"{name}: first time ~ {hi:.6f}  (rho,rhodot,r)={tuple(round(q,5) for q in jac(hi))}")
