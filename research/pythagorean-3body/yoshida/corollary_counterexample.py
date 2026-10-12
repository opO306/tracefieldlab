"""Yoshida(1972) Corollary 의 진술(조건 (24)-(27))만으로는 결론 r<m1 m2/K 가 보장되지 않음을 보이는 상태 점검.

반례 후보: 피타고라스 궤적의 t≈17.12 상태(배정밀도 DOP853). 궤적의 정확성과 무관하게,
'이 상태 자체'를 초기조건으로 보고 그 상태의 에너지 H 로 K 를 다시 계산해 조건을 검사한다.
결론 'r<m1 m2/K (t>=0)' 은 t=0 (=이 상태)에서 이미 r 값과 비교하면 된다.
mpmath 50자리로 float 상태값을 그대로 유리수처럼 취급해 계산(반올림 영향 배제).
"""
import numpy as np, mpmath as mp
from scipy.integrate import solve_ivp
import importlib.util, os
spec = importlib.util.spec_from_file_location("yc", os.path.join(os.path.dirname(os.path.abspath(__file__)), "yoshida_check.py"))
yc = importlib.util.module_from_spec(spec); spec.loader.exec_module(yc)

mp.mp.dps = 50
s = yc.integrate(1e-13, 1e-15, tmax=17.2)
y = s.sol(17.118)
Y = [mp.mpf(float(z)) for z in y]           # 이 float 값들을 정확한 초기조건으로 채택
mm = [mp.mpf(3), mp.mpf(4), mp.mpf(5)]
x = [(Y[0], Y[1]), (Y[2], Y[3]), (Y[4], Y[5])]
v = [(Y[6], Y[7]), (Y[8], Y[9]), (Y[10], Y[11])]
# 질량중심 운동 제거 여부와 무관하도록 상대좌표만 사용
T = sum(mm[b] * (v[b][0]**2 + v[b][1]**2) for b in range(3)) / 2
Mt = sum(mm)
Vc = [sum(mm[b] * v[b][c] for b in range(3)) / Mt for c in range(2)]
T -= Mt * (Vc[0]**2 + Vc[1]**2) / 2
U = -sum(mm[i] * mm[j] / mp.sqrt((x[i][0]-x[j][0])**2 + (x[i][1]-x[j][1])**2) for i, j in [(0,1),(0,2),(1,2)])
H = T + U; K = -H
m0, m1, m2 = mm
lam = m2/(m1+m2); nu = m1/(m1+m2); M = Mt
rv = (x[2][0]-x[1][0], x[2][1]-x[1][1]); r = mp.sqrt(rv[0]**2 + rv[1]**2)
G = [(m1*x[1][c] + m2*x[2][c])/(m1+m2) for c in range(2)]
Gd = [(m1*v[1][c] + m2*v[2][c])/(m1+m2) for c in range(2)]
pv = [x[0][c]-G[c] for c in range(2)]; pd = [v[0][c]-Gd[c] for c in range(2)]
rho = mp.sqrt(pv[0]**2 + pv[1]**2); rhodot = (pv[0]*pd[0] + pv[1]*pd[1]) / rho
d = m1*m2/K
a = (lam + rho/r) / 2                       # lam < a <= rho/r 가 되도록 중간값 선택
print("H =", mp.nstr(H, 20), " (참값 -769/60 =", mp.nstr(mp.mpf(-769)/60, 20), ")")
print("lam =", mp.nstr(lam, 12), " a =", mp.nstr(a, 12))
print("(24) d = m1 m2/K =", mp.nstr(d, 12))
print("(25) rho >= a d :", mp.nstr(rho, 12), ">=", mp.nstr(a*d, 12), "->", rho >= a*d)
print("(26) rho/r >= a :", mp.nstr(rho/r, 12), ">=", mp.nstr(a, 12), "->", rho/r >= a)
need = mp.sqrt(2*M*(lam/(rho+nu*d) + nu/(rho-lam*d)))
print("(27) rhodot > sqrt(...):", mp.nstr(rhodot, 12), ">", mp.nstr(need, 12), "->", rhodot > need)
print("결론 'r < m1 m2/K' at t=0 :", mp.nstr(r, 12), "<", mp.nstr(d, 12), "->", r < d)
