"""PROOF_STAR_FINAL_{EN,KO}.md 에 새로 들어간 식의 점검 (sympy/Fraction). 증명 아님.

F1  Levi-Civita 비고: z = u^2, dt/ds = |u|^2, u'' = (h/2) u + (1/2)|u|^2 conj(u) f,
    h = (2|u'|^2 - k)/|u|^2  이면  d^2z/dt^2 = -k z/|z|^3 + f.
F2  |u'|^2 = (k + h |u|^2)/2  (h 의 정의를 다시 쓴 것)
F3  Yoshida(1972) Theorem 1 의 조건 (15) 우변 2H/mu + 2 m1 m2/(mu d) 는 d = m1 m2/K 에서 0.
F4  문서에 적은 피타고라스 수치(정확 유리수와 소수).
실행: python -I -X utf8 final_doc_checks.py
"""
from fractions import Fraction

import sympy as sp

ok_all = True


def report(name, ok):
    global ok_all
    ok_all &= bool(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}")


# ---------------------------------------------------------------- F1, F2
s = sp.symbols("s", real=True)
k = sp.symbols("k", positive=True)
fr, fi = sp.symbols("f_r f_i", real=True)
a = sp.Function("a")(s)
b = sp.Function("b")(s)
u = a + sp.I * b
ub = a - sp.I * b
u1 = sp.diff(u, s)
ub1 = sp.diff(ub, s)
absu2 = u * ub
f = fr + sp.I * fi
h = (2 * u1 * ub1 - k) / absu2
upp_rule = h / 2 * u + sp.Rational(1, 2) * absu2 * ub * f          # u''
ubpp_rule = h / 2 * ub + sp.Rational(1, 2) * absu2 * u * (fr - sp.I * fi)  # conj(u'')
z = u ** 2
# d/dt = (1/|u|^2) d/ds
zdot = sp.diff(z, s) / absu2
zddot = sp.diff(zdot, s) / absu2
# a'' = Re(u''), b'' = Im(u'') 를 u'', conj(u'') 규칙으로 대입
A2 = (upp_rule + ubpp_rule) / 2
B2 = (upp_rule - ubpp_rule) / (2 * sp.I)
zddot = zddot.subs({sp.diff(a, s, 2): A2, sp.diff(b, s, 2): B2})
target = -k * u ** 2 / (absu2 ** 3) + f   # -k z/|z|^3 + f,  |z| = |u|^2
resid = sp.simplify(sp.together(sp.expand(zddot - target)))
report("F1 Levi-Civita: d2z/dt2 = -k z/|z|^3 + f", resid == 0)
zdot_abs2 = sp.simplify(zdot * (2 * ub1 / u))  # |zdot|^2 = (2u'/ub)(2ub'/u)
h_def = sp.Rational(1, 2) * zdot_abs2 - k / absu2
report("F2a h = |zdot|^2/2 - k/|z|", sp.simplify(h_def - h) == 0)
report("F2b |u'|^2 = (k + h|u|^2)/2", sp.simplify(u1 * ub1 - (k + h * absu2) / 2) == 0)

# 음성 대조: u'' 의 1/2 를 1 로 바꾸면 F1 이 깨져야 한다
upp_bad = h / 2 * u + absu2 * ub * f
ubpp_bad = h / 2 * ub + absu2 * u * (fr - sp.I * fi)
zb = sp.diff(sp.diff(z, s) / absu2, s) / absu2
zb = zb.subs({sp.diff(a, s, 2): (upp_bad + ubpp_bad) / 2, sp.diff(b, s, 2): (upp_bad - ubpp_bad) / (2 * sp.I)})
report("F1-neg (일부러 틀린 계수는 잡혀야 함)", sp.simplify(sp.together(sp.expand(zb - target))) != 0)

# ---------------------------------------------------------------- F3
m1, m2, mu, K = sp.symbols("m1 m2 mu K", positive=True)
H = -K
d1 = m1 * m2 / K
report("F3 2H/mu + 2 m1 m2/(mu d) = 0 at d = m1 m2/K", sp.simplify(2 * H / mu + 2 * m1 * m2 / (mu * d1)) == 0)

# ---------------------------------------------------------------- F4
M0, M1, M2 = Fraction(3), Fraction(4), Fraction(5)
M = M0 + M1 + M2
lam, nu = M2 / (M1 + M2), M1 / (M1 + M2)
Hp = -(M0 * M1 / 5 + M0 * M2 / 4 + M1 * M2 / 3)
Kp = -Hp
d1p = M1 * M2 / Kp
vals = {"M": M, "lam": lam, "nu": nu, "mu": M0 * (M1 + M2) / M, "m": M1 * M2 / (M1 + M2),
        "H": Hp, "K": Kp, "d1": d1p, "lam d1": lam * d1p, "nu d1": nu * d1p, "2M": 2 * M}
for kk, vv in vals.items():
    print(f"   {kk} = {vv} = {float(vv):.12f}")
expect = {"M": Fraction(12), "lam": Fraction(5, 9), "nu": Fraction(4, 9), "mu": Fraction(9, 4), "m": Fraction(20, 9),
          "H": Fraction(-769, 60), "K": Fraction(769, 60), "d1": Fraction(1200, 769),
          "lam d1": Fraction(2000, 2307), "nu d1": Fraction(1600, 2307), "2M": Fraction(24)}
report("F4 피타고라스 값(문서 기재값과 일치)", all(vals[kk] == expect[kk] for kk in expect))
# 초기 조건: 질량 3 (1,3), 4 (-2,-1), 5 (1,-1)
P = {3: (1, 3), 4: (-2, -1), 5: (1, -1)}
d2 = lambda i, j: (P[i][0] - P[j][0]) ** 2 + (P[i][1] - P[j][1]) ** 2
report("F4b 거리 |3-4|=5, |3-5|=4, |4-5|=3, 무게중심 원점",
       (d2(3, 4), d2(3, 5), d2(4, 5)) == (25, 16, 9) and sum(m * P[m][0] for m in P) == 0 and sum(m * P[m][1] for m in P) == 0)

print("요약:", "모두 PASS" if ok_all else "FAIL 있음")
