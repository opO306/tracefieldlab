"""(★) 증명 독립 검토 - 기호 대수 점검 (sympy).

이 스크립트는 증명이 아니다. 증명 문서(PROOF_REVIEW_KO.md)의 각 대수 단계가
기호 계산으로도 같은 결과를 내는지 확인한다. 실행: python symbolic_checks.py
"""
from fractions import Fraction

import sympy as sp

results = []


def check(name, expr_zero):
    """expr_zero가 기호적으로 0이면 PASS."""
    z = sp.simplify(sp.together(sp.expand(expr_zero)))
    ok = (z == 0)
    results.append((name, ok))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + ("" if ok else f"  residual={z}"))
    return ok


m0, m1, m2 = sp.symbols("m0 m1 m2", positive=True)
M = m0 + m1 + m2
lam = m2 / (m1 + m2)          # G-m1 거리 = lam*r
nu = m1 / (m1 + m2)           # G-m2 거리 = nu*r
mu = m0 * (m1 + m2) / M       # 바깥 환산질량
mb = m1 * m2 / (m1 + m2)      # 쌍성 환산질량

# ---------------------------------------------------------------- S1 야코비 분해
dim = 3
x = [sp.Matrix(sp.symbols(f"x{i}_0:{dim}", real=True)) for i in range(3)]
v = [sp.Matrix(sp.symbols(f"v{i}_0:{dim}", real=True)) for i in range(3)]
mass = [m0, m1, m2]
Gpos = (m1 * x[1] + m2 * x[2]) / (m1 + m2)
Gvel = (m1 * v[1] + m2 * v[2]) / (m1 + m2)
rvec = x[2] - x[1]
rdot = v[2] - v[1]
rho = x[0] - Gpos
rhodot = v[0] - Gvel
Vcm = (m0 * v[0] + m1 * v[1] + m2 * v[2]) / M

T_cart = sum(sp.Rational(1, 2) * mass[i] * v[i].dot(v[i]) for i in range(3))
T_jac = (sp.Rational(1, 2) * M * Vcm.dot(Vcm) + sp.Rational(1, 2) * mb * rdot.dot(rdot)
         + sp.Rational(1, 2) * mu * rhodot.dot(rhodot))
check("S1a 운동에너지 = 1/2 M Vcm^2 + 1/2 m |rdot_vec|^2 + 1/2 mu |rhodot_vec|^2 (벡터)", T_cart - T_jac)
check("S1b x0 - x2 = rho_vec - nu r_vec  (r1 = |m2 m0|)", sp.simplify((x[0] - x[2]) - (rho - nu * rvec)).norm())
check("S1c x0 - x1 = rho_vec + lam r_vec (r2 = |m1 m0|)", sp.simplify((x[0] - x[1]) - (rho + lam * rvec)).norm())
check("S1d x2 - G = nu r_vec, x1 - G = -lam r_vec",
      sp.simplify((x[2] - Gpos) - nu * rvec).norm() + sp.simplify((x[1] - Gpos) + lam * rvec).norm())

# ---------------------------------------------------------------- S2 운동방정식
# 거리의 세제곱을 독립 기호로 둔다(식이 1/R^3에 대해 선형이므로 충분하다).
R3, R13, R23 = sp.symbols("R3 R13 R23", positive=True)   # |r|^3, r1^3=|x0-x2|^3, r2^3=|x0-x1|^3
a0 = -m1 * (x[0] - x[1]) / R23 - m2 * (x[0] - x[2]) / R13
a1 = m0 * (x[0] - x[1]) / R23 + m2 * (x[2] - x[1]) / R3
a2 = m0 * (x[0] - x[2]) / R13 + m1 * (x[1] - x[2]) / R3
rhoddot_newton = a0 - (m1 * a1 + m2 * a2) / (m1 + m2)
rhoddot_formula = -M * (lam / R13 + nu / R23) * rho + M * lam * nu * (1 / R13 - 1 / R23) * rvec
check("S2a 야코비 운동방정식 rho_ddot (뉴턴 직교좌표에서 유도와 일치)",
      sum(sp.simplify(e) ** 2 for e in (rhoddot_newton - rhoddot_formula)))
check("S2b mu*M*lam = m0*m2, mu*M*nu = m0*m1", (mu * M * lam - m0 * m2) ** 2 + (mu * M * nu - m0 * m1) ** 2)

# ---------------------------------------------------------------- S3 rho_ddot 분해
t = sp.symbols("t", real=True)
for d_ in (2, 3):
    fs = [sp.Function(f"p{k}")(t) for k in range(d_)]
    pv = sp.Matrix(fs)
    pd = pv.diff(t)
    pdd = pd.diff(t)
    rr = sp.sqrt(pv.dot(pv))
    lhs = sp.diff(rr, t, 2)
    if d_ == 3:
        cross2 = pv.cross(pd).dot(pv.cross(pd))
    else:
        cross2 = (pv[0] * pd[1] - pv[1] * pd[0]) ** 2
    rhs = pv.dot(pdd) / rr + cross2 / rr ** 3
    check(f"S3 ({d_}D) rho_ddot = (rho.rho_ddot)/rho + |rho x rho_dot|^2/rho^3", lhs - rhs)

# ---------------------------------------------------------------- S4 F = -dP/drho, (rho.rhoddot)/rho = -M F
L, rr_, rh, c, d = sp.symbols("lambda r rho c d", positive=True)
N_ = 1 - L
r1 = sp.sqrt(rh ** 2 + N_ ** 2 * rr_ ** 2 - 2 * N_ * rh * rr_ * c)
r2 = sp.sqrt(rh ** 2 + L ** 2 * rr_ ** 2 + 2 * L * rh * rr_ * c)
P = L / r1 + N_ / r2
F = L * (rh - N_ * rr_ * c) / r1 ** 3 + N_ * (rh + L * rr_ * c) / r2 ** 3
check("S4a F := lam(rho-nu r c)/r1^3 + nu(rho+lam r c)/r2^3 = -dP/drho", -sp.diff(P, rh) - F)
Mm = sp.symbols("M", positive=True)
radial = -Mm * (L / r1 ** 3 + N_ / r2 ** 3) * rh + Mm * L * N_ * (1 / r1 ** 3 - 1 / r2 ** 3) * rr_ * c
check("S4b (rho_vec.rho_ddot_vec)/rho = -M F  (rho_vec.r_vec = rho r c)", radial + Mm * F)

# ---------------------------------------------------------------- S5 르장드르 전개
a = sp.symbols("a", positive=True)
Gf = (1 - 2 * a * c + a ** 2) ** sp.Rational(-1, 2)
g = (1 - a * c) * (1 - 2 * a * c + a ** 2) ** sp.Rational(-3, 2)
check("S5a g(a,c) = G + a dG/da  (G = 생성함수)", Gf + a * sp.diff(Gf, a) - g)
NMAX = 14
cser = sp.symbols("c")
ser_G = sp.series(Gf.subs(c, cser), a, 0, NMAX + 1).removeO()
ser_g = sp.series(g.subs(c, cser), a, 0, NMAX + 1).removeO()
okG = all(sp.expand(ser_G.coeff(a, n) - sp.legendre(n, cser)) == 0 for n in range(NMAX + 1))
okg = all(sp.expand(ser_g.coeff(a, n) - (n + 1) * sp.legendre(n, cser)) == 0 for n in range(NMAX + 1))
results.append((f"S5b G = sum a^n P_n(c), n<= {NMAX}", okG))
results.append((f"S5c g = sum (n+1) a^n P_n(c), n<= {NMAX}", okg))
print(f"[{'PASS' if okG else 'FAIL'}] S5b 생성함수 계수 = P_n (n<={NMAX})")
print(f"[{'PASS' if okg else 'FAIL'}] S5c g 계수 = (n+1) P_n (n<={NMAX})")

# rho^2 F = lam g(nu s, c) + nu g(lam s, -c),  rho P = lam G(nu s,c) + nu G(lam s,-c), s=r/rho
s = sp.symbols("s", positive=True)
Fn = (L * (1 - N_ * s * cser) / (1 - 2 * N_ * s * cser + N_ ** 2 * s ** 2) ** sp.Rational(3, 2)
      + N_ * (1 + L * s * cser) / (1 + 2 * L * s * cser + L ** 2 * s ** 2) ** sp.Rational(3, 2))
Pn = (L / sp.sqrt(1 - 2 * N_ * s * cser + N_ ** 2 * s ** 2)
      + N_ / sp.sqrt(1 + 2 * L * s * cser + L ** 2 * s ** 2))
check("S5d rho^2 F(rho,r,c) = 정규화식 (r=s rho)",
      (rh ** 2 * F).subs(rr_, s * rh).subs(c, cser) - Fn)
check("S5e rho P(rho,r,c) = 정규화식 (r=s rho)",
      (rh * P).subs(rr_, s * rh).subs(c, cser) - Pn)
ser_F = sp.series(Fn, s, 0, 11).removeO()
ser_P = sp.series(Pn, s, 0, 11).removeO()
okA = okB = True
for n in range(11):
    An = (n + 1) * (L * N_ ** n + N_ * (-L) ** n)
    Bn = (L * N_ ** n + N_ * (-L) ** n)
    okA &= sp.simplify(sp.expand(ser_F.coeff(s, n) - An * sp.legendre(n, cser))) == 0
    okB &= sp.simplify(sp.expand(ser_P.coeff(s, n) - Bn * sp.legendre(n, cser))) == 0
results.append(("S5f rho^2 F 계수 = A_n P_n(c), A_n=(n+1)(lam nu^n + nu(-lam)^n), n<=10", okA))
results.append(("S5g rho P 계수 = B_n P_n(c),  B_n=(lam nu^n + nu(-lam)^n), n<=10", okB))
print(f"[{'PASS' if okA else 'FAIL'}] S5f rho^2 F 의 s^n 계수 = (n+1)(lam nu^n + nu(-lam)^n) P_n(c), n<=10")
print(f"[{'PASS' if okB else 'FAIL'}] S5g rho P 의 s^n 계수 = (lam nu^n + nu(-lam)^n) P_n(c), n<=10")
# 부호: (-1)^n (lam nu^n + nu(-lam)^n) = (-1)^n lam nu^n + nu lam^n
#   n 짝수: lam nu^n + nu lam^n > 0.  n 홀수: nu lam^n - lam nu^n = lam nu (lam^(n-1) - nu^(n-1)) >= 0 (lam>=nu)
nn = sp.symbols("n", positive=True, integer=True)
ok_sign = True
for n in range(1, 30, 2):
    expr = (-1) ** n * (L * N_ ** n + N_ * (-L) ** n) - L * N_ * (L ** (n - 1) - N_ ** (n - 1))
    ok_sign &= sp.expand(expr) == 0
results.append(("S5h 홀수 n: (-1)^n B_n = lam nu (lam^(n-1)-nu^(n-1)), n<30", ok_sign))
print(f"[{'PASS' if ok_sign else 'FAIL'}] S5h 홀수 n<30: (-1)^n B_n = lam nu (lam^(n-1) - nu^(n-1))")

# ---------------------------------------------------------------- S6 c=-1 값, W 미분, 단조성
W = L / (rh + N_ * d) + N_ / (rh - L * d)
Fbound = L / (rh + N_ * d) ** 2 + N_ / (rh - L * d) ** 2
check("S6a dW/drho = -(lam/(rho+nu d)^2 + nu/(rho-lam d)^2)", sp.diff(W, rh) + Fbound)
# c=-1 에서 r1=rho+nu r, r2=rho-lam r (rho>lam r 일 때 양수). 직접 대입:
r1m = rh + N_ * rr_
r2m = rh - L * rr_
check("S6b1 c=-1: r1^2 = (rho+nu r)^2",
      (rh ** 2 + N_ ** 2 * rr_ ** 2 + 2 * N_ * rh * rr_) - r1m ** 2)
check("S6b2 c=-1: r2^2 = (rho-lam r)^2",
      (rh ** 2 + L ** 2 * rr_ ** 2 - 2 * L * rh * rr_) - r2m ** 2)
F_m1 = L * (rh + N_ * rr_) / r1m ** 3 + N_ * (rh - L * rr_) / r2m ** 3
check("S6c F(r,c=-1) = lam/(rho+nu r)^2 + nu/(rho-lam r)^2 = Fbound(d=r)", F_m1 - Fbound.subs(d, rr_))
check("S6d P(r,c=-1) = W(rho,r)", (L / r1m + N_ / r2m) - W.subs(d, rr_))
check("S6e dW/dd = lam nu [ (rho-lam d)^-2 - (rho+nu d)^-2 ]",
      sp.diff(W, d) - L * N_ * ((rh - L * d) ** -2 - (rh + N_ * d) ** -2))
check("S6f dFbound/dd = 2 lam nu [ (rho-lam d)^-3 - (rho+nu d)^-3 ]",
      sp.diff(Fbound, d) - 2 * L * N_ * ((rh - L * d) ** -3 - (rh + N_ * d) ** -3))
check("S6g (rho+nu d) - (rho-lam d) = d", (rh + N_ * d) - (rh - L * d) - d)

# ---------------------------------------------------------------- S7 theta=pi 와 theta=0 비교(보조)
F_p1 = L / (rh - N_ * d) ** 2 + N_ / (rh + L * d) ** 2      # c=+1, r=d
check("S7 Fbound - F(c=+1) = 4 rho d lam nu [ (rho^2-lam^2 d^2)^-2 - (rho^2-nu^2 d^2)^-2 ]",
      Fbound - F_p1 - 4 * rh * d * L * N_ * ((rh ** 2 - L ** 2 * d ** 2) ** -2 - (rh ** 2 - N_ ** 2 * d ** 2) ** -2))

# ---------------------------------------------------------------- S9 보조정리 A 의 두 번째(미분) 증명 경로
check("S9a dP/dc = lam nu rho r (r1^-3 - r2^-3)", sp.diff(P, c) - L * N_ * rh * rr_ * (r1 ** -3 - r2 ** -3))
check("S9b r2^2 - r1^2 = r((lam-nu) r + 2 rho c)", (r2 ** 2 - r1 ** 2) - rr_ * ((L - N_) * rr_ + 2 * rh * c))
P_p1 = L / (rh - N_ * rr_) + N_ / (rh + L * rr_)     # c=+1 (rho > lam r)
check("S9c P(r,-1) - P(r,+1) = 2 lam nu r [ (rho^2-lam^2 r^2)^-1 - (rho^2-nu^2 r^2)^-1 ]",
      (L / r1m + N_ / r2m) - P_p1 - 2 * L * N_ * rr_ * ((rh ** 2 - L ** 2 * rr_ ** 2) ** -1 - (rh ** 2 - N_ ** 2 * rr_ ** 2) ** -1))
# 닫힌꼴: sum (n+1) s^n |beta_n| = lam/(1+nu s)^2 + nu/(1-lam s)^2,  sum s^n |beta_n| = lam/(1+nu s) + nu/(1-lam s)
abs_series_F = sum((n + 1) * s ** n * ((-1) ** n * L * N_ ** n + N_ * L ** n) for n in range(12))
abs_series_P = sum(s ** n * ((-1) ** n * L * N_ ** n + N_ * L ** n) for n in range(12))
okc = (sp.expand(sp.series(L / (1 + N_ * s) ** 2 + N_ / (1 - L * s) ** 2, s, 0, 12).removeO() - abs_series_F) == 0
       and sp.expand(sp.series(L / (1 + N_ * s) + N_ / (1 - L * s), s, 0, 12).removeO() - abs_series_P) == 0)
results.append(("S9d sum|A_n| , sum|B_n| 닫힌꼴 (s^11 까지)", okc))
print(f"[{'PASS' if okc else 'FAIL'}] S9d sum (n+1)s^n|beta_n| = lam/(1+nu s)^2 + nu/(1-lam s)^2, sum s^n|beta_n| = lam/(1+nu s)+nu/(1-lam s) (s^11 까지)")

# ---------------------------------------------------------------- S8 피타고라스 상수(정확 유리수)
M0, M1, M2 = Fraction(3), Fraction(4), Fraction(5)
H = Fraction(1, 2) * 0 - (M0 * M1 / 5 + M0 * M2 / 4 + M1 * M2 / 3)  # 거리 3-4=5, 3-5=4, 4-5=3, 정지 출발
K = -H
d1 = M1 * M2 / K
lamP, nuP = M2 / (M1 + M2), M1 / (M1 + M2)
muP = M0 * (M1 + M2) / (M0 + M1 + M2)
mP = M1 * M2 / (M1 + M2)
print("\nS8 피타고라스 상수(정확):")
print(f"  H = {H} = {float(H):.10f},  K = {K}")
print(f"  d1 = m1 m2 / K = {d1} = {float(d1):.12f}")
print(f"  lam = {lamP}, nu = {nuP}, lam*d1 = {lamP * d1} = {float(lamP * d1):.12f}")
print(f"  mu = {muP}, m(binary) = {mP}, M = {M0 + M1 + M2}")
print(f"  W(rho,d1) = (5/9)/(rho + {nuP * d1}) + (4/9)/(rho - {lamP * d1})")
print(f"            = (5/9)/(rho + {float(nuP * d1):.12f}) + (4/9)/(rho - {float(lamP * d1):.12f})")
ok8 = (K == Fraction(769, 60) and d1 == Fraction(1200, 769))
results.append(("S8 K=769/60, d1=1200/769", ok8))
print(f"[{'PASS' if ok8 else 'FAIL'}] S8 K=769/60, d1=1200/769")
# 초기 위치로 거리·무게중심 재확인
pos = {0: (1, 3), 1: (-2, -1), 2: (1, -1)}
dist2 = lambda i, j: (pos[i][0] - pos[j][0]) ** 2 + (pos[i][1] - pos[j][1]) ** 2
ok8b = (dist2(0, 1), dist2(0, 2), dist2(1, 2)) == (25, 16, 9)
cm = (sum([3 * 1, 4 * -2, 5 * 1]), sum([3 * 3, 4 * -1, 5 * -1]))
ok8b &= cm == (0, 0)
results.append(("S8b 초기거리 (5,4,3), 무게중심 원점", ok8b))
print(f"[{'PASS' if ok8b else 'FAIL'}] S8b 초기 거리 3-4=5, 3-5=4, 4-5=3, 무게중심 (0,0)")

print("\n요약:", sum(ok for _, ok in results), "/", len(results), "PASS")

# ---------------------------------------------------------------- 음성 대조(이 점검기가 틀린 식을 FAIL로 잡는지)
print("\n음성 대조(아래는 일부러 틀린 식이다. 모두 'FAIL'이어야 점검기가 민감하다는 뜻):")
neg = []
def neg_check(name, expr):
    z = sp.simplify(sp.together(sp.expand(expr)))
    caught = (z != 0)
    neg.append(caught)
    print(f"  [{'잡음(FAIL)' if caught else '놓침(PASS)'}] {name}")
W_sw = N_ / (rh + L * d) + L / (rh - N_ * d)        # lam, nu 를 바꾼 W
neg_check("N1 lam<->nu 바꾼 W 가 P(r=d,c=-1)과 같다", (L / r1m + N_ / r2m).subs(rr_, d) - W_sw)
neg_check("N2 운동방정식에서 lam nu 항 부호를 뒤집음",
          sum(sp.simplify(e) ** 2 for e in (rhoddot_newton
              - (-M * (lam / R13 + nu / R23) * rho - M * lam * nu * (1 / R13 - 1 / R23) * rvec))))
neg_check("N3 g 의 계수를 n P_n 으로 둠(n=3)", ser_g.coeff(a, 3) - 3 * sp.legendre(3, cser))
neg_check("N4 rho_ddot 분해에서 각운동량 항 누락(3D)", lhs - pv.dot(pdd) / rr)
print("음성 대조:", sum(neg), "/", len(neg), "잡음")
