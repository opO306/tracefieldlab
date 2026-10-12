"""S3 정칙화 계에 대한 기호 검산 (sympy).

정칙화 계 (s = 가상 시간, ' = d/ds, 복소수 = 평면 벡터, z = u^2 = 상대위치 r 벡터):
  t' = |u|^2,  u' = w,  w' = (h/2) u + (1/2)|u|^2 conj(u) f,  h' = Re(conj(2uw) f),
  rho' = |u|^2 sigma,  sigma' = |u|^2 Phi
  f   = m0[(rho - nu z)/r1^3 - (rho + lam z)/r2^3],  r1 = |rho - nu z|, r2 = |rho + lam z|
  Phi = -M(lam/r1^3 + nu/r2^3) rho + lam nu M (1/r1^3 - 1/r2^3) z
r1, r2 는 기호 R1, R2 로 두고(미분은 연쇄법칙), 마지막에 R1^2 = |rho - nu z|^2, R2^2 = |rho + lam z|^2 로
나머지 연산(rem)해서 0 인지 본다. 출력이 전부 0 이면 항등식 성립.
"""
import sympy as sp

m0, m1, m2 = sp.symbols('m0 m1 m2', positive=True)
k = m1 + m2; M = m0 + k; lam = m2 / k; nu = m1 / k; m = m1 * m2 / k; mu = m0 * k / M
u1, u2, w1, w2, h, p1, p2, s1, s2, R1, R2 = sp.symbols('u1 u2 w1 w2 h p1 p2 s1 s2 R1 R2', real=True)

def cm(a, b): return (a[0]*b[0] - a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def cj(a): return (a[0], -a[1])
def sc(c, a): return (c*a[0], c*a[1])
def add(*vs): return (sum(v[0] for v in vs), sum(v[1] for v in vs))
def dot(a, b): return a[0]*b[0] + a[1]*b[1]          # = Re(conj(a) b)
def n2(a): return dot(a, a)

u = (u1, u2); w = (w1, w2); rho = (p1, p2); sig = (s1, s2)
z = cm(u, u); U2 = n2(u)
d1 = add(rho, sc(-nu, z)); d2 = add(rho, sc(lam, z))
f = add(sc(m0 / R1**3, d1), sc(-m0 / R2**3, d2))
Phi = add(sc(-M*(lam/R1**3 + nu/R2**3), rho), sc(lam*nu*M*(1/R1**3 - 1/R2**3), z))
P = lam / R1 + nu / R2

du = w
dw = add(sc(h/2, u), sc(U2/2, cm(cj(u), f)))
dh = 2 * dot(cm(u, w), f)                          # Re(conj(2uw) f)
drho = sc(U2, sig); dsig = sc(U2, Phi)
dz = sc(2, cm(u, w))
dR1 = dot(d1, add(drho, sc(-nu, dz))) / R1
dR2 = dot(d2, add(drho, sc(lam, dz))) / R2
state = [u1, u2, w1, w2, h, p1, p2, s1, s2, R1, R2]
field = [du[0], du[1], dw[0], dw[1], dh, drho[0], drho[1], dsig[0], dsig[1], dR1, dR2]
def Lie(F): return sum(sp.diff(F, x) * v for x, v in zip(state, field))

def is_zero(expr, extra=None):
    num = sp.numer(sp.together(sp.expand(expr)))
    if extra:
        num = sp.numer(sp.together(sp.expand(num.subs(extra))))
    num = sp.expand(num)
    num = sp.rem(sp.Poly(num, R1), sp.Poly(R1**2 - n2(d1), R1)).as_expr()
    num = sp.rem(sp.Poly(sp.expand(num), R2), sp.Poly(R2**2 - n2(d2), R2)).as_expr()
    return sp.expand(num) == 0

res = {}
res['a: 2|w|^2 - h|u|^2 conserved'] = is_zero(Lie(2*n2(w) - h*U2))
res['b: H = m h + mu|sigma|^2/2 - mu M P conserved'] = is_zero(Lie(m*h + mu*n2(sig)/2 - mu*M*P))
zdot = sc(1/U2, dz)                                # dz/dt = (dz/ds)/t'
zdd = [Lie(c) / U2 for c in zdot]
target = add(sc(-k / U2**3, z), f)                 # -k z/|z|^3 + f, |z| = |u|^2
I_sub = {h: (2*n2(w) - k) / U2}                   # 불변량 = k
res['c: z_tt = -k z/|z|^3 + f (x)'] = is_zero(zdd[0] - target[0], I_sub)
res['c: z_tt = -k z/|z|^3 + f (y)'] = is_zero(zdd[1] - target[1], I_sub)
res['d: h = |z_t|^2/2 - k/|z|'] = is_zero(n2(zdot)/2 - k/U2 - h, I_sub)
res['d2: rho_t = sigma'] = all(sp.simplify(drho[i]/U2 - sig[i]) == 0 for i in range(2))

# (e) 뉴턴 방정식 대조: 절대좌표에서 x0-x2 = rho - nu r, x0-x1 = rho + lam r, x2-x1 = r
X = sp.symbols('x0a x0b x1a x1b x2a x2b', real=True); Rr = sp.symbols('Rr', positive=True)
x0, x1, x2 = (X[0], X[1]), (X[2], X[3]), (X[4], X[5])
rN = add(x2, sc(-1, x1)); g = sc(1/k, add(sc(m1, x1), sc(m2, x2))); rhoN = add(x0, sc(-1, g))
geo = [sp.simplify(c) for c in add(rhoN, sc(-nu, rN), sc(-1, add(x0, sc(-1, x2))))] + \
      [sp.simplify(c) for c in add(rhoN, sc(lam, rN), sc(-1, add(x0, sc(-1, x1))))]
res['e0: Jacobi geometry'] = all(c == 0 for c in geo)
# 뉴턴 가속도 (거리 기호: |x0-x2| = R1, |x0-x1| = R2, |x2-x1| = Rr)
a0 = add(sc(m1 / R2**3, add(x1, sc(-1, x0))), sc(m2 / R1**3, add(x2, sc(-1, x0))))
a1 = add(sc(m0 / R2**3, add(x0, sc(-1, x1))), sc(m2 / Rr**3, add(x2, sc(-1, x1))))
a2 = add(sc(m0 / R1**3, add(x0, sc(-1, x2))), sc(m1 / Rr**3, add(x1, sc(-1, x2))))
rdd = add(a2, sc(-1, a1)); gdd = sc(1/k, add(sc(m1, a1), sc(m2, a2))); rhodd = add(a0, sc(-1, gdd))
fN = add(sc(m0 / R1**3, add(rhoN, sc(-nu, rN))), sc(-m0 / R2**3, add(rhoN, sc(lam, rN))))
kep = add(sc(-k / Rr**3, rN), fN)
PhiN = add(sc(-M*(lam/R1**3 + nu/R2**3), rhoN), sc(lam*nu*M*(1/R1**3 - 1/R2**3), rN))
res['e1: r_tt (Newton) = -k r/r^3 + f'] = all(sp.simplify(rdd[i] - kep[i]) == 0 for i in range(2))
res['e2: rho_tt (Newton) = Phi'] = all(sp.simplify(rhodd[i] - PhiN[i]) == 0 for i in range(2))
# (f) 에너지 식 (E) 대조: 뉴턴 에너지 = m h + mu|rho_t|^2/2 - mu M P (h 는 Kepler 에너지)
V = sp.symbols('v0a v0b v1a v1b v2a v2b', real=True)
v0, v1, v2 = (V[0], V[1]), (V[2], V[3]), (V[4], V[5])
mom = add(sc(m0, v0), sc(m1, v1), sc(m2, v2))
EN = (m0*n2(v0) + m1*n2(v1) + m2*n2(v2))/2 - m0*m1/R2 - m0*m2/R1 - m1*m2/Rr
rv = add(v2, sc(-1, v1)); rhov = add(v0, sc(-1, sc(1/k, add(sc(m1, v1), sc(m2, v2)))))
EJ = m*(n2(rv)/2 - k/Rr) + mu*n2(rhov)/2 - mu*M*(lam/R1 + nu/R2)
diffE = sp.expand(EN - EJ)
# 무게중심계(총운동량 0)에서 비교: v0 = -(m1 v1 + m2 v2)/m0
diffE0 = sp.simplify(diffE.subs({V[0]: -(m1*V[2] + m2*V[4])/m0, V[1]: -(m1*V[3] + m2*V[5])/m0}))
res['f: energy (E) in barycentric frame'] = diffE0 == 0
for kk, v in res.items():
    print('%-45s %s' % (kk, 'OK' if v else 'FAIL'))
print('ALL OK' if all(res.values()) else 'SOME FAIL')
