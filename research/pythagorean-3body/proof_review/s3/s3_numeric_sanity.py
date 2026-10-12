"""(비엄밀) 정칙화 계 (LC)를 피타고라스 초기조건으로 s-시간 적분하고, 뉴턴 방정식 t-시간 적분과 비교.
수식 세팅(부호, Jacobi 좌표, 초기값 사상)의 실수를 잡기 위한 점검이며 증명의 일부가 아니다."""
import numpy as np
from scipy.integrate import solve_ivp

m0, m1, m2 = 3.0, 4.0, 5.0            # m0 = 탈출 후보, m1 <= m2
x0, x1, x2 = np.array([1., 3.]), np.array([-2., -1.]), np.array([1., -1.])
k = m1 + m2; M = m0 + k; lam = m2 / k; nu = m1 / k

def newton(t, y):
    X = y[:6].reshape(3, 2); ms = (m0, m1, m2); a = np.zeros((3, 2))
    for i in range(3):
        for j in range(3):
            if i != j:
                d = X[j] - X[i]; a[i] += ms[j] * d / np.linalg.norm(d) ** 3
    return np.concatenate([y[6:], a.ravel()])

def jac(X):
    r = X[2] - X[1]; g = (m1 * X[1] + m2 * X[2]) / k; return r, X[0] - g

def fPhi(rho, z):
    d1 = rho - nu * z; d2 = rho + lam * z; R1 = np.linalg.norm(d1); R2 = np.linalg.norm(d2)
    f = m0 * (d1 / R1**3 - d2 / R2**3)
    Phi = -M * (lam / R1**3 + nu / R2**3) * rho + lam * nu * M * (1 / R1**3 - 1 / R2**3) * z
    return f, Phi

def lc(s, Y):
    t, u, w, h, rho, sig = Y[0], complex(Y[1], Y[2]), complex(Y[3], Y[4]), Y[5], Y[6:8], Y[8:10]
    z = u * u; f, Phi = fPhi(rho, np.array([z.real, z.imag])); fc = complex(*f); U2 = abs(u) ** 2
    dw = h / 2 * u + U2 / 2 * u.conjugate() * fc
    dh = 2 * ((u * w).conjugate() * fc).real
    return [U2, w.real, w.imag, dw.real, dw.imag, dh, *(U2 * sig), *(U2 * Phi)]

T = 30.0
y0 = np.concatenate([np.array([x0, x1, x2]).ravel(), np.zeros(6)])
nt = solve_ivp(newton, (0, T), y0, method='DOP853', rtol=1e-13, atol=1e-15)
XT = nt.y[:6, -1].reshape(3, 2); rN, rhoN = jac(XT)
dists = []
sol = solve_ivp(newton, (0, T), y0, method='DOP853', rtol=1e-13, atol=1e-15, dense_output=True)
ts = np.linspace(0, T, 300001); Xs = sol.sol(ts)[:6].reshape(3, 2, -1)
r12 = np.linalg.norm(Xs[2] - Xs[1], axis=0); print('min r(4-5) on [0,T] (Newton, sampled):', r12.min())

r0, rho0 = jac(np.array([x0, x1, x2]))
u0 = np.sqrt(complex(*r0)); w0 = u0.conjugate() * 0 / 2; h0 = -k / np.linalg.norm(r0)
Y0 = [0.0, u0.real, u0.imag, w0.real, w0.imag, h0, *rho0, 0.0, 0.0]
ev = lambda s, Y: Y[0] - T; ev.terminal = True
ls = solve_ivp(lc, (0, 1e4), Y0, method='DOP853', rtol=1e-13, atol=1e-15, events=ev)
Yf = ls.y_events[0][0]; uf = complex(Yf[1], Yf[2]); zf = uf * uf
print('t at end (LC):', Yf[0])
print('|r_LC - r_Newton|   =', np.linalg.norm(np.array([zf.real, zf.imag]) - rN))
print('|rho_LC - rho_Newton| =', np.linalg.norm(Yf[6:8] - rhoN))
I = 2 * (Yf[3]**2 + Yf[4]**2) - Yf[5] * abs(uf)**2
print('invariant I - k at end:', I - k)
