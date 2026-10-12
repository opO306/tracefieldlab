"""부동소수 사전 조사와 중심 다듬기. 비엄밀 — Krawczyk 중심 선택에만 쓰고 판정 근거가 아니다."""
import math

import numpy as np
from scipy.integrate import solve_ivp

PAIRS = [(0, 1), (0, 2), (1, 2)]


def _rhs(m):
    def f(t, u):
        p = u[:6].reshape(3, 2)
        a = np.zeros((3, 2))
        for i, j in PAIRS:
            r = p[j] - p[i]
            k = r / np.linalg.norm(r) ** 3
            a[i] += m[j] * k
            a[j] -= m[i] * k
        return np.concatenate([u[6:], a.ravel()])
    return f


def _solve(masses, x, y, tau, dense=False):
    u0 = np.array([-0.5, 0, 0.5, 0, x, y] + [0] * 6, float)
    return solve_ivp(_rhs(np.array(masses, float)), (0, tau), u0, method="DOP853",
                     rtol=1e-13, atol=1e-15, dense_output=dense)


def prescan(masses, x, y, T):
    sol = _solve(masses, x, y, T / 2, dense=True)
    ts = np.linspace(0, T / 2, 20001)
    P = sol.sol(ts)[:6].T.reshape(-1, 3, 2)
    dmin = min(float(np.min(np.linalg.norm(P[:, i] - P[:, j], axis=1))) for i, j in PAIRS)
    return {"min_dist": dmin, "vmax_half": float(np.max(np.abs(sol.y[6:, -1]))),
            "x2_minus_x3": float(sol.y[2, -1] - sol.y[4, -1]), "nfev": int(sol.nfev)}


def _resid(masses, z):
    s = _solve(masses, z[0], z[1], z[2])
    return s.y[6:9, -1]


def refine(masses, x, y, tau, iters=4, eps=1e-7):
    z = np.array([x, y, tau], float)
    for _ in range(iters):
        g = _resid(masses, z)
        J = np.empty((3, 3))
        for k in range(3):
            dz = np.zeros(3)
            dz[k] = eps
            J[:, k] = (_resid(masses, z + dz) - _resid(masses, z - dz)) / (2 * eps)
        z = z - np.linalg.solve(J, g)
    return {"x": float(z[0]), "y": float(z[1]), "tau": float(z[2]),
            "resid": float(np.max(np.abs(_resid(masses, z))))}


def choose_times(masses, x, y, tau, k, quantum=2.0 ** -16):
    """근접도 가중 ∫ r_min^(-3/2) dt 를 k등분하는 시각 k−1개 (quantum 배수로 내림). 비엄밀 — 역 위치 선택용."""
    sol = _solve(masses, x, y, tau, dense=True)
    ts = np.linspace(0, tau, 20001)
    P = sol.sol(ts)[:6].T.reshape(-1, 3, 2)
    d = np.min([np.linalg.norm(P[:, i] - P[:, j], axis=1) for i, j in PAIRS], axis=0)
    w = d ** -1.5
    cum = np.concatenate([[0.0], np.cumsum((w[1:] + w[:-1]) / 2 * np.diff(ts))])
    return [math.floor(float(np.interp(i / k * cum[-1], cum, ts)) / quantum) * quantum for i in range(1, k)]
