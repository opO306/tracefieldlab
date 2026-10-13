"""조준기(비엄밀): 엄밀 실행 전에 중심 (x, y, τ)를 고정밀로 맞춘다.
판정 근거가 아니다 — Krawczyk 중심 선택에만 쓴다(ff_prescan.refine 과 같은 지위, 더 정밀할 뿐).

점 궤도: decimal 테일러 급수 적분(기본 50자리, 차수 40). 걸음 h = 0.9·ρ·tol^(1/N),
ρ = min(‖c_{N−1}‖^(−1/(N−1)), ‖c_N‖^(−1/N)) (Jorba–Zou식 수렴 반지름 추정).
g = (v1x, v1y, v2x)(τ) — ff_brake 의 eq_comps (6, 7, 8) 와 같은 성분.
J: 같은 적분기의 전진 차분. 뉴턴 z ← z − J⁻¹ g(z) (J⁻¹ 곱은 부동소수, z 갱신은 decimal).
"""
import decimal
import time
from decimal import Decimal
from fractions import Fraction

import numpy as np

PAIRS = ((0, 1), (0, 2), (1, 2))


def _dec(q):
    q = Fraction(q)
    return Decimal(q.numerator) / Decimal(q.denominator)


def _series(ms, P, V, N):
    """위치 P[6], 속도 V[6] 에서의 테일러 계수 xs[c][k], vs[c][k] (k = 0..N)."""
    xs = [[p] for p in P]
    vs = [[v] for v in V]
    pd = [([], [], [], []) for _ in PAIRS]  # dx, dy, s=|r|², w=s^(−3/2) 의 계수
    for k in range(N):
        a = [Decimal(0)] * 6
        for (i, j), (dx, dy, s, w) in zip(PAIRS, pd):
            dx.append(xs[2 * j][k] - xs[2 * i][k])
            dy.append(xs[2 * j + 1][k] - xs[2 * i + 1][k])
            s.append(sum(dx[l] * dx[k - l] + dy[l] * dy[k - l] for l in range(k + 1)))
            if k == 0:
                w.append(1 / (s[0] * s[0].sqrt()))
            else:  # s·w' = (−3/2)·s'·w  →  k·s0·w_k = Σ_{l=1..k} (−(k + l/2))·s_l·w_{k−l}
                w.append(sum(-(k + Decimal(l) / 2) * s[l] * w[k - l] for l in range(1, k + 1)) / (k * s[0]))
            qx = sum(dx[l] * w[k - l] for l in range(k + 1))
            qy = sum(dy[l] * w[k - l] for l in range(k + 1))
            a[2 * i] += ms[j] * qx
            a[2 * i + 1] += ms[j] * qy
            a[2 * j] -= ms[i] * qx
            a[2 * j + 1] -= ms[i] * qy
        for c in range(6):
            xs[c].append(vs[c][k] / (k + 1))
            vs[c].append(a[c] / (k + 1))
    return xs, vs


def _radius(cs, N):
    rho = None
    for m in (N - 1, N):
        nrm = max(abs(c[m]) for c in cs)
        if nrm > 0:
            r = 10 ** (-float(nrm.log10()) / m)
            rho = r if rho is None else min(rho, r)
    return rho


def _horner(c, h):
    v = c[-1]
    for ck in reversed(c[:-1]):
        v = v * h + ck
    return v


def _flow(ms, P, V, tau, N, tol, max_steps):
    t, steps = Decimal(0), 0
    while t < tau:
        xs, vs = _series(ms, P, V, N)
        rho = _radius(xs + vs, N)
        h = Decimal(repr(0.9 * rho * tol ** (1.0 / N))) if rho is not None else tau - t
        last = t + h >= tau
        if last:
            h = tau - t
        P = [_horner(c, h) for c in xs]
        V = [_horner(c, h) for c in vs]
        t = tau if last else t + h
        steps += 1
        if steps > max_steps:
            raise RuntimeError("too many steps (%d) at t=%s" % (steps, t))
    return P, V, steps


def _ctx(dps):
    ctx = decimal.Context(prec=dps)
    return decimal.localcontext(ctx)


def flow(masses, x, y, tau, dps=50, N=40, tol=1e-42, max_steps=200000):
    """정지 출발 상태를 τ까지 적분: (위치 6, 속도 6, 걸음 수). 입력은 Fraction/int/str 가능."""
    with _ctx(dps):
        ms = [_dec(m) for m in masses]
        P = [Decimal(-1) / 2, Decimal(0), Decimal(1) / 2, Decimal(0), _dec(x), _dec(y)]
        return _flow(ms, P, [Decimal(0)] * 6, _dec(tau), N, tol, max_steps)


def energy(masses, P, V, dps=50):
    with _ctx(dps):
        ms = [_dec(m) for m in masses]
        P = [_dec(p) if not isinstance(p, Decimal) else p for p in P]
        V = [_dec(v) if not isinstance(v, Decimal) else v for v in V]
        kin = sum(ms[i] * (V[2 * i] ** 2 + V[2 * i + 1] ** 2) for i in range(3)) / 2
        pot = sum(ms[i] * ms[j] / ((P[2 * j] - P[2 * i]) ** 2 + (P[2 * j + 1] - P[2 * i + 1]) ** 2).sqrt()
                  for i, j in PAIRS)
        return kin - pot


def aim(masses, x, y, tau, iters=6, eps="1e-15", stop=1e-32, dps=50, N=40, tol=1e-42):
    """뉴턴으로 g(x, y, τ) = 0 을 고정밀로 푼다. 반환: 가장 작은 |g|의 z와 반복 기록 (비엄밀)."""
    t0 = time.time()
    with _ctx(dps):
        ms = [_dec(m) for m in masses]
        P0 = [Decimal(-1) / 2, Decimal(0), Decimal(1) / 2, Decimal(0)]
        eps = Decimal(eps)

        def g(z):
            _P, V, st = _flow(ms, P0 + [z[0], z[1]], [Decimal(0)] * 6, z[2], N, tol, 200000)
            return V[:3], st

        z = [_dec(x), _dec(y), _dec(tau)]
        best, hist, steps = None, [], 0
        for it in range(iters + 1):
            gz, steps = g(z)
            gn = float(max(abs(v) for v in gz))
            hist.append(gn)
            if best is None or gn < best[0]:
                best = (gn, list(z))
            if gn < stop or it == iters or (len(hist) >= 2 and gn >= hist[-2]):
                break
            J = np.empty((3, 3))
            for k in range(3):
                zk = list(z)
                zk[k] += eps
                gk, _ = g(zk)
                J[:, k] = [float((a - b) / eps) for a, b in zip(gk, gz)]
            dz = np.linalg.solve(J, [float(v) for v in gz])
            z = [zi - Decimal(repr(float(d))) for zi, d in zip(z, dz)]
        zb = best[1]
        return {"x": str(zb[0]), "y": str(zb[1]), "tau": str(zb[2]), "hist": hist, "g_best": best[0], "steps": steps,
                "sec": round(time.time() - t0, 1), "dps": dps, "order": N, "tol": tol}
