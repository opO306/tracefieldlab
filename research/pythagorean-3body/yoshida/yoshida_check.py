"""Yoshida(1972) 계열 탈출 충분조건의 독립 수치 점검 (비엄밀, 참고용).

이 스크립트는 저장소의 다른 코드에 의존하지 않는다.
1) 에너지 E = -769/60 을 Fraction으로 정확 계산
2) 보조정리(퍼텐셜 상계/가속도 하계)를 격자 전수 탐색으로 점검
3) DOP853 (rtol 1e-13, atol 1e-15) 로 t=0..70 적분, 공차를 바꿔 교차 비교
4) 독립 유도한 조건(= Yoshida Corollary / Theorem 1 형태)이 처음 성립하는 시각과 여유

실행: python -I yoshida_check.py   (단일 코어, 수십 초 ~ 수 분)
"""
from fractions import Fraction as Fr
import math
import numpy as np
from scipy.integrate import solve_ivp

# ---- 문제 정의 -------------------------------------------------------------
# m0 = 탈출체(질량 3), m1, m2 = 쌍성(질량 4, 5)  -- Yoshida 기호
m0, m1, m2 = 3.0, 4.0, 5.0
P0 = [(1, 3), (-2, -1), (1, -1)]  # 질량 3, 4, 5 의 초기 위치
masses = np.array([m0, m1, m2])
M = m0 + m1 + m2
lam = m2 / (m1 + m2)   # G-m1 거리 = lam*r
nu = m1 / (m1 + m2)    # G-m2 거리 = nu*r
mu = m0 * (m1 + m2) / M
mred = m1 * m2 / (m1 + m2)


def exact_energy():
    ms = [Fr(3), Fr(4), Fr(5)]
    U = Fr(0)
    for i in range(3):
        for j in range(i + 1, 3):
            dx = P0[i][0] - P0[j][0]
            dy = P0[i][1] - P0[j][1]
            d2 = dx * dx + dy * dy
            s = math.isqrt(d2)
            assert s * s == d2, "거리가 정수가 아님"
            U -= ms[i] * ms[j] / s
    return U  # 정지 출발이므로 E = U


E_exact = exact_energy()
K = -float(E_exact)
d1 = m1 * m2 / K


# ---- 보조정리 격자 점검 ----------------------------------------------------
def lemma_grid_check():
    """rho 고정, 0<=r<=d, 0<=theta<=pi 에서
    S1 = lam/r1 + nu/r2 의 sup 이 lam/(rho+nu d)+nu/(rho-lam d) (lam>=nu) 인지,
    F2 = -M(lam/r1^3+nu/r2^3) rho + lam nu M (1/r1^3-1/r2^3) r cos(theta) 의 inf 가
    -M{lam/(rho+nu d)^2 + nu/(rho-lam d)^2} 이상인지 격자로 확인."""
    worst1, worst2 = -1e9, -1e9
    for a in [0.6, 0.8, 1.0, 1.5, 3.0]:
        for d in [0.3, 1.0, 1.56]:
            rho = a * d
            r = np.linspace(0, d, 401)[:, None]
            th = np.linspace(0, math.pi, 721)[None, :]
            r1 = np.sqrt(rho**2 + nu**2 * r**2 - 2 * nu * r * rho * np.cos(th))
            r2 = np.sqrt(rho**2 + lam**2 * r**2 + 2 * lam * r * rho * np.cos(th))
            S1 = lam / r1 + nu / r2
            bound1 = lam / (rho + nu * d) + nu / (rho - lam * d)
            worst1 = max(worst1, (S1.max() - bound1) / bound1)
            F2 = -M * (lam / r1**3 + nu / r2**3) * rho + lam * nu * M * (1 / r1**3 - 1 / r2**3) * r * np.cos(th)
            bound2 = -M * (lam / (rho + nu * d) ** 2 + nu / (rho - lam * d) ** 2)
            worst2 = max(worst2, (bound2 - F2.min()) / abs(bound2))
    return worst1, worst2  # <= 0 (반올림 수준) 이면 상계/하계가 격자에서 성립


# ---- 운동방정식 ------------------------------------------------------------
def rhs(t, y):
    x = y[:6].reshape(3, 2)
    v = y[6:].reshape(3, 2)
    a = np.zeros((3, 2))
    for i in range(3):
        for j in range(i + 1, 3):
            dx = x[j] - x[i]
            r3 = (dx @ dx) ** 1.5
            a[i] += masses[j] * dx / r3
            a[j] -= masses[i] * dx / r3
    return np.concatenate([v.ravel(), a.ravel()])


def energy(y):
    x = y[:6].reshape(3, 2)
    v = y[6:].reshape(3, 2)
    T = 0.5 * sum(masses[i] * v[i] @ v[i] for i in range(3))
    U = 0.0
    for i in range(3):
        for j in range(i + 1, 3):
            U -= masses[i] * masses[j] / np.linalg.norm(x[i] - x[j])
    return T + U


def jacobi(y):
    x = y[:6].reshape(3, 2)
    v = y[6:].reshape(3, 2)
    rv = x[2] - x[1]; rd = v[2] - v[1]                    # r = m1->m2
    G = (m1 * x[1] + m2 * x[2]) / (m1 + m2)
    Gd = (m1 * v[1] + m2 * v[2]) / (m1 + m2)
    pv = x[0] - G; pd = v[0] - Gd                          # rho = G->m0
    rho = np.linalg.norm(pv); r = np.linalg.norm(rv)
    rhodot = pv @ pd / rho
    rdot = rv @ rd / r
    return rho, rhodot, r, rdot


def W(rho, d):
    """sup_{r<=d} (lam/r1 + nu/r2), lam>=nu (m1<=m2)."""
    return lam / (rho + nu * d) + nu / (rho - lam * d)


def corollary_margin(y):
    """Yoshida Corollary (d=d1=m1 m2/K) 형태. 반환: (성립여부, 여유들)."""
    rho, rhodot, r, _ = jacobi(y)
    if rho <= lam * d1:
        return False, None
    need = math.sqrt(2 * M * W(rho, d1))
    a_room = rho / max(d1, r)  # 이 값이 lam 보다 크면 lam<a<=a_room 인 a 존재
    ok = (a_room > lam) and (rhodot > need)
    return ok, dict(rho=rho, rhodot=rhodot, r=r, need=need, a_room=a_room,
                    vel_margin=rhodot - need, a_margin=a_room - lam)


def theorem1_best(y, dgrid):
    """Theorem 1 (임의 d): rho>=a d>a r, rhodot>sqrt(2M W), rhodot^2-2M W > (2H/mu+2 m1 m2/(mu d))."""
    rho, rhodot, r, _ = jacobi(y)
    best = None
    for d in dgrid:
        if d <= r or rho <= lam * d or rho <= lam * r:
            continue
        w = 2 * M * W(rho, d)
        c14 = rhodot - math.sqrt(w)
        c15 = (rhodot**2 - w) - (-2 * K / mu + 2 * m1 * m2 / (mu * d))
        m = min(c14, c15)
        if best is None or m > best[0]:
            best = (m, d, c14, c15)
    return best


def integrate(rtol, atol, tmax=70.0):
    y0 = np.array([1, 3, -2, -1, 1, -1, 0, 0, 0, 0, 0, 0], float)
    return solve_ivp(rhs, (0, tmax), y0, method="DOP853", rtol=rtol, atol=atol,
                     dense_output=True)


def main():
    print("E exact =", E_exact, "=", float(E_exact))
    print("lam, nu, M, mu, K, d1 =", lam, nu, M, mu, K, d1)
    w1, w2 = lemma_grid_check()
    print("lemma grid worst rel. excess: S1 sup %.2e, F2 inf %.2e (<=~0 이면 OK)" % (w1, w2))

    sols = {}
    for rt in [1e-12, 1e-13, 1e-14]:
        s = integrate(rt, rt * 1e-2)
        sols[rt] = s
        dE = max(abs(energy(s.y[:, k]) - float(E_exact)) for k in range(s.y.shape[1]))
        print(f"rtol={rt:g}: steps={s.t.size}, status={s.status}, max|dE|={dE:.2e}")

    s = sols[1e-13]
    # 최근접 조우 점검
    tt = np.linspace(15.5, 16.2, 70001)
    Y = s.sol(tt)
    d12 = np.hypot(Y[0] - Y[4], Y[1] - Y[5])  # 질량 3-5
    d01 = np.hypot(Y[0] - Y[2], Y[1] - Y[3])  # 3-4
    d23 = np.hypot(Y[2] - Y[4], Y[3] - Y[5])  # 4-5
    for name, dd in [("3-4", d01), ("3-5", d12), ("4-5", d23)]:
        k = dd.argmin()
        print(f"min dist {name} in [15.5,16.2]: {dd[k]:.3e} at t={tt[k]:.5f}")

    # 공차 간 상태 차이
    for t in [16.0, 30.0, 50.0, 59.72, 59.99, 70.0]:
        a = sols[1e-13].sol(t); b = sols[1e-14].sol(t); c = sols[1e-12].sol(t)
        print(f"t={t}: |y(1e-13)-y(1e-14)|={np.max(abs(a-b)):.2e}  |y(1e-12)-y(1e-14)|={np.max(abs(c-b)):.2e}")

    # Yoshida 논문의 Yuasa 데이터와 대조 (질량중심 좌표; 1=질량3,2=질량4,3=질량5 로 가정)
    for t in [59.72, 59.99]:
        y = sols[1e-14].sol(t)
        print(f"t={t} state (x,y per body):", np.round(y[:6], 5), "v:", np.round(y[6:], 5))
        print("   jacobi rho,rhodot,r =", [round(v, 5) for v in jacobi(y)[:3]])

    # 처음 성립 시각 스캔
    grid = np.arange(0.0, 70.0, 0.001)
    for rt in [1e-12, 1e-13, 1e-14]:
        sol = sols[rt].sol
        first = None
        for t in grid:
            ok, info = corollary_margin(sol(t))
            if ok:
                first = (t, info); break
        if first:
            t0, info = first
            # 이분법으로 정밀화
            lo, hi = t0 - 0.001, t0
            for _ in range(40):
                mid = 0.5 * (lo + hi)
                if corollary_margin(sol(mid))[0]: hi = mid
                else: lo = mid
            print(f"[Corollary d=d1] rtol={rt:g}: first t = {hi:.6f}; info at t0={t0:.3f}:",
                  {k: round(v, 5) for k, v in info.items()})
        # 이후 유지 여부
        bad = [t for t in np.arange(first[0][0] if False else (first[0] if first else 0), 70, 0.01)
               if not corollary_margin(sol(t))[0]] if first else []
        print(f"   corollary 미성립 시각 수(이후 0.01 간격): {len(bad)}" + (f", 마지막 {bad[-1]:.2f}" if bad else ""))

    sol = sols[1e-14].sol
    for t in [55, 58, 59, 59.5, 59.72, 60, 61, 62, 65, 70]:
        ok, info = corollary_margin(sol(t))
        print(f"t={t}: corollary ok={ok}", {k: round(v, 4) for k, v in (info or {}).items()})

    # Theorem 1 (d 최적화)
    dgrid = np.linspace(0.05, 5.0, 1000)
    firstT1 = None
    for t in grid:
        b = theorem1_best(sol(t), dgrid)
        if b and b[0] > 0:
            firstT1 = (t, b); break
    print("[Theorem 1, d 최적화(격자)] first t =", firstT1)


if __name__ == "__main__":
    main()
