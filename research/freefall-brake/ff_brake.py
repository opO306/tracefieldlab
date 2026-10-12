"""삼체 정지(free-fall) 문제: 별1 (−1/2,0), 별2 (1/2,0), 별3 (x,y), 모든 속도 0, G=1 (Li & Liao 2019 2절).

미지수 (x, y, τ); 방정식 v_E(τ) = 0 (속도 3성분). 나머지 3성분은 보존량 P=0(2식), L=0(1식)에서
A_R v_R = −A_E v_E 이므로, A_R(위치 구간)의 행렬식이 0을 포함하지 않으면 v_E=0 ⇒ v_R=0.
"""
from fractions import Fraction

from ff_common import DI, ts
from ff_problem import Problem


def _row_entries(masses, P):
    """보존식 3줄(Px, Py, L)의 속도 성분 6..11 계수 (구간)."""
    rows = [[DI(0)] * 12 for _ in range(3)]
    for i, mi in enumerate(masses):
        M = DI(Fraction(mi))
        x, y = P[2 * i], P[2 * i + 1]
        vx, vy = 6 + 2 * i, 7 + 2 * i
        rows[0][vx] = M
        rows[1][vy] = M
        rows[2][vx] = -(M * y)
        rows[2][vy] = M * x
    return rows


def _det3(a):
    return (a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
            - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
            + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]))


def rest_velocity_check(masses, eq_comps, P):
    rest = [c for c in range(6, 12) if c not in eq_comps]
    rows = _row_entries(masses, P)
    A_R = [[rows[k][c] for c in rest] for k in range(3)]
    d = _det3(A_R)
    ok = d.lo_q > 0 or d.hi_q < 0
    return ok, "det A_R in [%.3e, %.3e]" % (float(d.lo_q), float(d.hi_q))


def brake_problem(name, masses, h_max=Fraction(1, 8), eq_comps=(6, 7, 8), system=None):
    masses = [Fraction(m) for m in masses]
    sysm = system if system is not None else ts.make_nbody(masses)
    base = [Fraction(-1, 2), Fraction(0), Fraction(1, 2), Fraction(0), Fraction(0), Fraction(0)] + [Fraction(0)] * 6
    return Problem(name, sysm, base, [4, 5], list(eq_comps), [0, 0, 0], h_max, tol_control=True,
                   post_check=lambda P: rest_velocity_check(masses, eq_comps, P))
