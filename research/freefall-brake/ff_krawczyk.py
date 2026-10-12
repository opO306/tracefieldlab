"""Krawczyk 판정: K(X) = m − C·g(m) + (I − C·[Dg](X))·(X − m).

K(X) ⊆ X 이면 X 안에 g의 영점이 있다(존재; Schichl–Neumaier 2004 Thm 2.2, 출처는 C1_RESULT_KO.md).
여기서는 엄격한 포함 K(X) ⊂ int X 를 PASS 조건으로 쓴다. 유일성은 이 모듈이 아니라
ff_check 의 ‖I − C·[Dg]‖∞ < 1 정확 계산으로 따로 판정한다.
C는 근사 역행렬이면 충분하다(정리는 임의의 C에 대해 성립).
구간 연산은 dyad의 바깥 반올림이라 K는 참 Krawczyk 상자를 포함한다.
"""
from fractions import Fraction

import numpy as np

from ff_common import DI


def _sum(xs):
    acc = xs[0]
    for x in xs[1:]:
        acc = acc + x
    return acc


def imat_mul(A, B):
    k = len(B)
    cols = len(B[0]) if k else 0
    return [[_sum([A[i][l] * B[l][j] for l in range(k)]) for j in range(cols)] for i in range(len(A))]


def imat_vec(A, v):
    return [_sum([A[i][l] * v[l] for l in range(len(v))]) for i in range(len(A))]


def _mid(x):
    return (x.lo_q + x.hi_q) / 2


def approx_inverse(Dg):
    M = np.array([[float(_mid(e)) for e in row] for row in Dg], dtype=float)
    Ci = np.linalg.inv(M)
    return [[Fraction(float(x)) for x in row] for row in Ci]   # float는 정확한 이진 분수


def krawczyk(gm, Dg, X, m):
    q = len(X)
    for x, mi in zip(X, m):
        if not (x.lo_q <= mi <= x.hi_q):
            raise ValueError("centre m must lie in X")
    try:
        C = approx_inverse(Dg)
    except (np.linalg.LinAlgError, ValueError, OverflowError) as e:
        return {"verdict": "REFUSED", "K": None, "C": None,
                "reason": "singular or non-finite mid([Dg]): %s" % e}
    CD = [[DI(c) for c in row] for row in C]          # prec ≥ 53이면 정확 (작으면 바깥 반올림, 여전히 건전)
    CDg = imat_mul(CD, Dg)
    I_CDg = [[(1 if i == j else 0) - CDg[i][j] for j in range(q)] for i in range(q)]
    Xm = [x - DI(mi) for x, mi in zip(X, m)]
    Cg = imat_vec(CD, gm)
    corr = imat_vec(I_CDg, Xm)
    K = [DI(mi) - a + b for mi, a, b in zip(m, Cg, corr)]
    inside = all(k.lo_q > x.lo_q and k.hi_q < x.hi_q for k, x in zip(K, X))
    return {"verdict": "PASS" if inside else "REFUSED", "K": K, "C": C,
            "reason": "" if inside else "K(X) not inside int X"}


def q2s(q):
    q = Fraction(q)
    return "%d/%d" % (q.numerator, q.denominator)


def di2s(x):
    return [q2s(x.lo_q), q2s(x.hi_q)]
