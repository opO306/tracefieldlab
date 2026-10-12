import os, sys
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ff_common import dyad, DI
from ff_krawczyk import krawczyk, imat_mul, imat_vec

dyad.set_prec(128)


def test_sqrt2_pass_and_contains_root():
    X = [DI(F(181, 128), F(182, 128))]          # [1.4140625, 1.421875] ∋ √2
    m = [F(363, 256)]
    gm = [DI(m[0]) * DI(m[0]) - 2]
    Dg = [[X[0] * 2]]
    r = krawczyk(gm, Dg, X, m)
    assert r["verdict"] == "PASS"
    k = r["K"][0]
    assert k.lo_q * k.lo_q < 2 < k.hi_q * k.hi_q


def test_box_without_root_is_refused():
    X = [DI(F(3, 2), F(13, 8))]                  # no √2 inside
    m = [F(25, 16)]
    gm = [DI(m[0]) * DI(m[0]) - 2]
    Dg = [[X[0] * 2]]
    assert krawczyk(gm, Dg, X, m)["verdict"] == "REFUSED"


def test_linear_2d_exact_root():
    # A z = b, A=[[2,1],[1,3]], b=[3,5] -> z = (4/5, 7/5)
    A = [[DI(2), DI(1)], [DI(1), DI(3)]]
    X = [DI(F(1, 2), F(1)), DI(F(1), F(2))]
    m = [F(3, 4), F(3, 2)]
    gm = [a - b for a, b in zip(imat_vec(A, [DI(v) for v in m]), [DI(3), DI(5)])]
    r = krawczyk(gm, A, X, m)
    assert r["verdict"] == "PASS"
    assert r["K"][0].contains(F(4, 5)) and r["K"][1].contains(F(7, 5))


def test_imat_mul_shapes():
    A = [[DI(1), DI(2)]]                         # 1x2
    B = [[DI(3)], [DI(4)]]                       # 2x1
    assert imat_mul(A, B)[0][0].contains(11)
    assert imat_mul([[DI(1)], [DI(2)]], [[]]) == [[], []]   # n x 0 allowed


def test_singular_derivative_is_refused_not_crash():
    X = [DI(F(1), F(2))]
    m = [F(3, 2)]
    r = krawczyk([DI(F(1, 4))], [[DI(0)]], X, m)
    assert r["verdict"] == "REFUSED" and "singular" in r["reason"]
