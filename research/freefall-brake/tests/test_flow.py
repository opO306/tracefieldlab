import os, sys
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mpmath import mp
from ff_common import dyad, DI, ts, ls
from ff_flow import integrate, flow_on_tau_box

dyad.set_prec(128)
mp.dps = 60


def br(x, eps=F(1, 10 ** 50)):
    q = F(mp.nstr(x, 58))
    return q - eps, q + eps


def contains_real(d, x):
    lo, hi = br(x)
    return d.lo_q <= lo and hi <= d.hi_q


def test_point_integration_lands_exactly_and_encloses_cos_sin():
    sysm = ts.make_harmonic()
    st, _, n = integrate(sysm, ls.lohner_init(sysm, [F(1), F(0)]), F(1), F(1, 8))
    assert st.t == 1 and n >= 8
    box = ls.enclosure_box(st)
    assert contains_real(box[0], mp.cos(1)) and contains_real(box[1], -mp.sin(1))
    assert max(b.width_q() for b in box) < F(1, 10 ** 30)


def test_tau_box_derivative_harmonic():
    sysm = ts.make_harmonic()
    r = F(1, 2 ** 30)
    st0 = ls.LohnerState([F(1), F(0)], [[F(1), F(0)], [F(0), F(1)]], [DI(-r, r), DI(0)], F(0))
    S = [[DI(1)], [DI(0)]]                       # parameter = x0
    P, D, Fv = flow_on_tau_box(sysm, st0, S, F(1), F(1) + r, F(1, 8))
    assert contains_real(D[0][0], mp.cos(1)) and contains_real(D[1][0], -mp.sin(1))
    assert contains_real(Fv[1], -mp.cos(1))     # v' = -x, x(1)=cos 1 (x0=1)
    assert D[0][0].width_q() < F(1, 10 ** 6)
