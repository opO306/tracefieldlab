import math, os, sys
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mpmath import mp
from ff_common import dyad, ts, ls
from ff_flow import integrate, step_size, FlowRefused
from ff_problem import Problem, run_problem

dyad.set_prec(128)
mp.dps = 60


def contains_real(pair, x, eps=F(1, 10 ** 50)):
    q = F(mp.nstr(x, 58))
    return F(pair[0]) <= q - eps and q + eps <= F(pair[1])


def test_step_size_rule_matches_r4():
    assert step_size(F(1, 2 ** 112), F(0)) == F(1, 8)
    h = step_size(F(1, 2 ** 112), F(10))
    assert h.numerator == 1 and h.denominator & (h.denominator - 1) == 0 and h < F(1, 8)


def test_kepler_passes_with_tol_control_and_large_hmax():
    # C1에서 고정 h=1/16은 REFUSED였다. 보폭 제어를 켜면 h_max=1/8이어도 PASS.
    p = Problem("kep_tol", ts.make_kepler(1), [F(3, 2), F(0), F(0), F(0)], [3], [0, 1],
                [F(-1, 2), F(0)], F(1, 8), tol_control=True)
    res = run_problem(p, [F(1 / math.sqrt(3)), F(math.pi)], [F(1, 2 ** 40)] * 2)
    assert res["verdict"] == "PASS", res["reason"]
    assert contains_real(res["K"][0], 1 / mp.sqrt(3)) and contains_real(res["K"][1], mp.pi)


def test_tol_control_point_width_is_tight():
    s = ts.make_kepler(1)
    st, _, n = integrate(s, ls.lohner_init(s, [F(3, 2), F(0), F(0), F(577, 1000)]), F(3), F(1, 8),
                         tol_control=True)
    assert st.t == 3
    assert max(b.width_q() for b in ls.enclosure_box(st)) < F(1, 10 ** 25)
