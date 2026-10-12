import math, os, sys
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mpmath import mp
from ff_common import dyad, ts
from ff_problem import Problem
from ff_mshoot import run_ms

dyad.set_prec(128)
mp.dps = 60
R = F(1, 2 ** 40)


def contains_real(pair, x, eps=F(1, 10 ** 50)):
    q = F(mp.nstr(x, 58))
    return F(pair[0]) <= q - eps and q + eps <= F(pair[1])


def kepler():
    return Problem("kep_ms", ts.make_kepler(1), [F(3, 2), F(0), F(0), F(0)], [3], [0, 1],
                   [F(-1, 2), F(0)], F(1, 8), tol_control=True)


def test_kepler_three_segments_pass_contains_exact():
    res = run_ms(kepler(), [F(1), F(2)], [F(1 / math.sqrt(3))], F(math.pi), R, R, R)
    assert res["verdict"] == "PASS", res["reason"]
    assert res["k"] == 3 and len(res["X"]) == 1 + 4 * 2 + 1
    assert contains_real(res["K"][0], 1 / mp.sqrt(3)) and contains_real(res["K"][-1], mp.pi)


def test_single_segment_reduces_to_single_shooting():
    res = run_ms(kepler(), [], [F(1 / math.sqrt(3))], F(math.pi), R, R, R)
    assert res["verdict"] == "PASS" and res["k"] == 1 and len(res["X"]) == 2


import ff_check
from ff_prescan import choose_times


def test_retry_recovers_off_centre_and_radius_bounded():
    res = run_ms(kepler(), [F(1), F(2)], [F(1 / math.sqrt(3)) + F(1, 10 ** 11)], F(math.pi), R, R, R)
    assert res["attempts"][0]["verdict"] == "REFUSED"
    assert res["verdict"] == "PASS", res["attempts"]
    assert max(F(v) for v in res["r"]) <= 64 * R
    assert contains_real(res["K"][0], 1 / mp.sqrt(3))


def test_independent_check_on_ms_record():
    res = run_ms(kepler(), [F(1), F(2)], [F(1 / math.sqrt(3))], F(math.pi), R, R, R)
    chk = ff_check.check(res)
    assert chk["verdict"] == "PASS" and chk["unique"]


def test_choose_times_increasing_dyadic():
    ts_ = choose_times([1, 0.8, 0.8], 0.0009114239445, 0.3019805958, 0.91431242, 8)
    assert len(ts_) == 7 and all(a < b for a, b in zip(ts_, ts_[1:])) and 0 < ts_[0] and ts_[-1] < 0.9143
    assert all((t * 2 ** 16).is_integer() for t in ts_)


def test_retry_shrinks_radius_when_box_too_large():
    # 중심은 정확하지만 r이 너무 크면 K가 X를 못 넘는다 → 재시도는 r을 줄여야 한다(설계서 11절 새 해석)
    res = run_ms(kepler(), [], [F(1 / math.sqrt(3))], F(math.pi), F(1, 2 ** 8), F(1, 2 ** 8), F(1, 2 ** 8), retries=3)
    assert res["attempts"][0]["verdict"] == "REFUSED", res["attempts"]
    assert res["verdict"] == "PASS", res["attempts"]
    assert max(F(v) for v in res["r"]) < F(1, 2 ** 8)
