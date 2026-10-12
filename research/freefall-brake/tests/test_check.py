import copy, math, os, sys
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ff_common import dyad, ts
from ff_problem import Problem, run_problem
import ff_check

dyad.set_prec(128)


def _rec():
    p = Problem("C1a_harmonic", ts.make_harmonic(), [F(1), F(0)], [], [1], [F(0)], F(1, 8))
    return run_problem(p, [F(math.pi)], [F(1, 2 ** 40)])


def test_independent_check_agrees():
    rec = _rec()
    assert rec["verdict"] == "PASS"
    assert ff_check.check(rec)["verdict"] == "PASS"


def test_tampered_record_fails():
    rec = copy.deepcopy(_rec())
    lo, hi = (F(x) for x in rec["gm"][0])
    rec["gm"][0] = [str(lo + 1), str(hi + 1)]     # g(m) shifted by 1
    assert ff_check.check(rec)["verdict"] == "FAIL"


def test_widened_box_fails():
    rec = copy.deepcopy(_rec())
    rec["X"] = [[str(F(lo) - F(1, 1000)), str(F(hi) + F(1, 1000))] for lo, hi in rec["X"]]
    assert ff_check.check(rec)["verdict"] == "FAIL"


def test_uniqueness_norm_reported():
    r = ff_check.check(_rec())
    assert r["unique"] is True and F(r["norm_I_minus_CDg"]) < 1
