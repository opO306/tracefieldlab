"""계약 엄수 래퍼: fake_r3c_J_no_remainder (자동 생성; 수정하지 말 것)."""
import os, sys
from fractions import Fraction
sys.path.insert(0, '/home/zun24/r3c_check_v2_dev')
import dyad
from dyad import DI
import fake_r3c_J_no_remainder as _m
from fake_r3c_J_no_remainder import make_harmonic, make_kepler, make_nbody
import _fake_r3c_core as _core


def _dyadic_ok(q, prec):
    q = Fraction(q)
    d = q.denominator
    if d & (d - 1):
        return False
    n = abs(q.numerator)
    if n == 0:
        return True
    n >>= (n & -n).bit_length() - 1
    return n.bit_length() <= prec


def _check_state(state):
    prec = dyad.get_prec()
    entries = list(state.c) + [x for row in state.B for x in row] + [e for d in state.R for e in (d.lo_q, d.hi_q)]
    for q in entries:
        if not _dyadic_ok(q, prec):
            raise ValueError("state entry %s is not a dyadic with <= %d mantissa bits (R0 decision A: one precision per run)" % (Fraction(q), prec))


def lohner_init(system, x0):
    st = _m.lohner_init(system, x0)
    prec = dyad.get_prec()
    c = list(st.c)
    R = list(st.R)
    for i, q in enumerate(c):
        if not _dyadic_ok(q, prec):
            c2 = _core.round_dyadic(q, prec)
            R[i] = DI(Fraction(q) - c2)          # 오차를 R 로 바깥 반올림해 덮는다
            c[i] = c2
    st.c, st.R = c, R
    return st


def lohner_step(system, state, h, order):
    _check_state(state)
    return _m.lohner_step(system, state, h, order)


def enclosure_box(state):
    return _m.enclosure_box(state)
