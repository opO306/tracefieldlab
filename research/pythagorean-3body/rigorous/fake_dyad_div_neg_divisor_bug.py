# -*- coding: utf-8 -*-
"""fake_dyad_div_neg_divisor_bug.py -- check_dyad.py 자기시험용 가짜 구현.

음수 제수 구간에서 부호를 빠뜨림
(검사기 시험 전용. 제품 코드 아님.)
"""
from _fake_dyad_core import build

BUG_DESCRIPTION = '음수 제수 구간에서 부호를 빠뜨림'
EXPECTED_CAUGHT_BY = ('R0-S1', 'R0-S2')  # 이 결함을 잡아야 하는 검사 ID

_ns = build({'div_neg_divisor_bug'})
DI = _ns.DI
set_prec = _ns.set_prec
get_prec = _ns.get_prec
