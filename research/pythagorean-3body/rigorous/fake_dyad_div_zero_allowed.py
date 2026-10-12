# -*- coding: utf-8 -*-
"""fake_dyad_div_zero_allowed.py -- check_dyad.py 자기시험용 가짜 구현.

0을 걸친 제수 구간에서 ZeroDivisionError 대신 값을 돌려줌
(검사기 시험 전용. 제품 코드 아님.)
"""
from _fake_dyad_core import build

BUG_DESCRIPTION = '0을 걸친 제수 구간에서 ZeroDivisionError 대신 값을 돌려줌'
EXPECTED_CAUGHT_BY = ('R0-API',)  # 이 결함을 잡아야 하는 검사 ID(정답 구현은 빈 튜플)

_ns = build({'div_zero_allowed'})
DI = _ns.DI
set_prec = _ns.set_prec
get_prec = _ns.get_prec
