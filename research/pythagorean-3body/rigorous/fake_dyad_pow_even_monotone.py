# -*- coding: utf-8 -*-
"""fake_dyad_pow_even_monotone.py -- check_dyad.py 자기시험용 가짜 구현.

짝수 거듭제곱을 단조로 취급: 0을 걸친 구간의 하한이 0이 아님(불포함)
(검사기 시험 전용. 제품 코드 아님.)
"""
from _fake_dyad_core import build

BUG_DESCRIPTION = '짝수 거듭제곱을 단조로 취급: 0을 걸친 구간의 하한이 0이 아님(불포함)'
EXPECTED_CAUGHT_BY = ('R0-S1', 'R0-S2')  # 이 결함을 잡아야 하는 검사 ID(정답 구현은 빈 튜플)

_ns = build({'pow_even_monotone'})
DI = _ns.DI
set_prec = _ns.set_prec
get_prec = _ns.get_prec
