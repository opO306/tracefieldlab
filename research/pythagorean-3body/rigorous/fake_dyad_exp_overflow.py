# -*- coding: utf-8 -*-
"""fake_dyad_exp_overflow.py -- check_dyad.py 자기시험용 가짜 구현.

2^1000 초과 상한을 2^1000으로 포화(부동소수점식 오버플로)
(검사기 시험 전용. 제품 코드 아님.)
"""
from _fake_dyad_core import build

BUG_DESCRIPTION = '2^1000 초과 상한을 2^1000으로 포화(부동소수점식 오버플로)'
EXPECTED_CAUGHT_BY = ('R0-S1', 'R0-S2')  # 이 결함을 잡아야 하는 검사 ID(정답 구현은 빈 튜플)

_ns = build({'exp_overflow'})
DI = _ns.DI
set_prec = _ns.set_prec
get_prec = _ns.get_prec
