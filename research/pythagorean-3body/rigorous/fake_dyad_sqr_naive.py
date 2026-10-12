# -*- coding: utf-8 -*-
"""fake_dyad_sqr_naive.py -- check_dyad.py 자기시험용 가짜 구현.

sqr를 x*x로 구현: 0을 걸치면 하한이 음수(계약: 하한 0)
(검사기 시험 전용. 제품 코드 아님.)
"""
from _fake_dyad_core import build

BUG_DESCRIPTION = 'sqr를 x*x로 구현: 0을 걸치면 하한이 음수(계약: 하한 0)'
EXPECTED_CAUGHT_BY = ('R0-API',)  # 이 결함을 잡아야 하는 검사 ID(정답 구현은 빈 튜플)

_ns = build({'sqr_naive'})
DI = _ns.DI
set_prec = _ns.set_prec
get_prec = _ns.get_prec
