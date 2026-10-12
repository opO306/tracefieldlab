# -*- coding: utf-8 -*-
"""fake_dyad_sqrt_neg_clamp.py -- check_dyad.py 자기시험용 가짜 구현.

음수 하한 sqrt에서 ValueError 대신 하한을 0으로 잘라 계산
(검사기 시험 전용. 제품 코드 아님.)
"""
from _fake_dyad_core import build

BUG_DESCRIPTION = '음수 하한 sqrt에서 ValueError 대신 하한을 0으로 잘라 계산'
EXPECTED_CAUGHT_BY = ('R0-API',)  # 이 결함을 잡아야 하는 검사 ID(정답 구현은 빈 튜플)

_ns = build({'sqrt_neg_clamp'})
DI = _ns.DI
set_prec = _ns.set_prec
get_prec = _ns.get_prec
