# -*- coding: utf-8 -*-
"""fake_dyad_ctor_no_validate.py -- check_dyad.py 자기시험용 가짜 구현.

lo>hi 입력에서 ValueError 대신 조용히 맞바꿈
(검사기 시험 전용. 제품 코드 아님.)
"""
from _fake_dyad_core import build

BUG_DESCRIPTION = 'lo>hi 입력에서 ValueError 대신 조용히 맞바꿈'
EXPECTED_CAUGHT_BY = ('R0-API',)  # 이 결함을 잡아야 하는 검사 ID(정답 구현은 빈 튜플)

_ns = build({'ctor_no_validate'})
DI = _ns.DI
set_prec = _ns.set_prec
get_prec = _ns.get_prec
