# -*- coding: utf-8 -*-
"""fake_dyad_add_no_round.py -- check_dyad.py 자기시험용 가짜 구현.

덧셈 결과를 반올림하지 않음
(검사기 시험 전용. 제품 코드 아님.)
"""
from _fake_dyad_core import build

BUG_DESCRIPTION = '덧셈 결과를 반올림하지 않음'
EXPECTED_CAUGHT_BY = ('R0-S4',)  # 이 결함을 잡아야 하는 검사 ID(정답 구현은 빈 튜플)

_ns = build({'add_no_round'})
DI = _ns.DI
set_prec = _ns.set_prec
get_prec = _ns.get_prec
