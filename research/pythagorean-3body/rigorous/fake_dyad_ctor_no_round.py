# -*- coding: utf-8 -*-
"""fake_dyad_ctor_no_round.py -- check_dyad.py 자기시험용 가짜 구현.

생성자에서 반올림 누락(끝점 그대로 저장)
(검사기 시험 전용. 제품 코드 아님.)
"""
from _fake_dyad_core import build

BUG_DESCRIPTION = '생성자에서 반올림 누락(끝점 그대로 저장)'
EXPECTED_CAUGHT_BY = ('R0-S4', 'R0-S1', 'R0-S2')  # 이 결함을 잡아야 하는 검사 ID(정답 구현은 빈 튜플)

_ns = build({'ctor_no_round'})
DI = _ns.DI
set_prec = _ns.set_prec
get_prec = _ns.get_prec
