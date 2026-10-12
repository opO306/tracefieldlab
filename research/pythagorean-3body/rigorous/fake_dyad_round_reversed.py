# -*- coding: utf-8 -*-
"""fake_dyad_round_reversed.py -- check_dyad.py 자기시험용 가짜 구현.

반올림 방향 반대: lo를 올림, hi를 내림
(검사기 시험 전용. 제품 코드 아님.)
"""
from _fake_dyad_core import build

BUG_DESCRIPTION = '반올림 방향 반대: lo를 올림, hi를 내림'
EXPECTED_CAUGHT_BY = ('R0-S1', 'R0-S2')  # 이 결함을 잡아야 하는 검사 ID(정답 구현은 빈 튜플)

_ns = build({'round_reversed'})
DI = _ns.DI
set_prec = _ns.set_prec
get_prec = _ns.get_prec
