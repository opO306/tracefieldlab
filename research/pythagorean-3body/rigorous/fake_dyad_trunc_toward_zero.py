# -*- coding: utf-8 -*-
"""fake_dyad_trunc_toward_zero.py -- check_dyad.py 자기시험용 가짜 구현.

바깥쪽 대신 0 쪽으로 자름(음수 hi/양수 lo는 맞고 반대는 틀림)
(검사기 시험 전용. 제품 코드 아님.)
"""
from _fake_dyad_core import build

BUG_DESCRIPTION = '바깥쪽 대신 0 쪽으로 자름(음수 hi/양수 lo는 맞고 반대는 틀림)'
EXPECTED_CAUGHT_BY = ('R0-S1', 'R0-S2')  # 이 결함을 잡아야 하는 검사 ID(정답 구현은 빈 튜플)

_ns = build({'trunc_toward_zero'})
DI = _ns.DI
set_prec = _ns.set_prec
get_prec = _ns.get_prec
