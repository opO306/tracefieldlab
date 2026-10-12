# -*- coding: utf-8 -*-
"""fake_dyad_contains_open.py -- check_dyad.py 자기시험용 가짜 구현.

contains가 상한 끝점을 포함하지 않음(반열린 구간)
(검사기 시험 전용. 제품 코드 아님.)
"""
from _fake_dyad_core import build

BUG_DESCRIPTION = 'contains가 상한 끝점을 포함하지 않음(반열린 구간)'
EXPECTED_CAUGHT_BY = ('R0-API',)  # 이 결함을 잡아야 하는 검사 ID(정답 구현은 빈 튜플)

_ns = build({'contains_open'})
DI = _ns.DI
set_prec = _ns.set_prec
get_prec = _ns.get_prec
