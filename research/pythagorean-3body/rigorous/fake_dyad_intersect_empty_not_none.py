# -*- coding: utf-8 -*-
"""fake_dyad_intersect_empty_not_none.py -- check_dyad.py 자기시험용 가짜 구현.

교집합이 비어도 None 대신 뒤집힌 구간을 반환
(검사기 시험 전용. 제품 코드 아님.)
"""
from _fake_dyad_core import build

BUG_DESCRIPTION = '교집합이 비어도 None 대신 뒤집힌 구간을 반환'
EXPECTED_CAUGHT_BY = ('R0-API',)  # 이 결함을 잡아야 하는 검사 ID(정답 구현은 빈 튜플)

_ns = build({'intersect_empty_not_none'})
DI = _ns.DI
set_prec = _ns.set_prec
get_prec = _ns.get_prec
