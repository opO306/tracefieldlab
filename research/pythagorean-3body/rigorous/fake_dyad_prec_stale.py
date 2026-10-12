# -*- coding: utf-8 -*-
"""fake_dyad_prec_stale.py -- check_dyad.py 자기시험용 가짜 구현.

연산이 현재 정밀도 대신 생성 당시 정밀도로 반올림(정밀도 변경 후 가수 비트 위반)
(검사기 시험 전용. 제품 코드 아님.)
"""
from _fake_dyad_core import build

BUG_DESCRIPTION = '연산이 현재 정밀도 대신 생성 당시 정밀도로 반올림(정밀도 변경 후 가수 비트 위반)'
EXPECTED_CAUGHT_BY = ('R0-S4',)  # 이 결함을 잡아야 하는 검사 ID

_ns = build({'prec_stale'})
DI = _ns.DI
set_prec = _ns.set_prec
get_prec = _ns.get_prec
