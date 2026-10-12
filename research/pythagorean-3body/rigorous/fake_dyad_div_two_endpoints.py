# -*- coding: utf-8 -*-
"""fake_dyad_div_two_endpoints.py -- check_dyad.py 자기시험용 가짜 구현.

나눗셈에서 끝점 몫 4개 대신 두 개만 사용
(검사기 시험 전용. 제품 코드 아님.)
"""
from _fake_dyad_core import build

BUG_DESCRIPTION = '나눗셈에서 끝점 몫 4개 대신 두 개만 사용'
EXPECTED_CAUGHT_BY = ('R0-S1', 'R0-S2')  # 이 결함을 잡아야 하는 검사 ID(정답 구현은 빈 튜플)

_ns = build({'div_two_endpoints'})
DI = _ns.DI
set_prec = _ns.set_prec
get_prec = _ns.get_prec
