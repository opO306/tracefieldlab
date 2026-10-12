# -*- coding: utf-8 -*-
"""fake_dyad_carry_bug.py -- check_dyad.py 자기시험용 가짜 구현.

올림 시 자릿수가 올라가는 경우(가수 전부 1)에 올림수를 잃고 내림값 반환
(검사기 시험 전용. 제품 코드 아님.)
"""
from _fake_dyad_core import build

BUG_DESCRIPTION = '올림 시 자릿수가 올라가는 경우(가수 전부 1)에 올림수를 잃고 내림값 반환'
EXPECTED_CAUGHT_BY = ('R0-S1', 'R0-S2')  # 이 결함을 잡아야 하는 검사 ID

_ns = build({'carry_bug'})
DI = _ns.DI
set_prec = _ns.set_prec
get_prec = _ns.get_prec
