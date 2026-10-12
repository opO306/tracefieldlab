# -*- coding: utf-8 -*-
"""fake_dyad_hull_first_only.py -- check_dyad.py 자기시험용 가짜 구현.

hull의 하한을 첫 번째 구간 것만 사용
(검사기 시험 전용. 제품 코드 아님.)
"""
from _fake_dyad_core import build

BUG_DESCRIPTION = 'hull의 하한을 첫 번째 구간 것만 사용'
EXPECTED_CAUGHT_BY = ('R0-S1', 'R0-S2')  # 이 결함을 잡아야 하는 검사 ID

_ns = build({'hull_first_only'})
DI = _ns.DI
set_prec = _ns.set_prec
get_prec = _ns.get_prec
