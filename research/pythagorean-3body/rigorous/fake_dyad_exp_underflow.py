# -*- coding: utf-8 -*-
"""fake_dyad_exp_underflow.py -- check_dyad.py 자기시험용 가짜 구현.

2^-1000 미만 양수 상한을 0으로 플러시(부동소수점식 언더플로)
(검사기 시험 전용. 제품 코드 아님.)
"""
from _fake_dyad_core import build

BUG_DESCRIPTION = '2^-1000 미만 양수 상한을 0으로 플러시(부동소수점식 언더플로)'
EXPECTED_CAUGHT_BY = ('R0-S1', 'R0-S2')  # 이 결함을 잡아야 하는 검사 ID

_ns = build({'exp_underflow'})
DI = _ns.DI
set_prec = _ns.set_prec
get_prec = _ns.get_prec
