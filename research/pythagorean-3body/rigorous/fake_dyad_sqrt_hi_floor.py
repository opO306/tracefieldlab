# -*- coding: utf-8 -*-
"""fake_dyad_sqrt_hi_floor.py -- check_dyad.py 자기시험용 가짜 구현.

sqrt 상한을 올림 대신 내림
(검사기 시험 전용. 제품 코드 아님.)
"""
from _fake_dyad_core import build

BUG_DESCRIPTION = 'sqrt 상한을 올림 대신 내림'
EXPECTED_CAUGHT_BY = ('R0-S3', 'R0-S1', 'R0-S2')  # 이 결함을 잡아야 하는 검사 ID(정답 구현은 빈 튜플)

_ns = build({'sqrt_hi_floor'})
DI = _ns.DI
set_prec = _ns.set_prec
get_prec = _ns.get_prec
