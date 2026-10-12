# -*- coding: utf-8 -*-
"""fake_dyad_correct.py -- check_dyad.py 자기시험용 정답 구현.

정답 구현: 모든 끝점을 정확한 Fraction으로 계산 후 prec 비트로 바깥쪽 반올림 (결함 없음)
(검사기 시험 전용. 제품 코드 아님.)
"""
from _fake_dyad_core import build

BUG_DESCRIPTION = '정답 구현: 모든 끝점을 정확한 Fraction으로 계산 후 prec 비트로 바깥쪽 반올림 (결함 없음)'
EXPECTED_CAUGHT_BY = ()  # 이 결함을 잡아야 하는 검사 ID(정답 구현은 빈 튜플)

_ns = build(set())
DI = _ns.DI
set_prec = _ns.set_prec
get_prec = _ns.get_prec
