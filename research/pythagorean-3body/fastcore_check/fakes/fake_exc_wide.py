# -*- coding: utf-8 -*-
# 자동 생성(fakes/_gen_fakes.py) — 검사기 자기시험용 가짜 대상. 직접 고치지 말고 생성기를 고칠 것.
from _mutant import build

EXPECT = 'detect'
DESCRIPTION = '특이점에서 예외 대신 거대한 구간을 돌려줌'

from dyad import DI as _DI


def _wide(system, K):
    big = _DI._new((-1, 8192), (1, 8192))
    return [[big] * (K + 1) for _ in range(system.dim)]


_mutant, make_nbody_fast = build('exc_wide', on_singular=_wide)
