# -*- coding: utf-8 -*-
# 자동 생성(fakes/_gen_fakes.py) — 검사기 자기시험용 가짜 대상. 직접 고치지 말고 생성기를 고칠 것.
from _mutant import build

EXPECT = 'detect'
DESCRIPTION = '특이점에서 _Singular 대신 ZeroDivisionError'

def _zde(system, K):
    raise ZeroDivisionError('singular')


_mutant, make_nbody_fast = build('exc_type', on_singular=_zde)
