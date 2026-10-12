# -*- coding: utf-8 -*-
# 자동 생성(fakes/_gen_fakes.py) — 검사기 자기시험용 가짜 대상. 직접 고치지 말고 생성기를 고칠 것.
from _mutant import build

EXPECT = 'detect'
DESCRIPTION = '결과를 리스트가 아닌 튜플로 돌려줌(모양 차이만; 값은 같음)'

def _tup(res):
    return tuple(tuple(c) for c in res)


_mutant, make_nbody_fast = build('container_tuple', post=_tup)
