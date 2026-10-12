# -*- coding: utf-8 -*-
# 자동 생성(fakes/_gen_fakes.py) — 검사기 자기시험용 가짜 대상. 직접 고치지 말고 생성기를 고칠 것.
from _mutant import build

EXPECT = 'detect'
DESCRIPTION = '계산 뒤에 전역 정밀도를 바꿔 놓음'

import dyad as _dy


def _bump(res):
    _dy.set_prec(_dy.get_prec() + 1)
    return res


_mutant, make_nbody_fast = build('prec_side_effect', post=_bump)
