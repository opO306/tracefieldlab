# -*- coding: utf-8 -*-
# 자동 생성(fakes/_gen_fakes.py) — 검사기 자기시험용 가짜 대상. 직접 고치지 말고 생성기를 고칠 것.
from _mutant import build

EXPECT = 'detect'
DESCRIPTION = '결과 끝점의 지수를 부호 있는 32비트로 감아 돌려줌(어댑터가 지수를 int32 로 저장하는 실수) — 지수가 ±2^31 을 넘는 입력에서만 틀림'

import dyad as _dy
from lohner_step import _GD as _G


def _w(e):
    return ((e + 0x80000000) % 0x100000000) - 0x80000000


def _ep(t):
    return (t[0], _w(t[1])) if t[0] else t


def _wrap(x):
    if x is None:
        return None
    if type(x) is list:
        return [_wrap(i) for i in x]
    if type(x) is _dy.DI:
        return _dy.DI._new(_ep(x._l), _ep(x._h))
    return _G(_wrap(x.v), [_wrap(g) for g in x.g])


_mutant, make_nbody_fast = build('exp_wrap32', post=_wrap)
