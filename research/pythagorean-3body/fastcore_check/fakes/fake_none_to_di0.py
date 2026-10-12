# -*- coding: utf-8 -*-
# 자동 생성(fakes/_gen_fakes.py) — 검사기 자기시험용 가짜 대상. 직접 고치지 말고 생성기를 고칠 것.
from _mutant import build

EXPECT = 'detect'
DESCRIPTION = '_GD 입력의 None(정확히 0)을 DI(0) 으로 바꿔 계산 — 값은 같지만 None 위치가 달라짐'

from dyad import DI as _DI
from lohner_step import _GD as _G


def _dense(X):
    return [(_G(x.v, [_DI(0) if g is None else g for g in x.g]) if type(x) is _G else x)
            for x in X]


_mutant, make_nbody_fast = build('none_to_di0', pre=_dense)
