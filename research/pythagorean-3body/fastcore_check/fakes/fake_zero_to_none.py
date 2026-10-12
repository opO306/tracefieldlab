# -*- coding: utf-8 -*-
# 자동 생성(fakes/_gen_fakes.py) — 검사기 자기시험용 가짜 대상. 직접 고치지 말고 생성기를 고칠 것.
from _mutant import build

EXPECT = 'detect'
DESCRIPTION = '결과 기울기의 정확한 0 인 DI(0) 을 None 으로 바꿈'

from lohner_step import _GD as _G


def _sparsify(res):
    out = []
    for comp in res:
        row = []
        for e in comp:
            if type(e) is _G:
                e = _G(e.v, [None if (g is not None and g._l == (0, 0) and g._h == (0, 0)) else g
                             for g in e.g])
            row.append(e)
        out.append(row)
    return out


_mutant, make_nbody_fast = build('zero_to_none', post=_sparsify)
