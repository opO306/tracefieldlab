# -*- coding: utf-8 -*-
# 자동 생성(fakes/_gen_fakes.py) — 검사기 자기시험용 가짜 대상. 직접 고치지 말고 생성기를 고칠 것.
from _mutant import build

EXPECT = 'detect'
DESCRIPTION = '_GD.sqr() 의 미분에서 2배를 빠뜨림'

_mutant, make_nbody_fast = build('gd_sqr_factor', lohner_patches=[('        return _GD(av.sqr(), [None if a is None else (av * a) * 2 for a in self.g])\n', '        return _GD(av.sqr(), [None if a is None else (av * a) for a in self.g])\n')])
