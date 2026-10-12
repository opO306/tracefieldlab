# -*- coding: utf-8 -*-
# 자동 생성(fakes/_gen_fakes.py) — 검사기 자기시험용 가짜 대상. 직접 고치지 말고 생성기를 고칠 것.
from _mutant import build

EXPECT = 'detect'
DESCRIPTION = '_GD 의 (상수 / _GD) 미분에서 음수 부호를 빠뜨림'

_mutant, make_nbody_fast = build('gd_rtruediv_sign', lohner_patches=[('        return _GD(q, [None if a is None else -(qu * a) for a in self.g])\n', '        return _GD(q, [None if a is None else (qu * a) for a in self.g])\n')])
