# -*- coding: utf-8 -*-
# 자동 생성(fakes/_gen_fakes.py) — 검사기 자기시험용 가짜 대상. 직접 고치지 말고 생성기를 고칠 것.
from _mutant import build

EXPECT = 'detect'
DESCRIPTION = '_GD 의 (상수 - _GD) 에서 기울기 부호를 뒤집지 않음'

_mutant, make_nbody_fast = build('gd_rsub_sign', lohner_patches=[('        return _GD(o - self.v, [None if a is None else -a for a in self.g])\n', '        return _GD(o - self.v, [None if a is None else a for a in self.g])\n')])
