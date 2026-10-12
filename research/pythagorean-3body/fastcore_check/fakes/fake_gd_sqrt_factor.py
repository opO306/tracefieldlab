# -*- coding: utf-8 -*-
# 자동 생성(fakes/_gen_fakes.py) — 검사기 자기시험용 가짜 대상. 직접 고치지 말고 생성기를 고칠 것.
from _mutant import build

EXPECT = 'detect'
DESCRIPTION = '_GD.sqrt() 의 미분 ds/(2 sqrt s) 에서 2를 빠뜨림'

_mutant, make_nbody_fast = build('gd_sqrt_factor', lohner_patches=[('        r2 = r * 2\n', '        r2 = r\n')])
