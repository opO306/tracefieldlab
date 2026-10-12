# -*- coding: utf-8 -*-
# 자동 생성(fakes/_gen_fakes.py) — 검사기 자기시험용 가짜 대상. 직접 고치지 말고 생성기를 고칠 것.
from _mutant import build

EXPECT = 'detect'
DESCRIPTION = 'sqr() 를 x*x 로 대체 — 0 을 가로지르는 구간에서 하한이 음수로 나옴'

_mutant, make_nbody_fast = build('sqr_as_mul', dyad_patches=[('        return self ** 2\n', '        return self * self\n')])
