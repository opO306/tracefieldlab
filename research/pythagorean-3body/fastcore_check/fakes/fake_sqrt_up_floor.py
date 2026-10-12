# -*- coding: utf-8 -*-
# 자동 생성(fakes/_gen_fakes.py) — 검사기 자기시험용 가짜 대상. 직접 고치지 말고 생성기를 고칠 것.
from _mutant import build

EXPECT = 'detect'
DESCRIPTION = 'sqrt 상한이 완전제곱이 아닐 때 +1 을 빼먹음(올림 누락)'

_mutant, make_nbody_fast = build('sqrt_up_floor', dyad_patches=[('    if not exact:\n        s += 1\n', '    if False:\n        s += 1\n')])
