# -*- coding: utf-8 -*-
# 자동 생성(fakes/_gen_fakes.py) — 검사기 자기시험용 가짜 대상. 직접 고치지 말고 생성기를 고칠 것.
from _mutant import build

EXPECT = 'detect'
DESCRIPTION = '_up(상한 올림)이 올림 대신 내림'

_mutant, make_nbody_fast = build('up_floor', dyad_patches=[('        m = -((-m) >> s)\n', '        m = m >> s\n')])
