# -*- coding: utf-8 -*-
# 자동 생성(fakes/_gen_fakes.py) — 검사기 자기시험용 가짜 대상. 직접 고치지 말고 생성기를 고칠 것.
from _mutant import build

EXPECT = 'detect'
DESCRIPTION = '_dn(하한 내림)이 음수에서 0 쪽으로 자름(floor 대신 truncate) — GMP tdiv 를 쓴 이식 실수'

_mutant, make_nbody_fast = build('dn_trunc', dyad_patches=[('        m >>= s\n        e += s\n', '        m = (m >> s) if m >= 0 else -((-m) >> s)\n        e += s\n')])
