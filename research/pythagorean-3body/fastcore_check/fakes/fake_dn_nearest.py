# -*- coding: utf-8 -*-
# 자동 생성(fakes/_gen_fakes.py) — 검사기 자기시험용 가짜 대상. 직접 고치지 말고 생성기를 고칠 것.
from _mutant import build

EXPECT = 'detect'
DESCRIPTION = '_dn 이 내림 대신 가장 가까운 값으로 반올림(반올림 동률 포함)'

_mutant, make_nbody_fast = build('dn_nearest', dyad_patches=[('        m >>= s\n        e += s\n', '        m = (m + (1 << (s - 1))) >> s\n        e += s\n')])
