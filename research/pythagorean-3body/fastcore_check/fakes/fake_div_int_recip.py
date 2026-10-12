# -*- coding: utf-8 -*-
# 자동 생성(fakes/_gen_fakes.py) — 검사기 자기시험용 가짜 대상. 직접 고치지 말고 생성기를 고칠 것.
from _mutant import build

EXPECT = 'detect'
DESCRIPTION = 'a/(k+1) 를 a * DI(1/(k+1)) 로(역수를 먼저 반올림)'

_mutant, make_nbody_fast = build('div_int_recip', taylor_patches=[('    return a / n\n', '    return a * DI(Fraction(1, n))\n')])
