# -*- coding: utf-8 -*-
# 자동 생성(fakes/_gen_fakes.py) — 검사기 자기시험용 가짜 대상. 직접 고치지 말고 생성기를 고칠 것.
from _mutant import build

EXPECT = 'detect'
DESCRIPTION = 'w_[k] 의 나눗셈을 한 번 대신 두 번으로 쪼갬(반올림 두 번)'

_mutant, make_nbody_fast = build('pow_div_twice', taylor_patches=[('    return -(acc / (s0 * (2 * k)))\n', '    return -(acc / s0 / (2 * k))\n')])
