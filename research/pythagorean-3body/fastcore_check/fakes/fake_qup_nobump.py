# -*- coding: utf-8 -*-
# 자동 생성(fakes/_gen_fakes.py) — 검사기 자기시험용 가짜 대상. 직접 고치지 말고 생성기를 고칠 것.
from _mutant import build

EXPECT = 'detect'
DESCRIPTION = '나눗셈 상한 몫에서 나머지가 있어도 +1 하지 않음'

_mutant, make_nbody_fast = build('qup_nobump', dyad_patches=[('    if r:\n        q += 1\n    return _up(q, ea - eb - k, p)\n', '    return _up(q, ea - eb - k, p)\n')])
