# -*- coding: utf-8 -*-
# 자동 생성(fakes/_gen_fakes.py) — 검사기 자기시험용 가짜 대상. 직접 고치지 말고 생성기를 고칠 것.
from _mutant import build

EXPECT = 'detect'
DESCRIPTION = '_canon 이 꼬리 0 비트를 안 지움 — 같은 값이지만 (가수, 지수) 표현이 다름'

_mutant, make_nbody_fast = build('canon_off', dyad_patches=[('    t = (m & -m).bit_length() - 1\n    if t:\n        return (m >> t, e + t)\n    return (m, e)\n', '    return (m, e)\n')])
