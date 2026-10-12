# -*- coding: utf-8 -*-
# 자동 생성(fakes/_gen_fakes.py) — 검사기 자기시험용 가짜 대상. 직접 고치지 말고 생성기를 고칠 것.
from _mutant import build

EXPECT = 'detect'
DESCRIPTION = '곱셈 결합 순서 변경 (s*w)*c -> s*(w*c)'

_mutant, make_nbody_fast = build('assoc_order', taylor_patches=[('    acc = (s[1] * w[k - 1]) * (1 + 2 * k)\n', '    acc = s[1] * (w[k - 1] * (1 + 2 * k))\n')])
