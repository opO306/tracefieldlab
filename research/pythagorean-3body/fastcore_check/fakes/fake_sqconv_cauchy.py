# -*- coding: utf-8 -*-
# 자동 생성(fakes/_gen_fakes.py) — 검사기 자기시험용 가짜 대상. 직접 고치지 말고 생성기를 고칠 것.
from _mutant import build

EXPECT = 'detect'
DESCRIPTION = '제곱 코시 곱에서 대칭 쌍을 한 번 계산해 2배 하지 않고 두 번 곱해 모두 더함'

_mutant, make_nbody_fast = build('sqconv_cauchy', taylor_patches=[('def _sq_conv(a, k):\n', 'def _sq_conv(a, k):\n    return _conv(a, a, k)\n\n\ndef _sq_conv_unused(a, k):\n')])
