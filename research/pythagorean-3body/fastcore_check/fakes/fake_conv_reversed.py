# -*- coding: utf-8 -*-
# 자동 생성(fakes/_gen_fakes.py) — 검사기 자기시험용 가짜 대상. 직접 고치지 말고 생성기를 고칠 것.
from _mutant import build

EXPECT = 'detect'
DESCRIPTION = '코시 곱 합산 순서를 j=k..0 으로 뒤집음'

_mutant, make_nbody_fast = build('conv_reversed', taylor_patches=[('def _conv(a, b, k):\n', 'def _conv(a, b, k):\n    acc = a[k] * b[0]\n    for j in range(k - 1, -1, -1):\n        acc = acc + a[j] * b[k - j]\n    return acc\n\n\ndef _conv_unused(a, b, k):\n')])
