# -*- coding: utf-8 -*-
# 자동 생성(fakes/_gen_fakes.py) — 검사기 자기시험용 가짜 대상. 직접 고치지 말고 생성기를 고칠 것.
from _mutant import build

EXPECT = 'detect'
DESCRIPTION = '_sq_conv 의 가운데 항을 sqr 대신 곱셈으로'

_mutant, make_nbody_fast = build('sqconv_mid_mul', taylor_patches=[('        m = a[k // 2].sqr()\n', '        m = a[k // 2] * a[k // 2]\n')])
