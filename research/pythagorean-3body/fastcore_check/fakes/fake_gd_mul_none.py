# -*- coding: utf-8 -*-
# 자동 생성(fakes/_gen_fakes.py) — 검사기 자기시험용 가짜 대상. 직접 고치지 말고 생성기를 고칠 것.
from _mutant import build

EXPECT = 'detect'
DESCRIPTION = '_GD 곱셈에서 한쪽 기울기가 None 일 때 다른 쪽 항을 버림'

_mutant, make_nbody_fast = build('gd_mul_none', lohner_patches=[('                    g.append(None if b is None else av * b)\n', '                    g.append(None)\n')])
