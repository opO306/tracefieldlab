# -*- coding: utf-8 -*-
# 자동 생성(fakes/_gen_fakes.py) — 검사기 자기시험용 가짜 대상. 직접 고치지 말고 생성기를 고칠 것.
from _mutant import build

EXPECT = 'detect'
DESCRIPTION = '_GD 나눗셈 미분 (da - q*db)/b 를 da/b - (q*db)/b 로 쪼갬'

_mutant, make_nbody_fast = build('gd_div_split', lohner_patches=[('                    g.append((a - q * b) / bv)\n', '                    g.append(a / bv - (q * b) / bv)\n')])
