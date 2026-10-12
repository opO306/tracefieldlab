# -*- coding: utf-8 -*-
# 자동 생성(fakes/_gen_fakes.py) — 검사기 자기시험용 가짜 대상. 직접 고치지 말고 생성기를 고칠 것.
from _mutant import build

EXPECT = 'detect'
DESCRIPTION = '정밀도가 128비트를 넘으면 128로 잘라서 계산(고정 폭 한계 이식 실수)'

_mutant, make_nbody_fast = build('prec_clamp128', prec_fn=lambda p: min(p, 128))
