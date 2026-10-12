# -*- coding: utf-8 -*-
# 자동 생성(fakes/_gen_fakes.py) — 검사기 자기시험용 가짜 대상. 직접 고치지 말고 생성기를 고칠 것.
from _mutant import build

EXPECT = 'detect'
DESCRIPTION = '내부 계산을 요구 정밀도보다 1비트 높게'

_mutant, make_nbody_fast = build('prec_plus_one', prec_delta=1)
