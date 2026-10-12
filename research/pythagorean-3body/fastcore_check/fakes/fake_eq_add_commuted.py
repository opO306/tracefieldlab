# -*- coding: utf-8 -*-
# 자동 생성(fakes/_gen_fakes.py) — 검사기 자기시험용 가짜 대상. 직접 고치지 말고 생성기를 고칠 것.
from _mutant import build

EXPECT = 'equivalent'
DESCRIPTION = '[동치] 덧셈의 피연산자 순서 교환(구간 덧셈은 교환법칙이 비트 단위로 성립)'

_mutant, make_nbody_fast = build('eq_add_commuted', taylor_patches=[('                acc[2 * a] = acc[2 * a] + M[b] * ux\n', '                acc[2 * a] = M[b] * ux + acc[2 * a]\n')])
