# -*- coding: utf-8 -*-
# 자동 생성(fakes/_gen_fakes.py) — 검사기 자기시험용 가짜 대상. 직접 고치지 말고 생성기를 고칠 것.
from _mutant import build

EXPECT = 'equivalent'
DESCRIPTION = '[동치] 2배를 acc*2 대신 acc+acc 로 (2의 거듭제곱 곱은 반올림과 교환되므로 비트 동일)'

_mutant, make_nbody_fast = build('eq_acc_plus_acc', taylor_patches=[('        acc = acc * 2', '        acc = acc + acc')])
