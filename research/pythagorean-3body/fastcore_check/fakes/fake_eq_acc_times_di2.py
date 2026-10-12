# -*- coding: utf-8 -*-
# 자동 생성(fakes/_gen_fakes.py) — 검사기 자기시험용 가짜 대상. 직접 고치지 말고 생성기를 고칠 것.
from _mutant import build

EXPECT = 'equivalent'
DESCRIPTION = '[동치] 2배를 정수 2 대신 DI(2) 로'

_mutant, make_nbody_fast = build('eq_acc_times_di2', taylor_patches=[('        acc = acc * 2', '        acc = acc * DI(2)')])
