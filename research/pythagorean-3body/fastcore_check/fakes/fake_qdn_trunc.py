# -*- coding: utf-8 -*-
# 자동 생성(fakes/_gen_fakes.py) — 검사기 자기시험용 가짜 대상. 직접 고치지 말고 생성기를 고칠 것.
from _mutant import build

EXPECT = 'detect'
DESCRIPTION = '나눗셈 하한 몫을 내림 대신 0 쪽으로 자름(GMP tdiv 이식 실수) — 분자가 음수일 때만 틀림'

_mutant, make_nbody_fast = build('qdn_trunc', dyad_patches=[('    q = (ma << k) // mb\n', '    q = (abs(ma << k) // mb) * (1 if ma >= 0 else -1)\n')])
