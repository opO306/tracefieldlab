# -*- coding: utf-8 -*-
# 자동 생성(fakes/_gen_fakes.py) — 검사기 자기시험용 가짜 대상. 직접 고치지 말고 생성기를 고칠 것.
from _mutant import build

EXPECT = 'detect'
DESCRIPTION = '나눗셈 하한 몫을 내림 대신 올림으로 계산(반올림 방향 효과가 드문 변형)'

_mutant, make_nbody_fast = build('qdn_ceil', dyad_patches=[('    q = (ma << k) // mb\n', '    q = -((-(ma << k)) // mb)\n')])
