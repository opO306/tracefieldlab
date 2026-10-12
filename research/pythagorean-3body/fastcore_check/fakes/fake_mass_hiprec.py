# -*- coding: utf-8 -*-
# 자동 생성(fakes/_gen_fakes.py) — 검사기 자기시험용 가짜 대상. 직접 고치지 말고 생성기를 고칠 것.
from _mutant import build

EXPECT = 'detect'
DESCRIPTION = '질량 상수 DI(m) 를 호출 정밀도보다 64비트 높게 만듦'

_mutant, make_nbody_fast = build('mass_hiprec', taylor_patches=[('        M = [DI(m) for m in masses]', '        _q = dyad.get_prec(); dyad.set_prec(_q + 64); M = [DI(m) for m in masses]; dyad.set_prec(_q)')])
