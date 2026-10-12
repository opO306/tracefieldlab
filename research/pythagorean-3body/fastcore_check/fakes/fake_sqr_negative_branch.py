# -*- coding: utf-8 -*-
# 자동 생성(fakes/_gen_fakes.py) — 검사기 자기시험용 가짜 대상. 직접 고치지 말고 생성기를 고칠 것.
from _mutant import build

EXPECT = 'detect'
DESCRIPTION = 'x**2 에서 구간이 전부 음수일 때 끝점을 맞바꿈(하한이 상한보다 커짐)'

_mutant, make_nbody_fast = build('sqr_negative_branch', dyad_patches=[('        return DI._new(_dn(*_powx(hi, n), p), _up(*_powx(lo, n), p))\n', '        return DI._new(_dn(*_powx(lo, n), p), _up(*_powx(hi, n), p))\n')])
