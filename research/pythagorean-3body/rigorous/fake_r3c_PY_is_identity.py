"""R3c 검사기 자기시험 가짜 구현: PY_is_identity - 결함 아님(유효한 변형): 중복인 나머지 미분 몫에서 PY 대신 단위행렬 사용.  (_fake_r3c_core.py 위에서 fault 하나만 고정)"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _fake_r3c_core as core
from _fake_r3c_core import make_harmonic, make_kepler, make_nbody

FAULT = 'PY_is_identity'


def lohner_init(system, x0):
    return core.lohner_init(system, x0)


def lohner_step(system, state, h, order):
    return core.lohner_step(system, state, h, order, fault=FAULT)


def enclosure_box(state):
    return core.enclosure_box(state, fault=FAULT)
