"""R3c 검사기 자기시험 가짜 구현: inv_float - [B'^-1] 을 부동소수(double) 역행렬의 점 값으로만 사용(정확한 유리수 역행렬+바깥 반올림 생략).  (_fake_r3c_core.py 위에서 fault 하나만 고정)"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _fake_r3c_core as core
from _fake_r3c_core import make_harmonic, make_kepler, make_nbody

FAULT = 'inv_float'


def lohner_init(system, x0):
    return core.lohner_init(system, x0)


def lohner_step(system, state, h, order):
    return core.lohner_step(system, state, h, order, fault=FAULT)


def enclosure_box(state):
    return core.enclosure_box(state, fault=FAULT)
