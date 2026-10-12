"""R3c 검사기 자기시험 가짜 구현: no_remainder - 테일러 나머지 Z 생략.  (_fake_r3c_core.py 위에서 fault 하나만 고정)"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _fake_r3c_core as core
from _fake_r3c_core import make_harmonic, make_kepler, make_nbody

FAULT = 'no_remainder'


def lohner_init(system, x0):
    return core.lohner_init(system, x0)


def lohner_step(system, state, h, order):
    return core.lohner_step(system, state, h, order, fault=FAULT)


def enclosure_box(state):
    return core.enclosure_box(state, fault=FAULT)
