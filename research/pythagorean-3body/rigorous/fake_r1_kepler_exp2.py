"""R1 검사기 자기시험용 가짜 구현: kepler_exp2

결함: 케플러 힘 법칙을 r^-3 대신 r^-2로 계산한다(원 궤도 r=1에서는 가속도 같음).
(_fake_r1_core.py의 fault 이름 하나만 고정. 정답 구현과 나머지는 동일.)
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _fake_r1_core as core
from _fake_r1_core import make_harmonic, make_kepler, make_nbody

FAULT = 'kepler_exp2'


def taylor_step(system, X0, h, order):
    return core.taylor_step(system, X0, h, order, fault=FAULT)
