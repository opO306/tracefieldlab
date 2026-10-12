"""R1 검사기 자기시험용 가짜 구현: float_poly

결함: 다항식 부분을 double 부동소수점으로 계산해 점 구간으로 만든다(반올림 오차 무시).
(_fake_r1_core.py의 fault 이름 하나만 고정. 정답 구현과 나머지는 동일.)
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _fake_r1_core as core
from _fake_r1_core import make_harmonic, make_kepler, make_nbody

FAULT = 'float_poly'


def taylor_step(system, X0, h, order):
    return core.taylor_step(system, X0, h, order, fault=FAULT)
