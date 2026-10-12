"""R1 검사기 자기시험용 가짜 구현: picard_half_interval

결함: 사전 상자 검사를 [0,h] 대신 [0,h/2] 구간으로 수행한다(그러고도 h 전체를 보고한다).
(_fake_r1_core.py의 fault 이름 하나만 고정. 정답 구현과 나머지는 동일.)
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _fake_r1_core as core
from _fake_r1_core import make_harmonic, make_kepler, make_nbody

FAULT = 'picard_half_interval'


def taylor_step(system, X0, h, order):
    return core.taylor_step(system, X0, h, order, fault=FAULT)
