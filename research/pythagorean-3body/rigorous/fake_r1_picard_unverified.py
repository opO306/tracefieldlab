"""R1 검사기 자기시험용 가짜 구현: picard_unverified

결함: 포함 검사 없이 한 번짜리 상자 X0 ∪ (X0+[0,h]F(X0))를 B로 쓴다.
(_fake_r1_core.py의 fault 이름 하나만 고정. 정답 구현과 나머지는 동일.)
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _fake_r1_core as core
from _fake_r1_core import make_harmonic, make_kepler, make_nbody

FAULT = 'picard_unverified'


def taylor_step(system, X0, h, order):
    return core.taylor_step(system, X0, h, order, fault=FAULT)
