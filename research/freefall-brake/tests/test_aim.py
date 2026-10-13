import json, os, sys
from decimal import Decimal
from fractions import Fraction as F
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ff_aim import aim, flow, energy
from ff_prescan import _solve

CAT = {c["id"]: c for c in json.load(open(os.path.join(ROOT, "results", "li_liao_2019_catalog.json"), encoding="utf-8"))}
MS2 = os.path.join(ROOT, "results", "ms2")


def _start(cid):
    c = CAT[cid]
    return [F(m) for m in c["masses"]], F(c["x"]), F(c["y"]), F(c["T"]) / 2


def test_flow_matches_float_solver():
    masses, x, y, tau = _start("F1(1,0.8,0.8)")
    P, V, _ = flow(masses, x, y, tau)
    ref = _solve([float(m) for m in masses], float(x), float(y), float(tau)).y[:, -1]
    got = np.array([float(v) for v in P + V])
    assert np.max(np.abs(got - ref)) < 1e-10


def test_energy_conserved():
    masses, x, y, tau = _start("F2(1,1,1)")
    P, V, _ = flow(masses, x, y, tau)
    P0 = [F(-1, 2), 0, F(1, 2), 0, x, y]
    assert abs(energy(masses, P, V) - energy(masses, P0, [0] * 6)) < Decimal("1e-35")


def _inside_proven_K(rec_name, res):
    K = json.load(open(os.path.join(MS2, rec_name), encoding="utf-8"))["K"]
    for (lo, hi), v in zip([K[0], K[1], K[-1]], [res["x"], res["y"], res["tau"]]):
        assert F(lo) <= F(v) <= F(hi), (lo, v, hi)


def test_aim_n1_converges_into_proven_box():
    res = aim(*_start("F1(1,0.8,0.8)"))
    assert res["g_best"] < 1e-30, res["hist"]
    _inside_proven_K("ms_n1_f1_1_08_08_k8_r50.json", res)


def test_aim_c2_c3_converge_into_proven_box():
    for cid, rec in (("F2(1,1,1)", "ms_c2_f2_111_k8_r50.json"), ("F1(0.6,0.8,1)", "ms_c3_f1_0608_1_k8_r50.json")):
        res = aim(*_start(cid))
        assert res["g_best"] < 1e-30, (cid, res["hist"])
        _inside_proven_K(rec, res)
