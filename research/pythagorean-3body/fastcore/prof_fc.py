"""prof_fc.py - Gate FC1 FC-5 보조: 원본/FC system 으로 lohner_step 몇 걸음을 cProfile 로 잰다(측정만).

r3c_run 의 보폭 규칙(J1~J3)을 그대로 쓰기 위해 r3c_run 모듈의 step_size 를 import 해서 쓴다.
사용: python3 -I -B prof_fc.py --system fc|orig [--steps 20] [--prec 128] [--top 25]
"""
import argparse
import cProfile
import os
import pstats
import sys
import time
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fast_system as FS  # noqa: E402
import dyad  # noqa: E402
import taylor_step as ts  # noqa: E402
import lohner_step as ls  # noqa: E402
import r3c_run as R  # noqa: E402


def run(system_kind, steps, prec):
    dyad.set_prec(prec)
    system = FS.make_nbody_fast(R.MASSES) if system_kind == "fc" else ts.make_nbody(R.MASSES)
    state = ls.lohner_init(system, R.INIT)
    tol_base = Fraction(1, 2 ** (prec - 16))
    for _ in range(steps):
        prep = ls.lohner_prepare(system, state, R.ORDER)
        cabs = max(abs(x) for x in state.c)
        tol = tol_base * max(Fraction(1), cabs)
        h = R.step_size(tol, prep.center_coeff_mag)
        while True:
            new_state, info = ls.lohner_step(system, state, h, R.ORDER, diag={}, prep=prep, pivot="width")
            if max(info["width_remainder"]) > tol:
                h = info["h_used"] / 2
                continue
            break
        state = new_state
    return state


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--system", choices=["fc", "orig"], required=True)
    ap.add_argument("--steps", type=int, default=20)
    ap.add_argument("--prec", type=int, default=128)
    ap.add_argument("--top", type=int, default=25)
    ap.add_argument("--noprof", action="store_true", help="cProfile 없이 걸음 시간만")
    a = ap.parse_args()
    if a.noprof:
        t0 = time.perf_counter()
        st = run(a.system, a.steps, a.prec)
        dt = time.perf_counter() - t0
        print("system=%s steps=%d wall=%.3fs per_step=%.4fs t=%s" % (a.system, a.steps, dt, dt / a.steps, float(st.t)))
        return 0
    pr = cProfile.Profile()
    t0 = time.perf_counter()
    pr.enable()
    st = run(a.system, a.steps, a.prec)
    pr.disable()
    print("system=%s steps=%d wall=%.2fs t=%s" % (a.system, a.steps, time.perf_counter() - t0, float(st.t)))
    ps = pstats.Stats(pr).sort_stats("tottime")
    ps.print_stats(a.top)


if __name__ == "__main__":
    sys.exit(main())
