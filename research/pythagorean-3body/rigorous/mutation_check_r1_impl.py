"""taylor_step.py 변이 시험 (구현자 작성): 자체 시험이 너무 쉽게 통과하지 않는지 확인.

taylor_step.py 복사본에 일부러 틀린 곳을 하나씩 넣고(임시 폴더), test_taylor_impl.py의
시험 모음을 그 복사본에 대해 실행한다(첫 실패에서 멈춤). 변이가 "잡혔다" = 시험 실패/오류.
원본 taylor_step.py는 읽기만 한다.

실행:  python mutation_check_r1_impl.py
"""

import importlib.util
import os
import sys
import tempfile
import time
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import test_taylor_impl as T  # noqa: E402

MUTANTS = [
    ("pow_coef_i_plus_k", "* (i + 2 * k)", "* (i + k)"),
    ("pow_first_term", "(1 + 2 * k)", "(2 * k)"),
    ("w0_is_inv_s", "return 1 / (s0 * s0.sqrt())", "return 1 / s0"),
    ("sqconv_no_double", "acc = acc * 2  #", "acc = acc  #"),
    ("remainder_degree_p", "return [ci[order + 1] * Hp1 for ci in c]",
     "return [ci[order] * Hp1 for ci in c]"),
    ("remainder_half", "return [ci[order + 1] * Hp1 for ci in c]",
     "return [ci[order + 1] * Hp1 * Fraction(1, 2) for ci in c]"),
    ("remainder_on_X0", "R = _remainder(system, B, order, h_used)",
     "R = _remainder(system, X0, order, h_used)"),
    ("remainder_h_pow_p2", "Hp1 = DI(h ** (order + 1))", "Hp1 = DI(h ** (order + 2))"),
    ("subset_always_true", "return b.lo_q <= z.lo_q and z.hi_q <= b.hi_q", "return True"),
    ("apriori_returns_X0", "return Z, h_try, iters", "return list(X0), h_try, iters"),
    ("apriori_half_time", "T = DI(Fraction(0), h_try)", "T = DI(Fraction(0), h_try / 2)"),
    ("h_used_not_reduced", "return Z, h_try, iters", "return Z, h, iters"),
    ("kepler_sign", "ax = -(MU * _conv(x, w, k))", "ax = (MU * _conv(x, w, k))"),
    ("nbody_reaction_sign", "acc[2 * b] = acc[2 * b] - M[a] * ux",
     "acc[2 * b] = acc[2 * b] + M[a] * ux"),
    ("nbody_dy_mixup", "dy[p].append(pos[2 * b + 1][k] - pos[2 * a + 1][k])",
     "dy[p].append(pos[2 * b + 1][k] - pos[2 * a][k])"),
    ("harmonic_sign", "v.append(_div_int(-x[k], k + 1))", "v.append(_div_int(x[k], k + 1))"),
    ("harmonic_factorial", "x.append(_div_int(v[k], k + 1))", "x.append(_div_int(v[k], k + 2))"),
    ("poly_drop_top", "acc = ci[order] * Hk[order - 1]", "acc = ci[order] * Hk[order - 1] * 0"),
    ("prec_check_off", "if not _bits_ok(X0, prec):", "if False:"),
]


def main():
    src = open(os.path.join(HERE, "taylor_step.py"), encoding="utf-8").read()
    orig_ts = T.ts
    results = []
    with tempfile.TemporaryDirectory() as td:
        for i, (name, old, new) in enumerate(MUTANTS):
            n = src.count(old)
            if n != 1:
                results.append((name, "BAD_PATTERN(count=%d)" % n, 0.0, ""))
                print("%-22s %-9s %6.1fs  %s" % results[-1], flush=True)
                continue
            path = os.path.join(td, "ts_mut_%d.py" % i)
            with open(path, "w", encoding="utf-8") as fh:
                fh.write(src.replace(old, new))
            spec = importlib.util.spec_from_file_location("ts_mut_%d" % i, path)
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            T.ts = mod
            t0 = time.time()
            suite = unittest.defaultTestLoader.loadTestsFromModule(T)
            res = unittest.TextTestRunner(stream=open(os.devnull, "w"), failfast=True).run(suite)
            dt = time.time() - t0
            killed = bool(res.failures or res.errors)
            first = ""
            if res.failures:
                first = res.failures[0][0].id().split(".", 1)[1]
            elif res.errors:
                first = res.errors[0][0].id().split(".", 1)[1]
            results.append((name, "KILLED" if killed else "SURVIVED", dt, first))
            print("%-22s %-9s %6.1fs  %s" % results[-1], flush=True)
    T.ts = orig_ts
    n_k = sum(r[1] == "KILLED" for r in results)
    print("killed %d / %d" % (n_k, len(results)))


if __name__ == "__main__":
    main()
