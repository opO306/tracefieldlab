#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""check_fc1_unit_tests.py - check_fc1.py 의 비교기·입력 생성기 단위 시험 (stdlib unittest).

검사기 자신의 핵심 부품(비교기)이 맞는지 확인한다. 구현(fastcore/)과 무관하다.
실행:  python -I fastcore_check/check_fc1_unit_tests.py
"""
import os
import sys
import unittest
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
for p in (HERE, os.path.join(ROOT, "rigorous")):
    if p not in sys.path:
        sys.path.insert(0, p)

import dyad  # noqa: E402
from dyad import DI  # noqa: E402
import lohner_step  # noqa: E402
import taylor_step  # noqa: E402
import check_fc1 as C  # noqa: E402
import fc1_cases as G  # noqa: E402

GD = lohner_step._GD


def kinds(diffs):
    return sorted(set(d[0] for d in diffs))


class CompareDI(unittest.TestCase):
    def setUp(self):
        dyad.set_prec(128)

    def test_equal_di(self):
        a = DI(Fraction(1, 3), Fraction(1, 2))
        b = DI._new(a._l, a._h)
        self.assertEqual(C.compare_result(a, b), [])

    def test_mantissa_differs(self):
        a = DI(Fraction(1, 3), Fraction(1, 2))
        b = DI._new((a._l[0] + 2, a._l[1]), a._h)
        d = C.compare_result(a, b)
        self.assertEqual(kinds(d), ["endpoint_value"])

    def test_exponent_differs(self):
        a = DI(3, 5)
        b = DI._new((a._l[0], a._l[1] + 1), a._h)
        self.assertEqual(kinds(C.compare_result(a, b)), ["endpoint_value"])

    def test_hi_side_differs(self):
        a = DI(Fraction(1, 3), Fraction(1, 2))
        b = DI._new(a._l, (a._h[0], a._h[1] - 1))
        d = C.compare_result(a, b)
        self.assertEqual(len(d), 1)
        self.assertIn("_h", d[0][1])

    def test_noncanonical_same_value_is_flagged(self):
        a = DI(3, 3)  # (3, 0)
        b = DI._new((6, a._l[1] - 1), a._h)  # 같은 값, 표준형 아님
        d = C.compare_result(a, b)
        self.assertEqual(kinds(d), ["endpoint_repr"])

    def test_bool_or_float_endpoint_is_flagged(self):
        a = DI(1, 1)  # (1, 0)
        b = DI._new((True, 0), a._h)
        self.assertEqual(kinds(C.compare_result(a, b)), ["endpoint_form"])
        c = DI._new((1.0, 0), a._h)
        self.assertEqual(kinds(C.compare_result(a, c)), ["endpoint_form"])

    def test_list_instead_of_tuple_endpoint(self):
        a = DI(1, 1)
        b = DI._new([1, 0], a._h)
        self.assertEqual(kinds(C.compare_result(a, b)), ["endpoint_form"])

    def test_wrong_class(self):
        a = DI(1, 2)

        class Fake(object):
            _l = a._l
            _h = a._h
        self.assertEqual(kinds(C.compare_result(a, Fake())), ["type"])

    def test_subclass_rejected(self):
        a = DI(1, 2)
        sub = type("Sub", (DI,), {"__slots__": ()})
        b = sub._new(a._l, a._h)
        self.assertEqual(kinds(C.compare_result(a, b)), ["type"])


class CompareStructure(unittest.TestCase):
    def setUp(self):
        dyad.set_prec(128)

    def test_list_length(self):
        a = [[DI(1)], [DI(2)]]
        b = [[DI(1)]]
        self.assertEqual(kinds(C.compare_result(a, b)), ["length"])

    def test_inner_list_length(self):
        a = [[DI(1), DI(2)]]
        b = [[DI(1)]]
        self.assertEqual(kinds(C.compare_result(a, b)), ["length"])

    def test_container_type_is_form(self):
        a = [[DI(1)]]
        b = [(DI(1),)]
        self.assertEqual(kinds(C.compare_result(a, b)), ["form"])

    def test_tuple_outer_is_form(self):
        a = [[DI(1)]]
        b = ([DI(1)],)
        self.assertEqual(kinds(C.compare_result(a, b)), ["form"])

    def test_none_vs_di0_in_gradient(self):
        a = GD(DI(1), [None, DI(1)])
        b = GD(DI(1), [DI(0), DI(1)])
        d = C.compare_result(a, b)
        self.assertEqual(kinds(d), ["none_position"])
        d2 = C.compare_result(b, a)
        self.assertEqual(kinds(d2), ["none_position"])

    def test_gd_equal(self):
        a = GD(DI(1), [None, DI(1)])
        b = GD(DI._new(DI(1)._l, DI(1)._h), [None, DI(1)])
        self.assertEqual(C.compare_result(a, b), [])

    def test_gd_value_differs(self):
        a = GD(DI(1), [None, DI(1)])
        b = GD(DI(2), [None, DI(1)])
        self.assertEqual(kinds(C.compare_result(a, b)), ["endpoint_value"])

    def test_gd_gradient_length(self):
        a = GD(DI(1), [None, DI(1)])
        b = GD(DI(1), [None])
        self.assertEqual(kinds(C.compare_result(a, b)), ["length"])

    def test_gd_vs_di(self):
        a = GD(DI(1), [None])
        self.assertEqual(kinds(C.compare_result(a, DI(1))), ["type"])
        self.assertEqual(kinds(C.compare_result(DI(1), a)), ["type"])

    def test_gd_gradient_tuple_is_form(self):
        a = GD(DI(1), [None, DI(1)])
        b = GD(DI(1), (None, DI(1)))
        self.assertEqual(kinds(C.compare_result(a, b)), ["form"])

    def test_none_leaf_other_object(self):
        # 오라클이 None 인 자리에 다른 것이 오면 none_position
        a = [None]
        b = [DI(0)]
        self.assertEqual(kinds(C.compare_result(a, b)), ["none_position"])

    def test_diff_count_and_cap(self):
        dyad.set_prec(64)
        a = [[DI(i) for i in range(1, 11)]]
        b = [[DI(i + 100) for i in range(1, 11)]]
        d = C.compare_result(a, b, max_details=3)
        self.assertEqual(len(d), 3)


class CompareOutcome(unittest.TestCase):
    def test_same_exception(self):
        o = ("exc", taylor_step._Singular())
        t = ("exc", taylor_step._Singular())
        self.assertEqual(C.compare_outcome(o, t), [])

    def test_exception_type_differs(self):
        o = ("exc", taylor_step._Singular())
        t = ("exc", ZeroDivisionError("x"))
        self.assertEqual(kinds(C.compare_outcome(o, t)), ["exception_type"])

    def test_exception_subclass_rejected(self):
        class Sub(taylor_step._Singular):
            pass
        o = ("exc", taylor_step._Singular())
        t = ("exc", Sub())
        self.assertEqual(kinds(C.compare_outcome(o, t)), ["exception_type"])

    def test_exception_args_differ(self):
        o = ("exc", ValueError("a"))
        t = ("exc", ValueError("b"))
        self.assertEqual(kinds(C.compare_outcome(o, t)), ["exception_args"])

    def test_result_vs_exception(self):
        dyad.set_prec(64)
        o = ("ok", [[DI(1)]])
        t = ("exc", taylor_step._Singular())
        self.assertEqual(kinds(C.compare_outcome(o, t)), ["result_vs_exception"])
        self.assertEqual(kinds(C.compare_outcome(t, o)), ["result_vs_exception"])

    def test_ok_ok_delegates(self):
        dyad.set_prec(64)
        o = ("ok", [[DI(1)]])
        t = ("ok", [[DI(2)]])
        self.assertEqual(kinds(C.compare_outcome(o, t)), ["endpoint_value"])


class LoaderTests(unittest.TestCase):
    def _mk(self, tmp, with_init):
        pkg = os.path.join(tmp, "fcpkg_%d" % (1 if with_init else 0))
        os.makedirs(pkg)
        open(os.path.join(pkg, "_helper.py"), "w").write("VALUE = 7\n")
        body = ("import taylor_step\n%s\n"
                "def make_nbody_fast(masses):\n    return taylor_step.make_nbody(masses)\n")
        rel = "from . import _helper\n" if with_init else "import _helper\n"
        open(os.path.join(pkg, "fast_system.py"), "w").write(body % rel)
        if with_init:
            open(os.path.join(pkg, "__init__.py"), "w").write("")
        return os.path.join(pkg, "fast_system.py")

    def test_plain_directory_target_with_sibling_import(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            t = C.load_target(self._mk(tmp, False), 90)
            self.assertTrue(callable(t.make))

    def test_package_directory_target_with_relative_import(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            t = C.load_target(self._mk(tmp, True), 91)
            self.assertTrue(callable(t.make))

    def test_target_without_api_is_rejected(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            p = os.path.join(tmp, "noapi.py")
            open(p, "w").write("X = 1\n")
            with self.assertRaises(RuntimeError):
                C.load_target(p, 92)


class SnapshotTests(unittest.TestCase):
    def setUp(self):
        dyad.set_prec(64)

    def test_snapshot_detects_list_mutation(self):
        x = [DI(1), DI(2)]
        s = C.snapshot(x)
        x[0], x[1] = x[1], x[0]
        self.assertNotEqual(s, C.snapshot(x))

    def test_snapshot_detects_gd_gradient_mutation(self):
        g = GD(DI(1), [None, DI(1)])
        x = [g]
        s = C.snapshot(x)
        g.g[0] = DI(0)
        self.assertNotEqual(s, C.snapshot(x))

    def test_snapshot_stable(self):
        x = [GD(DI(1), [None, DI(1)]), DI(3)]
        self.assertEqual(C.snapshot(x), C.snapshot(list(x)))


class GeneratorTests(unittest.TestCase):
    def test_rng_below_uniform_range(self):
        r = G.Rng(1)
        vals = [r.below(7) for _ in range(2000)]
        self.assertEqual(sorted(set(vals)), list(range(7)))

    def test_determinism(self):
        a = G.generate("smoke", 5)
        b = G.generate("smoke", 5)
        self.assertEqual(a, b)
        c = G.generate("smoke", 6)
        self.assertNotEqual(a, c)

    def test_boxes_deterministic_and_valid(self):
        specs = G.generate("smoke", 3)
        self.assertGreater(len(specs), 100)
        for s in specs:
            p1 = G.gen_box(s)
            p2 = G.gen_box(s)
            self.assertEqual(p1["pairs"], p2["pairs"])
            self.assertEqual(len(p1["pairs"]), 4 * s["n"])
            for lo, hi in p1["pairs"]:
                self.assertLessEqual(lo, hi, s)
            self.assertEqual(len(p1["shifts"]), 4 * s["n"])

    def test_categories_covered_in_smoke(self):
        specs = G.generate("smoke", 3)
        cats = set(s["cat"] for s in specs)
        for c in G.ALL_CATS_REQUIRED:
            self.assertIn(c, cats, c)

    def test_masses_parse(self):
        for s in G.generate("smoke", 3):
            ms = [Fraction(x) for x in s["masses"]]
            self.assertEqual(len(ms), s["n"])

    def test_plan_profile_meets_plan_grid(self):
        specs = G.generate("plan", 1)
        # 계획 4절 FC-1/FC-2: 정밀도 x 차수 x 질량 x 체수 조합마다 >= 50 상자 (DI·GD 각각)
        from collections import Counter
        cnt = Counter()
        for s in specs:
            if s["blk"] == "A":
                cnt[(s["path"], s["prec"], s["K"], s["ms"], s["n"])] += 1
        want = [(path, p, K, ms, n)
                for path in ("DI", "GD")
                for p in (96, 112, 128, 144, 160, 192, 256)
                for K in (1, 5, 20, 30, 40)
                for ms in ("m345", "rand")
                for n in (2, 3)]
        self.assertEqual(len(want), 280)
        for k in want:
            self.assertGreaterEqual(cnt[k], 50, k)

    def test_pyth_pairwise_squared_distances_are_perfect_squares(self):
        from math import isqrt
        seen = 0
        for sd in range(40):
            spec = {"blk": "B", "tier": "plan", "cat": "pyth", "path": "DI", "n": 3, "ms": "m345",
                    "masses": ["3", "4", "5"], "prec": 128, "K": 1, "seed": G.sub_seed("t", sd),
                    "reuse": False, "mkprec": 128, "i": 0}
            pairs = G.gen_box(spec)["pairs"]
            if not all(lo == hi for lo, hi in pairs[:6]):
                continue        # 폭이 있는 상자는 비대칭이라 중점이 시작 배치가 아님 — 점 구간만 확인
            c = [lo for lo, hi in pairs]
            for a in range(3):
                for b in range(a + 1, 3):
                    q = (c[2 * a] - c[2 * b]) ** 2 + (c[2 * a + 1] - c[2 * b + 1]) ** 2
                    self.assertEqual(isqrt(q.numerator) ** 2, q.numerator)
                    self.assertEqual(isqrt(q.denominator) ** 2, q.denominator)
                    seen += 1
        self.assertGreaterEqual(seen, 30)

    def test_mass_pool_repeats_inexact_masses(self):
        from collections import Counter
        specs = G.generate("selftest", 1)
        cnt = Counter(tuple(s["masses"]) for s in specs if s["ms"] == "rand")
        repeated = [m for m, c in cnt.items() if c >= 3 and any("/" in x and not
                    (Fraction(x).denominator & (Fraction(x).denominator - 1)) == 0 for x in m)]
        self.assertGreaterEqual(len(repeated), 3)

    def test_exp_scaled_stays_within_oracle_range(self):
        for sd in range(200):
            spec = {"blk": "B", "tier": "plan", "cat": "exp_scaled", "path": "DI", "n": 3, "ms": "m345",
                    "masses": ["3", "4", "5"], "prec": 64, "K": 1, "seed": G.sub_seed("e", sd),
                    "reuse": False, "mkprec": 64, "i": 0}
            bx = G.gen_box(spec)
            self.assertLessEqual(abs(bx["E"]), 1 << 27)
            self.assertEqual(bx["E"] % 2, 0)

    def test_di_and_gd_share_boxes_in_block_a(self):
        specs = G.generate("smoke", 2)
        by = {}
        for s in specs:
            if s["blk"] == "A":
                by.setdefault((s["prec"], s["K"], s["ms"], s["n"], s["seed"]), set()).add(s["path"])
        self.assertTrue(by)
        for k, v in by.items():
            self.assertEqual(v, {"DI", "GD"}, k)


if __name__ == "__main__":
    unittest.main(verbosity=2)
