#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""check_dyad.py -- Gate R0 독립 검사기 (표준 라이브러리만).

계획 GATE_R0_PLAN_KO.md 2절 API(DI, set_prec, get_prec ...)만 근거로 만들었다.
검사 대상 모듈의 내부 코드는 보지 않고, 정확한 유리수(fractions.Fraction)를 기준값으로 쓴다.
실행기(dyad.py)의 계산 코드를 재사용하지 않는다.

사용법
  python check_dyad.py MODULE.py [--seed N] [--trials N] [--out report.json]
  python check_dyad.py --selftest        # fake_dyad_*.py 로 검사기 자체를 시험

검사 ID
  R0-S1  무작위 포함 시험 (정확한 값이 출력 구간 밖이면 위반)
  R0-S2  경계 사례 목록(코드 안 S2_* 상수) + 생성자의 "정확히 표현되는 값은 그대로"
  R0-S3  sqrt 포함을 제곱 비교로 정확히 확인
  R0-S4  끝점 가수 비트 수 <= prec (분모가 2의 거듭제곱인지, 가수 비트는 몇인지 직접 계산)
  R0-API 예외 계약, 타입/불변식, contains/width_q/hull/intersect 의미, 입력 불변성
위반이 0이어도 "수학적 증명"이 아니라 "이 시험 범위 안에서 위반 없음"이다.
"""
import argparse
import hashlib
import importlib.util
import json
import operator
import os
import random
import subprocess
import sys
import time
from fractions import Fraction

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)  # 2^10000 같은 큰 정수의 십진 변환 제한 해제

CHECKER_VERSION = "1.1"
PRECS = (24, 53, 128, 300)
REQUIRED_S1_TRIALS = 200000
IDS = ("R0-API", "R0-S1", "R0-S2", "R0-S3", "R0-S4")


# ----------------------------------------------------------------------------
# 정확한 유리수 도우미 (모두 정수 연산, 부동소수점 없음)
# ----------------------------------------------------------------------------
def mant_bits(q):
    """q가 m*2^e 꼴(분모가 2의 거듭제곱)이면 홀수 가수의 비트 수, 아니면 None."""
    n, d = q.numerator, q.denominator
    if d & (d - 1):
        return None
    if n == 0:
        return 0
    n = abs(n)
    if d == 1:
        n >>= (n & -n).bit_length() - 1  # 끝의 0 비트 제거(지수로 이동)
    return n.bit_length()


def floor_log2_pos(q):
    n, d = q.numerator, q.denominator
    e = n.bit_length() - d.bit_length()
    if e >= 0:
        lhs, rhs = n, d << e
    else:
        lhs, rhs = n << (-e), d
    if lhs < rhs:
        e -= 1
    return e


def grid_round(q, p, up):
    """q를 가수 p비트 격자로 바깥(up=True: 올림, False: 내림) 반올림한 정확한 기준값."""
    if q == 0:
        return Fraction(0)
    k = floor_log2_pos(abs(q)) - p + 1
    if k >= 0:
        num, den = q.numerator, q.denominator << k
    else:
        num, den = q.numerator << (-k), q.denominator
    r = -((-num) // den) if up else num // den
    return Fraction(r) * (Fraction(2) ** k)


def qsize(q):
    return q.numerator.bit_length() + q.denominator.bit_length()


def istr(n):
    """정수를 문자열로. 아주 큰 정수는 16진(0x...)으로(십진 변환 비용·제한 회피)."""
    return str(n) if abs(n).bit_length() <= 3000 else hex(n)


def fmt(q):
    if isinstance(q, bool):
        return str(q)
    if isinstance(q, int):
        q = Fraction(q)
    if not isinstance(q, Fraction):
        return repr(q)
    n, d = q.numerator, q.denominator
    if d == 1:
        if abs(n).bit_length() <= 130:
            return str(n)
        e = (abs(n) & -abs(n)).bit_length() - 1
        m = n >> e
        if abs(m).bit_length() <= 400:
            return "%s*2^%d" % (istr(m), e)
        return "int(bits=%d,hex_head=%s...)" % (abs(n).bit_length(), hex(abs(n) >> (abs(n).bit_length() - 32)))
    if d & (d - 1) == 0:
        k = d.bit_length() - 1
        if abs(n).bit_length() <= 400:
            return "%s*2^-%d" % (istr(n), k)
    if abs(n).bit_length() <= 400 and d.bit_length() <= 400:
        return "%s/%s" % (istr(n), istr(d))
    return "rational(num_bits=%d,den_bits=%d)" % (abs(n).bit_length(), d.bit_length())


def walk_q(v):
    if isinstance(v, Fraction):
        yield v
    elif isinstance(v, int) and not isinstance(v, bool):
        yield Fraction(v)
    elif isinstance(v, (list, tuple)):
        for x in v:
            yield from walk_q(x)
    elif isinstance(v, dict):
        for x in v.values():
            yield from walk_q(x)


def enc(v):
    if isinstance(v, (Fraction, int)) and not isinstance(v, bool):
        return fmt(v)
    if isinstance(v, (list, tuple)):
        return [enc(x) for x in v]
    if isinstance(v, dict):
        return {str(k): enc(x) for k, x in v.items()}
    return v if isinstance(v, (str, bool, type(None))) else repr(v)


def raw(v):
    if isinstance(v, (Fraction, int)) and not isinstance(v, bool):
        q = Fraction(v)
        return [istr(q.numerator), istr(q.denominator)]
    if isinstance(v, (list, tuple)):
        return [raw(x) for x in v]
    if isinstance(v, dict):
        return {str(k): raw(x) for k, x in v.items()}
    return v if isinstance(v, (str, bool, type(None))) else repr(v)


# ----------------------------------------------------------------------------
# 정확한 기준 연산 (구간 끝점은 모두 정확한 Fraction)
# ----------------------------------------------------------------------------
def ex_range(op, x, y=None, n=None):
    """입력 구간 x=(lo,hi), y=(lo,hi) 위에서 op의 정확한 값 집합의 닫힌 범위."""
    if op == "add":
        return (x[0] + y[0], x[1] + y[1])
    if op == "sub":
        return (x[0] - y[1], x[1] - y[0])
    if op == "mul":
        ps = [a * b for a in x for b in y]
        return (min(ps), max(ps))
    if op == "div":  # y는 0을 포함하지 않는다고 가정
        qs = [a / b for a in x for b in y]
        return (min(qs), max(qs))
    if op == "neg":
        return (-x[1], -x[0])
    if op == "pow":
        if n == 0:
            return (Fraction(1), Fraction(1))
        a, b = x[0] ** n, x[1] ** n
        if n % 2 == 1:
            return (a, b)
        if x[0] <= 0 <= x[1]:
            return (Fraction(0), max(a, b))
        return (min(a, b), max(a, b))
    raise ValueError(op)


def ex_val(op, a, b=None, n=None):
    if op == "add":
        return a + b
    if op == "sub":
        return a - b
    if op == "mul":
        return a * b
    if op == "div":
        return a / b
    if op == "neg":
        return -a
    if op == "pow":
        return a ** n
    raise ValueError(op)


PYOPS = {"add": operator.add, "sub": operator.sub, "mul": operator.mul, "div": operator.truediv}


# ----------------------------------------------------------------------------
# 결과 기록
# ----------------------------------------------------------------------------
class Report:
    def __init__(self):
        self.cases = {i: 0 for i in IDS}
        self.by = {i: {} for i in IDS}
        self.warn = {}
        self.point_checks = 0
        self.trials = {}          # 이름 -> 평가 횟수(전 단계)
        self.rand_trials = 0      # 무작위 단계 연산 평가 횟수
        self.t1 = {}              # (op,prec) -> [ratio,...]

    def case(self, cid, n=1):
        self.cases[cid] += n

    def trial(self, name, rand):
        self.trials[name] = self.trials.get(name, 0) + 1
        if rand:
            self.rand_trials += 1

    def viol(self, cid, op, prec, inputs, outputs, detail):
        size = sum(qsize(q) for q in walk_q(inputs)) + prec // 8
        ex = {
            "check": cid, "op": op, "prec": prec, "detail": detail,
            "inputs": enc(inputs), "outputs": enc(outputs), "_size": size,
            "_raw": [raw(inputs), raw(outputs)],
        }
        b = self.by[cid].setdefault(op, {"count": 0, "first": None, "smallest": []})
        b["count"] += 1
        if b["first"] is None:
            b["first"] = ex
        sm = b["smallest"]
        if len(sm) < 3 or size < sm[-1]["_size"]:
            sm.append(ex)
            sm.sort(key=lambda e: e["_size"])
            del sm[3:]

    def warning(self, name, op, prec, inputs, outputs, detail):
        w = self.warn.setdefault(name, {"count": 0, "first": None})
        w["count"] += 1
        if w["first"] is None:
            w["first"] = {"op": op, "prec": prec, "detail": detail,
                          "inputs": enc(inputs), "outputs": enc(outputs)}

    def t1_add(self, op, p, ratio):
        self.t1.setdefault("%s@%d" % (op, p), []).append(ratio)


def strip_private(ex):
    return {k: v for k, v in ex.items() if not k.startswith("_")}


# ----------------------------------------------------------------------------
# 무작위 입력 생성
# ----------------------------------------------------------------------------
def rand_val(rng, p):
    r = rng.random()
    if r < 0.22:
        return Fraction(rng.randint(-20, 20))
    if r < 0.34:
        return Fraction(rng.randint(-64, 64), 1 << rng.randint(0, 8))
    if r < 0.58:
        b = rng.choice((4, 20, p, p + 10, 2 * p))
        return Fraction(rng.randint(-(1 << b), 1 << b), rng.randint(1, 1 << b))
    if r < 0.72:
        b = rng.choice((3, p, p + 1, p + 5))
        n = rng.randint(1, 1 << b) * rng.choice((-1, 1))
        e = rng.randint(-2200, 2200)
        return Fraction(n << e) if e >= 0 else Fraction(n, 1 << (-e))
    if r < 0.86:
        base = (1 << rng.choice((p, p + 1, p + 2))) + rng.randint(-3, 3)
        base *= rng.choice((-1, 1))
        e = rng.choice((0, rng.randint(-2100, 2100), -p, -p - 3))
        return Fraction(base << e) if e >= 0 else Fraction(base, 1 << (-e))
    if r < 0.93:
        return Fraction(rng.randint(-100, 100), rng.randint(1, 10 ** rng.randint(1, 40)))
    return Fraction(0)


def rand_iv(rng, p):
    """요청 구간 (lo, hi) 정확한 Fraction, lo <= hi."""
    r = rng.random()
    v = rand_val(rng, p)
    if r < 0.12:
        return v, v
    if r < 0.30:
        scale = max(abs(v), Fraction(1, 1 << 8))
        delta = scale * Fraction(1, 1 << rng.randint(1, p + 20))
        return v, v + delta
    if r < 0.45:
        f = Fraction(rng.randint(1, 1000), rng.randint(1, 1000))
        return v, v + abs(v) * f
    w = rand_val(rng, p)
    return (v, w) if v <= w else (w, v)


def pts(rng, lo, hi, with_zero=False):
    out = [lo]
    if hi != lo:
        out.append(hi)
        out.append((lo + hi) / 2)
        out.append(lo + (hi - lo) * Fraction(rng.getrandbits(24), 1 << 24))
        if with_zero and lo < 0 < hi:
            out.append(Fraction(0))
    return out


def as_scalar(rng, q):
    return int(q) if q.denominator == 1 and rng.random() < 0.5 else q


# ----------------------------------------------------------------------------
# 검사기 본체
# ----------------------------------------------------------------------------
class Tester:
    def __init__(self, mod, rep, rng):
        self.m = mod
        self.DI = mod.DI
        self.rep = rep
        self.rng = rng
        self.phase = "R0-S1"
        self.rand = True

    # -- 보조 --
    def set_phase(self, tag):
        self.phase = tag
        self.rand = tag == "R0-S1"

    def is_di(self, o):
        return hasattr(o, "lo_q") and hasattr(o, "hi_q")

    def validate(self, d, op, p, ctx, bits=True):
        rep = self.rep
        try:
            lo, hi = d.lo_q, d.hi_q
        except Exception as e:  # noqa
            rep.viol("R0-API", op, p, ctx, None, "lo_q/hi_q 읽기 실패: %s: %s" % (type(e).__name__, e))
            return None
        if not (isinstance(lo, Fraction) and isinstance(hi, Fraction)):
            rep.viol("R0-API", op, p, ctx, [repr(type(lo)), repr(type(hi))], "lo_q/hi_q가 Fraction이 아님")
            return None
        rep.case("R0-API")
        if lo > hi:
            rep.viol("R0-API", op, p, ctx, [lo, hi], "불변식 lo_q <= hi_q 위반")
        if bits:
            rep.case("R0-S4")
            for name, q in (("lo", lo), ("hi", hi)):
                mb = mant_bits(q)
                if mb is None:
                    rep.viol("R0-S4", op, p, ctx, [lo, hi], "%s_q 분모가 2의 거듭제곱이 아님: %s" % (name, fmt(q)))
                elif mb > p:
                    rep.viol("R0-S4", op, p, ctx, [lo, hi],
                             "%s_q 가수 %d비트 > prec %d: %s" % (name, mb, p, fmt(q)))
        return lo, hi

    def snapshot(self, objs):
        return [(o.lo_q, o.hi_q) for o in objs if self.is_di(o)]

    def check_mutation(self, op, p, objs, snap, ctx):
        self.rep.case("R0-API")
        if self.snapshot(objs) != snap:
            self.rep.viol("R0-API", op, p, ctx, None, "연산이 입력 DI의 끝점을 바꿈")

    # -- 생성자 --
    def mk(self, lo, hi, p):
        """요청 구간 [lo,hi](정확한 Fraction)로 DI 생성 + 생성자 계약 검사. 실패하면 None."""
        rep, rng = self.rep, self.rng
        point = lo == hi and rng.random() < 0.5
        ctx = {"req_lo": lo, "req_hi": hi}
        try:
            d = self.DI(as_scalar(rng, lo)) if point else self.DI(as_scalar(rng, lo), as_scalar(rng, hi))
        except Exception as e:  # noqa
            rep.viol("R0-API", "ctor", p, ctx, None, "생성자 예외 %s: %s" % (type(e).__name__, e))
            return None
        r = self.validate(d, "ctor", p, ctx)
        if r is None:
            return None
        L, H = r
        rep.case(self.phase)
        if not (L <= lo and hi <= H):
            rep.viol(self.phase, "ctor", p, ctx, [L, H], "생성자 구간이 요청 구간을 포함하지 않음")
        # 정확히 표현되는 값은 그대로(바깥쪽 반올림 = 격자 내림/올림이므로 고정점)
        rep.case("R0-S2")
        for name, req, got in (("lo", lo, L), ("hi", hi, H)):
            mb = mant_bits(req)
            if mb is not None and mb <= p and got != req:
                rep.viol("R0-S2", "ctor_exact", p, ctx, [L, H],
                         "정확히 표현되는 %s=%s가 %s로 바뀜" % (name, fmt(req), fmt(got)))
        # 경고: 가장 촘촘한 바깥쪽 반올림보다 느슨한가
        if L < grid_round(lo, p, False) or H > grid_round(hi, p, True):
            rep.warning("ctor_not_tightest", "ctor", p, ctx, [L, H], "격자 내림/올림보다 넓은 끝점")
        return d

    # -- 포함 판정 공통 --
    def contain(self, name, p, inputs, out, erange, pointvals):
        rep = self.rep
        lo, hi = out
        bad = []
        for v in pointvals:
            rep.point_checks += 1
            if not (lo <= v <= hi):
                bad.append(v)
        rng_bad = not (lo <= erange[0] and erange[1] <= hi)
        rep.case(self.phase)
        if bad or rng_bad:
            rep.viol(self.phase, name, p, inputs,
                     {"out": [lo, hi], "exact_range": list(erange), "failing_points": bad[:2]},
                     "정확한 값이 출력 구간 밖(점 위반 %d개, 정확한 범위 %s)" % (len(bad), "위반" if rng_bad else "포함"))
            return False
        return True

    def t1(self, name, p, erange, out):
        """T1 조임 정도 = 출력 폭 / max(정확한 폭, 2^-prec*|값|).  (계획 3절: 폭이 0이면 2^-prec*|값| 기준.
        폭이 반올림 단위보다 작은 구간까지 같은 기준을 쓰도록 max로 확장한 것은 검사기의 해석이다.)"""
        if not self.rand:
            return
        val = max(abs(erange[0]), abs(erange[1]))
        base = max(erange[1] - erange[0], val / (1 << p))
        if base == 0:
            return
        try:
            self.rep.t1_add(name, p, float((out[1] - out[0]) / base))
        except OverflowError:
            self.rep.t1_add(name, p, float("inf"))

    # -- 이항 연산 --
    def t_bin(self, op, p, x, y, xiv, yiv, kind, mixed=False):
        """kind: DD(DI,DI) DS(DI,스칼라) SD(스칼라,DI)."""
        rep = self.rep
        name = op if kind == "DD" else "%s_%s" % (op, kind)
        rep.trial(name, self.rand)
        ctx = {"x": list(xiv), "y": list(yiv)}
        objs = [x, y]
        snap = self.snapshot(objs)
        expect_zero = op == "div" and yiv[0] <= 0 <= yiv[1]
        try:
            out = PYOPS[op](x, y)
        except ZeroDivisionError:
            rep.case("R0-API")
            if not expect_zero:
                rep.viol("R0-API", name, p, ctx, None, "0을 포함하지 않는 제수에서 ZeroDivisionError")
            return
        except TypeError as e:
            if kind == "SD":
                rep.warning("reflected_op_unsupported", name, p, ctx, None,
                            "스칼라가 왼쪽인 연산 미지원(%s) -- 계약 문구가 모호해 경고로만 기록" % e)
            else:
                rep.viol("R0-API", name, p, ctx, None, "TypeError: %s" % e)
            return
        except Exception as e:  # noqa
            rep.viol("R0-API" if expect_zero else self.phase, name, p, ctx, None,
                     "예외 %s: %s (기대: %s)" % (type(e).__name__, e, "ZeroDivisionError" if expect_zero else "정상 결과"))
            return
        if expect_zero:
            rep.case("R0-API")
            rep.viol("R0-API", name, p, ctx, repr(out), "0을 포함하는 제수인데 ZeroDivisionError가 나오지 않음")
            return
        r = self.validate(out, name, p, ctx, bits=True)
        self.check_mutation(name, p, objs, snap, ctx)
        if r is None:
            return
        er = ex_range(op, xiv, yiv)
        pv = [ex_val(op, a, b) for a in pts(self.rng, *xiv) for b in pts(self.rng, *yiv)]
        if self.contain(name, p, ctx, r, er, pv):
            self.t1(name, p, er, r)

    # -- 단항 연산 --
    def t_unary(self, op, p, x, xiv, n=None, bits=True):
        rep = self.rep
        name = op if op != "pow" else "pow%d" % n
        rep.trial(name, self.rand)
        ctx = {"x": list(xiv)}
        snap = self.snapshot([x])
        try:
            if op == "neg":
                out = -x
            elif op == "sqr":
                out = x.sqr()
            else:
                out = x ** n
        except Exception as e:  # noqa
            rep.viol(self.phase, name, p, ctx, None, "예외 %s: %s" % (type(e).__name__, e))
            return
        r = self.validate(out, name, p, ctx, bits=bits)
        self.check_mutation(name, p, [x], snap, ctx)
        if r is None:
            return
        if op == "sqr":
            er = ex_range("pow", xiv, n=2)
            pv = [a * a for a in pts(self.rng, *xiv, with_zero=True)]
            rep.case("R0-API")
            if xiv[0] <= 0 <= xiv[1] and r[0] != 0:
                rep.viol("R0-API", name, p, ctx, list(r), "0을 포함하는 입력인데 sqr 하한이 0이 아님(계약: 하한 0)")
            elif r[0] < 0 and not (xiv[0] <= 0 <= xiv[1]):
                rep.warning("sqr_negative_lower", name, p, ctx, list(r), "0을 포함하지 않는데 음수 하한(건전하지만 느슨)")
        elif op == "pow":
            er = ex_range("pow", xiv, n=n)
            pv = [a ** n for a in pts(self.rng, *xiv, with_zero=True)]
        else:
            er = ex_range("neg", xiv)
            pv = [-a for a in pts(self.rng, *xiv)]
        if self.contain(name, p, ctx, r, er, pv) and op in ("sqr", "pow"):
            self.t1(name, p, er, r)

    # -- sqrt --
    def t_sqrt(self, p, x, xiv):
        rep = self.rep
        rep.trial("sqrt", self.rand)
        lo, hi = xiv
        ctx = {"x": [lo, hi]}
        snap = self.snapshot([x])
        expect_err = lo < 0
        try:
            out = x.sqrt()
        except ValueError:
            rep.case("R0-API")
            if not expect_err:
                rep.viol("R0-API", "sqrt", p, ctx, None, "하한 >= 0인데 ValueError")
            return
        except Exception as e:  # noqa
            rep.viol("R0-API" if expect_err else "R0-S3", "sqrt", p, ctx, None,
                     "예외 %s: %s" % (type(e).__name__, e))
            return
        if expect_err:
            rep.case("R0-API")
            rep.viol("R0-API", "sqrt", p, ctx, repr(out), "하한 < 0인데 ValueError가 나오지 않음")
            return
        r = self.validate(out, "sqrt", p, ctx, bits=True)
        self.check_mutation("sqrt", p, [x], snap, ctx)
        if r is None:
            return
        a, b = r
        rep.case("R0-S3")
        fails = []
        if a < 0:
            fails.append("a<0")
        elif a * a > lo:
            fails.append("a^2 > inf(입력)")
        if b < 0 or b * b < hi:
            fails.append("b^2 < sup(입력)")
        if fails:
            rep.viol("R0-S3", "sqrt", p, ctx, {"out": [a, b], "a^2": a * a, "b^2": b * b},
                     "sqrt 포함 위반: " + ", ".join(fails))
        rep.case(self.phase)
        bad = []
        for v in pts(self.rng, lo, hi):
            rep.point_checks += 1
            if not ((a <= 0 or a * a <= v) and (b >= 0 and v <= b * b)):
                bad.append(v)
        if bad:
            rep.viol(self.phase, "sqrt", p, ctx, {"out": [a, b], "failing_points": bad[:2]},
                     "sqrt(x)가 출력 구간 밖(제곱 비교)")

    # -- hull / intersect --
    def t_hull(self, p, A, B, aiv, biv, bits=False):
        rep = self.rep
        rep.trial("hull", self.rand)
        ctx = {"a": list(aiv), "b": list(biv)}
        snap = self.snapshot([A, B])
        try:
            out = A.hull(B)
        except Exception as e:  # noqa
            rep.viol(self.phase, "hull", p, ctx, None, "예외 %s: %s" % (type(e).__name__, e))
            return
        r = self.validate(out, "hull", p, ctx, bits=bits)
        self.check_mutation("hull", p, [A, B], snap, ctx)
        if r is None:
            return
        er = (min(aiv[0], biv[0]), max(aiv[1], biv[1]))
        pv = list(aiv) + list(biv)
        self.contain("hull", p, ctx, r, er, pv)
        if r != er:
            rep.warning("hull_not_exact", "hull", p, ctx, list(r), "hull이 정확한 최소/최대 끝점이 아님")

    def t_intersect(self, p, A, B, aiv, biv, bits=False):
        rep = self.rep
        rep.trial("intersect", self.rand)
        ctx = {"a": list(aiv), "b": list(biv)}
        snap = self.snapshot([A, B])
        lo, hi = max(aiv[0], biv[0]), min(aiv[1], biv[1])
        empty = lo > hi
        try:
            out = A.intersect(B)
        except Exception as e:  # noqa
            rep.viol(self.phase, "intersect", p, ctx, None, "예외 %s: %s" % (type(e).__name__, e))
            return
        self.check_mutation("intersect", p, [A, B], snap, ctx)
        rep.case("R0-API")
        if empty:
            if out is not None:
                rep.viol("R0-API", "intersect", p, ctx, repr(out), "교집합이 비었는데 None이 아님")
            return
        rep.case(self.phase)
        if out is None:
            rep.viol(self.phase, "intersect", p, ctx, None, "교집합이 비어 있지 않은데 None 반환(값 손실)")
            return
        r = self.validate(out, "intersect", p, ctx, bits=bits)
        if r is None:
            return
        if not (r[0] <= lo and hi <= r[1]):
            rep.viol(self.phase, "intersect", p, ctx, {"out": list(r), "exact": [lo, hi]}, "정확한 교집합을 포함하지 않음")
        elif not (aiv[0] <= r[0] and r[1] <= aiv[1] and biv[0] <= r[0] and r[1] <= biv[1]):
            rep.viol("R0-API", "intersect", p, ctx, {"out": list(r), "exact": [lo, hi]}, "결과가 입력 구간 밖으로 넓어짐")

    # -- contains / width_q --
    def t_contains(self, p, A, aiv, qs):
        rep = self.rep
        rep.trial("contains", self.rand)
        for q in qs:
            ctx = {"x": list(aiv), "q": q}
            try:
                got = A.contains(q)
            except Exception as e:  # noqa
                rep.viol("R0-API", "contains", p, ctx, None, "예외 %s: %s" % (type(e).__name__, e))
                continue
            rep.case("R0-API")
            if not isinstance(got, bool) or got != (aiv[0] <= q <= aiv[1]):
                rep.viol("R0-API", "contains", p, ctx, got, "contains가 정확한 판정(닫힌 구간, bool)과 다름")

    def t_width(self, p, A, aiv):
        rep = self.rep
        rep.trial("width_q", self.rand)
        ctx = {"x": list(aiv)}
        try:
            w = A.width_q()
        except Exception as e:  # noqa
            rep.viol("R0-API", "width_q", p, ctx, None, "예외 %s: %s" % (type(e).__name__, e))
            return
        rep.case("R0-API")
        if not isinstance(w, Fraction) or w != aiv[1] - aiv[0]:
            rep.viol("R0-API", "width_q", p, ctx, w, "width_q가 hi_q-lo_q와 다름")


# ----------------------------------------------------------------------------
# R0-S2: 경계 사례 목록 (명시)
# ----------------------------------------------------------------------------
def s2_scalars(p):
    F = Fraction
    two_big = 1 << 2000
    vals = [
        # 0, 부호, 작은 정수
        F(0), F(1), F(-1), F(2), F(-2), F(3), F(5), F(10),
        # 정확히 표현되지 않는 유리수
        F(1, 3), F(-1, 3), F(2, 3), F(1, 10), F(-1, 10), F(1, 7), F(10 ** 30), F(1, 10 ** 30),
        # 매우 큰/작은 지수 (2^±2000)
        F(two_big), F(-two_big), F(1, two_big), F(-1, two_big),
        F(3 * (two_big >> 1)), F(1, 3 * two_big), F(two_big, 3),
        # 반올림 경계: 가수가 정확히 p, p+1, p+2 비트
        F(1 << p), F((1 << p) - 1), F((1 << p) + 1), F((1 << (p + 1)) - 1), F((1 << (p + 1)) + 1),
        F((1 << (p + 1)) + 2), F((1 << (p + 2)) + 3), F(-((1 << p) + 1)),
        F((1 << p) + 1, 1 << (p + 5)), F(-((1 << p) + 1), 1 << (p + 5)),
        F((1 << p) + 1) * two_big, F((1 << p) + 1, two_big), F(-((1 << p) + 1)) * two_big,
        # 1 근방
        1 + F(1, 1 << p), 1 + F(1, 1 << (p - 1)), 1 - F(1, 1 << p), 1 - F(1, 1 << (p + 1)),
        -1 - F(1, 1 << p),
    ]
    out, seen = [], set()
    for v in vals:
        if v not in seen:
            seen.add(v)
            out.append(v)
    return out


def s2_ranges(p):
    F = Fraction
    tb = 1 << 2000
    eps = F(1, 1 << (p + 40))
    return [
        (F(-1), F(1)), (F(0), F(1)), (F(-1), F(0)), (F(0), F(0)), (F(-2), F(3)), (F(-3), F(-2)), (F(2), F(3)),
        (F(-1, 3), F(1, 3)), (F(1, 3), F(2, 3)), (F(-2, 3), F(-1, 3)),
        (F(-1, tb), F(1, tb)), (F(-tb), F(tb)), (F(1, tb), F(tb)), (F(-tb), F(-1, tb)), (F(tb), F(2 * tb)),
        (F(1, 2 * tb), F(1, tb)),
        (F(1), F(1) + eps), (F(1) - eps, F(1)), (F(1, 3), F(1, 3) + eps), (-F(1, 3) - eps, -F(1, 3)),
        (F((1 << p) + 1), F((1 << p) + 3)), (F(-((1 << p) + 3)), F(-((1 << p) + 1))),
        (F(0), F(1, 3)), (F(-1, 3), F(0)), (F(0), F(tb)), (F(-tb), F(0)),
        (F(1), F(10 ** 30)), (F(-10 ** 30), F(10 ** 30)), (F(4), F(9)), (F(0), F(4)), (F(2), F(2)),
    ]


def s2_exception_cases(p):
    F = Fraction
    tb = 1 << 2000
    ctor_bad = [
        (F(1), F(0)), (F(0), F(-1)), (F(2), F(1)), (F(1, 3) + F(1, 10 ** 40), F(1, 3)),
        (F(1), F(1) - F(1, 1 << 3000)), (F(tb), F(-tb)), (F(1, tb), F(-1, tb)), (F(-1), F(-2)),
    ]
    div_zero_den = [(F(0), F(0)), (F(-1), F(1)), (F(0), F(1)), (F(-1), F(0)), (F(-1, tb), F(1, tb)),
                    (F(-tb), F(0)), (F(0), F(1, tb)), (F(-1, 3), F(1, 3))]
    sqrt_bad = [(F(-1), F(1)), (F(-1), F(-1, 2)), (F(-1, tb), F(1)), (F(-1, 3), F(-1, 3)),
                (F(-tb), F(tb)), (F(-1, 1 << (p + 40)), F(0))]
    sqrt_ok = [(F(0), F(0)), (F(0), F(1)), (F(1, tb), F(tb)), (F(2), F(2)), (F(0), F(1, tb))]
    return ctor_bad, div_zero_den, sqrt_bad, sqrt_ok


def run_exceptions(T, p):
    rep, DI = T.rep, T.DI
    ctor_bad, div_zero_den, sqrt_bad, sqrt_ok = s2_exception_cases(p)
    for lo, hi in ctor_bad:
        ctx = {"lo": lo, "hi": hi}
        for form in (0, 1):
            a, b = (as_scalar(T.rng, lo), as_scalar(T.rng, hi)) if form == 0 else (lo, hi)
            rep.case("R0-API")
            try:
                d = DI(a, b)
                rep.viol("R0-API", "ctor_lo_gt_hi", p, ctx, repr(d), "lo>hi인데 ValueError가 나오지 않음")
            except ValueError:
                pass
            except Exception as e:  # noqa
                rep.viol("R0-API", "ctor_lo_gt_hi", p, ctx, None, "ValueError가 아닌 예외 %s" % type(e).__name__)
    # 0을 포함하는 제수: DI/DI, DI/스칼라, 스칼라/DI
    num = T.mk(Fraction(1), Fraction(2), p)
    for lo, hi in div_zero_den:
        den = T.mk(lo, hi, p)
        if den is None or num is None:
            continue
        if not (den.lo_q <= 0 <= den.hi_q):
            continue
        for kind, f in (("DD", lambda: num / den), ("SD", lambda: Fraction(1) / den)):
            rep.case("R0-API")
            ctx = {"divisor": [den.lo_q, den.hi_q], "kind": kind}
            try:
                out = f()
                rep.viol("R0-API", "div_by_zero_interval_" + kind, p, ctx, repr(out), "ZeroDivisionError가 나오지 않음")
            except ZeroDivisionError:
                pass
            except TypeError as e:
                if kind == "SD":
                    rep.warning("reflected_op_unsupported", "div_SD", p, ctx, None, str(e))
                else:
                    rep.viol("R0-API", "div_by_zero_interval_" + kind, p, ctx, None, "TypeError %s" % e)
            except Exception as e:  # noqa
                rep.viol("R0-API", "div_by_zero_interval_" + kind, p, ctx, None, "ZeroDivisionError가 아닌 예외 %s" % type(e).__name__)
    if num is not None:
        for z in (0, Fraction(0)):
            rep.case("R0-API")
            ctx = {"divisor_scalar": z}
            try:
                out = num / z
                rep.viol("R0-API", "div_by_zero_scalar", p, ctx, repr(out), "스칼라 0으로 나눴는데 ZeroDivisionError가 아님")
            except ZeroDivisionError:
                pass
            except Exception as e:  # noqa
                rep.viol("R0-API", "div_by_zero_scalar", p, ctx, None, "ZeroDivisionError가 아닌 예외 %s" % type(e).__name__)
    for lo, hi in sqrt_bad:
        x = T.mk(lo, hi, p)
        if x is not None:
            T.t_sqrt(p, x, (x.lo_q, x.hi_q))
    for lo, hi in sqrt_ok:
        x = T.mk(lo, hi, p)
        if x is not None:
            T.t_sqrt(p, x, (x.lo_q, x.hi_q))


def run_s2(T, p, ivs_limit=None):
    """한 정밀도 p에서 경계 사례 전부."""
    T.set_phase("R0-S2")
    m = T.m
    m.set_prec(p)
    rep = T.rep
    rep.case("R0-API")
    if m.get_prec() != p:
        rep.viol("R0-API", "set_prec/get_prec", p, {"set": p}, m.get_prec(), "get_prec가 set_prec 값과 다름")
    items = []
    for v in s2_scalars(p):
        items.append((v, v))
    for lo, hi in s2_ranges(p):
        items.append((lo, hi))
    ivs = []
    for lo, hi in items:
        d = T.mk(lo, hi, p)
        if d is not None:
            ivs.append((d, (d.lo_q, d.hi_q)))
    scal = [Fraction(0), Fraction(1), Fraction(-1), Fraction(1, 3), Fraction(-2, 3), Fraction(3),
            Fraction(1 << 2000), Fraction(1, 1 << 2000), Fraction((1 << p) + 1), Fraction(10 ** 30)]
    # 이항(DI,DI)
    for A, aiv in ivs:
        for B, biv in ivs:
            for op in ("add", "sub", "mul", "div"):
                T.t_bin(op, p, A, B, aiv, biv, "DD")
            T.t_hull(p, A, B, aiv, biv, bits=True)
            T.t_intersect(p, A, B, aiv, biv, bits=True)
    # 스칼라 혼합
    for A, aiv in ivs:
        for s in scal:
            sc = as_scalar(T.rng, s)
            siv = (s, s)
            for op in ("add", "sub", "mul", "div"):
                T.t_bin(op, p, A, sc, aiv, siv, "DS")
                T.t_bin(op, p, sc, A, siv, aiv, "SD")
    # 단항, contains, width
    for A, aiv in ivs:
        T.t_unary("neg", p, A, aiv)
        T.t_unary("sqr", p, A, aiv)
        for n in range(6):
            T.t_unary("pow", p, A, aiv, n=n)
        T.t_sqrt(p, A, aiv)
        lo, hi = aiv
        tiny = Fraction(1, 1 << 4000)
        qs = [lo, hi, (lo + hi) / 2, lo - tiny, hi + tiny, Fraction(0), Fraction(1), Fraction(-1)]
        T.t_contains(p, A, aiv, qs)
        T.t_width(p, A, aiv)
    run_exceptions(T, p)
    return ivs


def run_s2_mixed(T):
    """서로 다른 정밀도: 생성은 p_in, 연산은 p_op. 포함만 엄밀 검사(복사형 연산의 비트 수는 건너뜀)."""
    T.set_phase("R0-S2")
    m = T.m
    for p_in, p_op in ((300, 24), (24, 300), (128, 53), (53, 128)):
        m.set_prec(p_in)
        items = [(v, v) for v in s2_scalars(p_in)[::2]] + s2_ranges(p_in)[::2]
        ivs = []
        for lo, hi in items:
            d = T.mk(lo, hi, p_in)
            if d is not None:
                ivs.append((d, (d.lo_q, d.hi_q)))
        m.set_prec(p_op)
        for A, aiv in ivs:
            for B, biv in ivs[::3]:
                for op in ("add", "sub", "mul", "div"):
                    T.t_bin(op, p_op, A, B, aiv, biv, "DD", mixed=True)
                T.t_hull(p_op, A, B, aiv, biv, bits=False)
                T.t_intersect(p_op, A, B, aiv, biv, bits=False)
            T.t_unary("neg", p_op, A, aiv, bits=False)
            T.t_unary("sqr", p_op, A, aiv)
            T.t_unary("pow", p_op, A, aiv, n=3)
            T.t_unary("pow", p_op, A, aiv, n=4)
            T.t_sqrt(p_op, A, aiv)


# ----------------------------------------------------------------------------
# R0-S1: 무작위 단계
# ----------------------------------------------------------------------------
TRIAL_KINDS = (
    [("bin", op) for op in ("add", "sub", "mul", "div") for _ in range(4)]
    + [("ds", op) for op in ("add", "sub", "mul", "div")]
    + [("sd", op) for op in ("add", "sub", "mul", "div")]
    + [("neg", None)] + [("sqr", None)] * 2
    + [("pow", n) for n in range(6)]
    + [("sqrt", None)] * 4 + [("hull", None)] * 2 + [("intersect", None)] * 3
    + [("contains", None), ("width", None)]
)


def run_s1(T, trials, deadline):
    T.set_phase("R0-S1")
    rng, m, rep = T.rng, T.m, T.rep
    done = 0
    while done < trials:
        if done % 512 == 0 and time.time() > deadline:
            return done, False
        p = rng.choice(PRECS)
        m.set_prec(p)
        kind, arg = rng.choice(TRIAL_KINDS)
        done += 1
        if kind == "bin" or kind in ("ds", "sd"):
            lo, hi = rand_iv(rng, p)
            A = T.mk(lo, hi, p)
            if A is None:
                continue
            aiv = (A.lo_q, A.hi_q)
            if kind == "bin":
                for _ in range(5):
                    lo2, hi2 = rand_iv(rng, p)
                    if arg == "div" and lo2 <= 0 <= hi2 and rng.random() < 0.8:
                        continue
                    break
                B = T.mk(lo2, hi2, p)
                if B is None:
                    continue
                T.t_bin(arg, p, A, B, aiv, (B.lo_q, B.hi_q), "DD")
            else:
                s = rand_val(rng, p)
                if arg == "div" and kind == "ds" and rng.random() < 0.9:
                    while s == 0:  # 제수(스칼라)는 대부분 0이 아니게
                        s = rand_val(rng, p)
                sc = as_scalar(rng, s)
                if kind == "ds":
                    T.t_bin(arg, p, A, sc, aiv, (s, s), "DS")
                else:
                    T.t_bin(arg, p, sc, A, (s, s), aiv, "SD")
        elif kind in ("neg", "sqr", "pow"):
            lo, hi = rand_iv(rng, p)
            A = T.mk(lo, hi, p)
            if A is None:
                continue
            T.t_unary(kind, p, A, (A.lo_q, A.hi_q), n=arg)
        elif kind == "sqrt":
            lo, hi = rand_iv(rng, p)
            if lo < 0 and rng.random() < 0.85:
                lo, hi = min(abs(lo), abs(hi)), max(abs(lo), abs(hi))
            A = T.mk(lo, hi, p)
            if A is None:
                continue
            T.t_sqrt(p, A, (A.lo_q, A.hi_q))
        elif kind in ("hull", "intersect"):
            lo, hi = rand_iv(rng, p)
            A = T.mk(lo, hi, p)
            if rng.random() < 0.6 and A is not None:   # 겹치는 구간을 자주 만든다
                c = (A.lo_q + A.hi_q) / 2
                w = max(A.hi_q - A.lo_q, Fraction(1, 1 << 8))
                lo2 = c - w * Fraction(rng.randint(0, 100), 100)
                hi2 = c + w * Fraction(rng.randint(0, 100), 100)
                if rng.random() < 0.3:
                    lo2, hi2 = A.hi_q, A.hi_q + w
            else:
                lo2, hi2 = rand_iv(rng, p)
            B = T.mk(lo2, hi2, p)
            if A is None or B is None:
                continue
            aiv, biv = (A.lo_q, A.hi_q), (B.lo_q, B.hi_q)
            if kind == "hull":
                T.t_hull(p, A, B, aiv, biv, bits=True)
            else:
                T.t_intersect(p, A, B, aiv, biv, bits=True)
        elif kind == "contains":
            lo, hi = rand_iv(rng, p)
            A = T.mk(lo, hi, p)
            if A is None:
                continue
            aiv = (A.lo_q, A.hi_q)
            qs = [aiv[0], aiv[1]] + pts(rng, *aiv)[2:] + [rand_val(rng, p)]
            T.t_contains(p, A, aiv, qs)
        else:  # width
            lo, hi = rand_iv(rng, p)
            A = T.mk(lo, hi, p)
            if A is not None:
                T.t_width(p, A, (A.lo_q, A.hi_q))
    return done, True


# ----------------------------------------------------------------------------
# 실행/보고
# ----------------------------------------------------------------------------
def load_module(path):
    path = os.path.abspath(path)
    sys.path.insert(0, os.path.dirname(path))
    spec = importlib.util.spec_from_file_location("dyad_under_test", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["dyad_under_test"] = mod
    spec.loader.exec_module(mod)
    return mod, path


def quantiles(xs):
    xs = sorted(xs)
    n = len(xs)
    if not n:
        return None
    pick = lambda f: xs[min(n - 1, int(f * n))]  # noqa
    return {"n": n, "min": xs[0], "median": pick(0.5), "p90": pick(0.9), "p99": pick(0.99), "max": xs[-1]}


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        h.update(f.read())
    return h.hexdigest()


def run_check(path, seed, trials, time_limit):
    t0 = time.time()
    result = {
        "tool": "check_dyad.py", "checker_version": CHECKER_VERSION,
        "module_path": os.path.abspath(path), "seed": seed, "trials_target": trials,
        "precisions": list(PRECS), "python": sys.version.split()[0],
        "claim": "시험 범위 안에서 위반 없음일 뿐 수학적 증명이 아님(계획 7절).",
    }
    try:
        result["module_sha256"] = sha256_file(path)
        mod, _ = load_module(path)
        for nm in ("DI", "set_prec", "get_prec"):
            if not hasattr(mod, nm):
                raise AttributeError("모듈에 %s 없음" % nm)
    except Exception as e:  # noqa
        result.update({"verdict": "ERROR", "error": "%s: %s" % (type(e).__name__, e)})
        return result
    rep = Report()
    rng = random.Random(seed)
    T = Tester(mod, rep, rng)
    try:
        # 기본 정밀도 128 (모듈을 막 불러온 상태)
        rep.case("R0-API")
        if mod.get_prec() != 128:
            rep.viol("R0-API", "default_prec", 128, {}, mod.get_prec(), "기본 정밀도가 128이 아님")
        for p in PRECS:                      # S2
            run_s2(T, p)
        run_s2_mixed(T)
        s2_trials = sum(rep.trials.values())
        deadline = t0 + time_limit
        done, completed = run_s1(T, trials, deadline)  # S1
        mod.set_prec(128)
    except Exception as e:  # noqa
        import traceback
        result.update({"verdict": "ERROR", "error": "검사 중 예외 %s: %s" % (type(e).__name__, e),
                       "traceback": traceback.format_exc()[-2000:]})
        return result
    result["elapsed_s"] = round(time.time() - t0, 1)
    result["s1_random_trials_completed"] = done
    result["s2_phase_op_evaluations"] = s2_trials
    result["point_membership_checks"] = rep.point_checks
    result["op_evaluations_by_name"] = dict(sorted(rep.trials.items()))
    checks = {}
    any_viol = False
    incomplete = False
    desc = {
        "R0-API": "예외 계약·타입·불변식·contains/width_q/hull/intersect 의미(보조, 계획 2절)",
        "R0-S1": "무작위 포함 시험",
        "R0-S2": "경계 사례 포함 + 정확히 표현되는 값 보존",
        "R0-S3": "sqrt 포함(제곱 비교)",
        "R0-S4": "끝점 가수 비트 수 <= prec",
    }
    for cid in IDS:
        nviol = sum(b["count"] for b in rep.by[cid].values())
        status = "PASS" if nviol == 0 else "FAIL"
        note = None
        if cid == "R0-S1" and nviol == 0 and (done < REQUIRED_S1_TRIALS or not completed):
            status, note = "INSUFFICIENT", "무작위 시험 %d건 < 요구 %d건" % (done, REQUIRED_S1_TRIALS)
            incomplete = True
        allex = [e for b in rep.by[cid].values() for e in b["smallest"]]
        mn = min(allex, key=lambda e: e["_size"]) if allex else None
        checks[cid] = {
            "description": desc[cid], "status": status, "cases_checked": rep.cases[cid],
            "violations": nviol, "note": note,
            "minimal_counterexample": strip_private(mn) if mn else None,
            "minimal_counterexample_raw": mn["_raw"] if mn and len(json.dumps(mn["_raw"])) < 4000 else None,
            "violations_by_op": {
                op: {"count": b["count"], "first": strip_private(b["first"]),
                     "smallest": [strip_private(e) for e in b["smallest"]]}
                for op, b in sorted(rep.by[cid].items())},
        }
        any_viol = any_viol or nviol > 0
    result["checks"] = checks
    result["warnings"] = rep.warn
    result["t1_tightness_ratio_outwidth_over_max_exactwidth_or_ulp"] = {
        k: quantiles(v) for k, v in sorted(rep.t1.items())}
    result["not_checked"] = [
        "R0-P1 속도(계획상 측정 보고 항목, 이 검사기 범위 밖)",
        "sqrt의 T1 조임 정도(정확한 폭이 무리수라 생략)",
        "음수 지수 pow, 비정수 지수, set_prec(<=0) 등 계약에 정의되지 않은 동작",
        "다른 정밀도로 만든 입력의 neg/hull/intersect/pow1 결과 가수 비트 수(계약 모호)",
        "스레드 안전성, 매우 긴 연쇄 연산(ODE 단계)에서의 누적 거동",
    ]
    result["verdict"] = "FAIL" if any_viol else ("INCOMPLETE" if incomplete else "PASS")
    return result


# ----------------------------------------------------------------------------
# 자기시험
# ----------------------------------------------------------------------------
def run_selftest(args):
    here = os.path.dirname(os.path.abspath(__file__))
    fakes = sorted(f for f in os.listdir(here) if f.startswith("fake_dyad_") and f.endswith(".py"))
    me = os.path.abspath(__file__)
    jobs = []
    for f in fakes:
        is_correct = f == "fake_dyad_correct.py"
        trials = args.correct_trials if is_correct else args.fake_trials
        out = os.path.join(args.workdir, f.replace(".py", ".check.json"))
        cmd = [sys.executable, "-I", me, os.path.join(here, f), "--seed", str(args.seed),
               "--trials", str(trials), "--out", out]
        jobs.append((f, is_correct, trials, out, subprocess.Popen(
            cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)))
    rows = []
    for f, is_correct, trials, out, proc in jobs:
        _, err = proc.communicate()
        try:
            with open(out, encoding="utf-8") as fh:
                res = json.load(fh)
        except Exception as e:  # noqa
            rows.append({"fake": f, "error": "결과 파일 없음: %s %s" % (e, err.decode("utf-8", "replace")[-300:])})
            continue
        spec = importlib.util.spec_from_file_location("_x_" + f[:-3], os.path.join(here, f))
        sys.path.insert(0, here)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        expected = list(getattr(mod, "EXPECTED_CAUGHT_BY", ()))
        hit = [cid for cid, c in res.get("checks", {}).items() if c["violations"] > 0]
        rows.append({
            "fake": f, "kind": "correct" if is_correct else "buggy",
            "bug": getattr(mod, "BUG_DESCRIPTION", ""),
            "trials": trials, "verdict": res.get("verdict"), "elapsed_s": res.get("elapsed_s"),
            "ids_with_violations": hit, "expected_ids": expected,
            "caught": (res.get("verdict") == "FAIL") if not is_correct else None,
            "expected_id_hit": bool(set(hit) & set(expected)) if expected else None,
            "false_positive": (res.get("verdict") != "PASS") if is_correct else None,
            "violation_counts": {cid: c["violations"] for cid, c in res.get("checks", {}).items()},
            "warnings": {k: v["count"] for k, v in res.get("warnings", {}).items()},
            "module_sha256": res.get("module_sha256"),
        })
    buggy = [r for r in rows if r.get("kind") == "buggy"]
    correct = [r for r in rows if r.get("kind") == "correct"]
    summary = {
        "buggy_total": len(buggy), "buggy_caught": sum(1 for r in buggy if r.get("caught")),
        "buggy_expected_id_hit": sum(1 for r in buggy if r.get("expected_id_hit")),
        "correct_total": len(correct), "correct_false_positives": sum(1 for r in correct if r.get("false_positive")),
        "errors": [r for r in rows if "error" in r],
    }
    return {"selftest": True, "seed": args.seed, "summary": summary, "rows": rows}


def main():
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    ap = argparse.ArgumentParser(description="Gate R0 독립 검사기 (Fraction 기준)")
    ap.add_argument("module", nargs="?", help="검사할 모듈 경로(DI, set_prec, get_prec 제공)")
    ap.add_argument("--seed", type=int, default=20261011)
    ap.add_argument("--trials", type=int, default=REQUIRED_S1_TRIALS)
    ap.add_argument("--time-limit", type=float, default=3300.0, help="무작위 단계 시간 상한(초)")
    ap.add_argument("--out", help="JSON 결과 파일(생략하면 표준출력만)")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--fake-trials", type=int, default=20000)
    ap.add_argument("--correct-trials", type=int, default=REQUIRED_S1_TRIALS)
    ap.add_argument("--workdir", default=None, help="자기시험 중간 JSON 폴더")
    args = ap.parse_args()
    if args.selftest:
        import tempfile
        args.workdir = args.workdir or tempfile.mkdtemp(prefix="check_dyad_selftest_")
        os.makedirs(args.workdir, exist_ok=True)
        res = run_selftest(args)
        text = json.dumps(res, ensure_ascii=False, indent=1)
        if args.out:
            with open(args.out, "w", encoding="utf-8") as f:
                f.write(text)
        print(json.dumps(res["summary"], ensure_ascii=False, indent=1))
        s = res["summary"]
        return 0 if (s["buggy_caught"] == s["buggy_total"] and s["correct_false_positives"] == 0 and not s["errors"]) else 1
    if not args.module:
        ap.error("module 경로가 필요합니다")
    res = run_check(args.module, args.seed, args.trials, args.time_limit)
    text = json.dumps(res, ensure_ascii=False, indent=1)
    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            f.write(text)
    print(text)
    return {"PASS": 0, "FAIL": 1, "INCOMPLETE": 3, "ERROR": 2}.get(res["verdict"], 2)


if __name__ == "__main__":
    sys.exit(main())
