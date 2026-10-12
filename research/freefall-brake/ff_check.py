"""독립 검사기: Krawczyk 상자를 정확한 유리수 구간으로 다시 계산 (표준 라이브러리만, dyad 미사용).

입력 레코드의 g(m), [Dg] 구간은 그대로 믿는다(적분기 엄밀성은 별도 검사됨).
잡는 것: Krawczyk 식·행렬곱·포함 비교의 구현 오류, 포함을 깨뜨리는 파일 손상·수정,
X가 기록된 m ± r 와 맞지 않는 수정(상자를 넓히는 등).
못 잡는 것: g(m)·[Dg]를 좁히는 등 포함을 유지하는 수정 — 이를 막으려면 적분을 다시 해야 한다.
유일성: ‖I − C·[Dg]‖∞ < 1 이면 [Dg]의 모든 행렬이 가역 → 볼록한 X에서 g가 단사(평균값 정리)
→ 영점이 많아야 하나. 이 노름을 정확한 유리수로 계산해 "unique"로 보고한다.
"""
import json
import sys
from fractions import Fraction as F


class Q:
    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi=None):
        self.lo = F(lo)
        self.hi = F(lo if hi is None else hi)
        if self.lo > self.hi:
            raise ValueError("lo > hi")

    def __add__(self, o):
        o = _q(o)
        return Q(self.lo + o.lo, self.hi + o.hi)

    def __sub__(self, o):
        o = _q(o)
        return Q(self.lo - o.hi, self.hi - o.lo)

    def __mul__(self, o):
        o = _q(o)
        p = (self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi)
        return Q(min(p), max(p))


def _q(x):
    return x if isinstance(x, Q) else Q(x)


def _iv(pair):
    return Q(F(pair[0]), F(pair[1]))


def check(rec):
    X = [_iv(x) for x in rec["X"]]
    m = [F(v) for v in rec["m"]]
    C = [[F(v) for v in row] for row in rec["C"]]
    gm = [_iv(g) for g in rec["gm"]]
    Dg = [[_iv(e) for e in row] for row in rec["Dg"]]
    q = len(X)
    if not all(x.lo <= mi <= x.hi for x, mi in zip(X, m)):
        return {"verdict": "FAIL", "K": [], "note": "m not in X"}
    if "r" not in rec:
        return {"verdict": "FAIL", "K": [], "note": "record has no radii r"}
    r = [F(v) for v in rec["r"]]
    # X = 바깥 반올림한 [m − r, m + r]: 그 구간을 포함하고, 2r 이상 벗어나지 않아야 한다.
    if not all(x.lo <= mi - ri and mi + ri <= x.hi and mi - 2 * ri <= x.lo and x.hi <= mi + 2 * ri
               for x, mi, ri in zip(X, m, r)):
        return {"verdict": "FAIL", "K": [], "note": "X does not match m +- r"}
    K = []
    for i in range(q):
        acc = Q(m[i])
        for l in range(q):
            acc = acc - Q(C[i][l]) * gm[l]
        for j in range(q):
            a = Q(1 if i == j else 0)
            for l in range(q):
                a = a - Q(C[i][l]) * Dg[l][j]
            acc = acc + a * (X[j] - Q(m[j]))
        K.append(acc)
    norm = F(0)
    for i in range(q):
        row = F(0)
        for j in range(q):
            a = Q(1 if i == j else 0)
            for l in range(q):
                a = a - Q(C[i][l]) * Dg[l][j]
            row += max(abs(a.lo), abs(a.hi))
        norm = max(norm, row)
    inside = all(k.lo > x.lo and k.hi < x.hi for k, x in zip(K, X))
    return {"verdict": "PASS" if inside else "FAIL",
            "K": [["%d/%d" % (k.lo.numerator, k.lo.denominator),
                   "%d/%d" % (k.hi.numerator, k.hi.denominator)] for k in K],
            "unique": inside and norm < 1,
            "norm_I_minus_CDg": "%d/%d" % (norm.numerator, norm.denominator),
            "note": "algebra + inclusion + uniqueness norm; g(m), [Dg] trusted from sealed integrator"}


if __name__ == "__main__":
    bad = 0
    for path in sys.argv[1:]:
        with open(path, encoding="utf-8") as f:
            r = check(json.load(f))
        print(path, r["verdict"])
        bad += r["verdict"] != "PASS"
    sys.exit(1 if bad else 0)
