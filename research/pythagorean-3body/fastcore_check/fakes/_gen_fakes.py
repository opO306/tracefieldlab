# -*- coding: utf-8 -*-
"""_gen_fakes.py - fakes/fake_*.py 를 표에서 생성한다 (재현용). 실행: python -I _gen_fakes.py

각 가짜는 원본 소스의 **한 군데**만 바꾼다(anchor 는 정확히 1번 나와야 함 — _mutant.py 가 확인).
EXPECT:
    detect      검사기가 반드시 잡아야 함
    equivalent  비트 단위로 같은 결과를 내는 '동치 변형' — 못 잡는 것이 정상(잡으면 오히려 놀랄 일)
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))

D = "dyad_patches"
T = "taylor_patches"
L = "lohner_patches"

FAKES = [
    # ---- dyad (구간 산술 기본 연산) -------------------------------------------------
    dict(tag="dn_trunc", expect="detect", kw={D: [(
        "        m >>= s\n        e += s\n",
        "        m = (m >> s) if m >= 0 else -((-m) >> s)\n        e += s\n")]},
        desc="_dn(하한 내림)이 음수에서 0 쪽으로 자름(floor 대신 truncate) — GMP tdiv 를 쓴 이식 실수"),
    dict(tag="up_floor", expect="detect", kw={D: [(
        "        m = -((-m) >> s)\n", "        m = m >> s\n")]},
        desc="_up(상한 올림)이 올림 대신 내림"),
    dict(tag="dn_nearest", expect="detect", kw={D: [(
        "        m >>= s\n        e += s\n",
        "        m = (m + (1 << (s - 1))) >> s\n        e += s\n")]},
        desc="_dn 이 내림 대신 가장 가까운 값으로 반올림(반올림 동률 포함)"),
    dict(tag="canon_off", expect="detect", kw={D: [(
        "    t = (m & -m).bit_length() - 1\n    if t:\n        return (m >> t, e + t)\n    return (m, e)\n",
        "    return (m, e)\n")]},
        desc="_canon 이 꼬리 0 비트를 안 지움 — 같은 값이지만 (가수, 지수) 표현이 다름"),
    dict(tag="mul_minmax", expect="detect", kw={D: [
        ("            if _lt(t, lo):\n", "            if t[0] < lo[0]:\n"),
        ("            if _lt(hi, t):\n", "            if hi[0] < t[0]:\n")]},
        desc="구간 곱의 최소·최대 선택이 지수 정렬 없이 가수만 비교"),
    dict(tag="sqrt_up_floor", expect="detect", kw={D: [(
        "    if not exact:\n        s += 1\n", "    if False:\n        s += 1\n")]},
        desc="sqrt 상한이 완전제곱이 아닐 때 +1 을 빼먹음(올림 누락)"),
    dict(tag="sqrt_exact_ignored", expect="detect", kw={D: [(
        "    return s, s * s == n, e // 2 - k\n", "    return s, False, e // 2 - k\n")]},
        desc="sqrt 의 '정확히 떨어짐' 판정을 빼먹음 — 완전제곱(예: 3-4-5 시작 배치)에서만 상한이 한 칸 넓어짐"),
    dict(tag="exp_wrap32", expect="detect", kw={},
        prelude="import dyad as _dy\nfrom lohner_step import _GD as _G\n\n\ndef _w(e):\n"
                "    return ((e + 0x80000000) % 0x100000000) - 0x80000000\n\n\ndef _ep(t):\n"
                "    return (t[0], _w(t[1])) if t[0] else t\n\n\ndef _wrap(x):\n"
                "    if x is None:\n        return None\n    if type(x) is list:\n"
                "        return [_wrap(i) for i in x]\n    if type(x) is _dy.DI:\n"
                "        return _dy.DI._new(_ep(x._l), _ep(x._h))\n"
                "    return _G(_wrap(x.v), [_wrap(g) for g in x.g])\n",
        extra="post=_wrap",
        desc="결과 끝점의 지수를 부호 있는 32비트로 감아 돌려줌(어댑터가 지수를 int32 로 저장하는 실수) — 지수가 ±2^31 을 넘는 입력에서만 틀림"),
    dict(tag="qdn_trunc", expect="detect", kw={D: [(
        "    q = (ma << k) // mb\n", "    q = (abs(ma << k) // mb) * (1 if ma >= 0 else -1)\n")]},
        desc="나눗셈 하한 몫을 내림 대신 0 쪽으로 자름(GMP tdiv 이식 실수) — 분자가 음수일 때만 틀림"),
    dict(tag="sqr_negative_branch", expect="detect", kw={D: [(
        "        return DI._new(_dn(*_powx(hi, n), p), _up(*_powx(lo, n), p))\n",
        "        return DI._new(_dn(*_powx(lo, n), p), _up(*_powx(hi, n), p))\n")]},
        desc="x**2 에서 구간이 전부 음수일 때 끝점을 맞바꿈(하한이 상한보다 커짐)"),
    dict(tag="qdn_ceil", expect="detect", kw={D: [(
        "    q = (ma << k) // mb\n", "    q = -((-(ma << k)) // mb)\n")]},
        desc="나눗셈 하한 몫을 내림 대신 올림으로 계산(반올림 방향 효과가 드문 변형)"),
    dict(tag="qup_nobump", expect="detect", kw={D: [(
        "    if r:\n        q += 1\n    return _up(q, ea - eb - k, p)\n",
        "    return _up(q, ea - eb - k, p)\n")]},
        desc="나눗셈 상한 몫에서 나머지가 있어도 +1 하지 않음"),
    dict(tag="sqr_as_mul", expect="detect", kw={D: [(
        "        return self ** 2\n", "        return self * self\n")]},
        desc="sqr() 를 x*x 로 대체 — 0 을 가로지르는 구간에서 하한이 음수로 나옴"),
    # ---- taylor_step (계수 재귀의 연산 순서) -------------------------------------------
    dict(tag="sqconv_cauchy", expect="detect", kw={T: [(
        "def _sq_conv(a, k):\n",
        "def _sq_conv(a, k):\n    return _conv(a, a, k)\n\n\ndef _sq_conv_unused(a, k):\n")]},
        desc="제곱 코시 곱에서 대칭 쌍을 한 번 계산해 2배 하지 않고 두 번 곱해 모두 더함"),
    dict(tag="sqconv_mid_mul", expect="detect", kw={T: [(
        "        m = a[k // 2].sqr()\n", "        m = a[k // 2] * a[k // 2]\n")]},
        desc="_sq_conv 의 가운데 항을 sqr 대신 곱셈으로"),
    dict(tag="conv_reversed", expect="detect", kw={T: [(
        "def _conv(a, b, k):\n",
        "def _conv(a, b, k):\n    acc = a[k] * b[0]\n    for j in range(k - 1, -1, -1):\n"
        "        acc = acc + a[j] * b[k - j]\n    return acc\n\n\ndef _conv_unused(a, b, k):\n")]},
        desc="코시 곱 합산 순서를 j=k..0 으로 뒤집음"),
    dict(tag="pow_div_twice", expect="detect", kw={T: [(
        "    return -(acc / (s0 * (2 * k)))\n", "    return -(acc / s0 / (2 * k))\n")]},
        desc="w_[k] 의 나눗셈을 한 번 대신 두 번으로 쪼갬(반올림 두 번)"),
    dict(tag="assoc_order", expect="detect", kw={T: [(
        "    acc = (s[1] * w[k - 1]) * (1 + 2 * k)\n", "    acc = s[1] * (w[k - 1] * (1 + 2 * k))\n")]},
        desc="곱셈 결합 순서 변경 (s*w)*c -> s*(w*c)"),
    dict(tag="mass_hiprec", expect="detect", kw={T: [(
        "        M = [DI(m) for m in masses]",
        "        _q = dyad.get_prec(); dyad.set_prec(_q + 64); M = [DI(m) for m in masses]; dyad.set_prec(_q)")]},
        desc="질량 상수 DI(m) 를 호출 정밀도보다 64비트 높게 만듦"),
    dict(tag="singular_lt", expect="detect", kw={T: [(
        "        if s0.lo_q <= 0:\n", "        if s0.lo_q < 0:\n")]},
        desc="특이점 판정 경계 오류: lo <= 0 대신 lo < 0 (거리 하한이 정확히 0 인 접촉 배치)"),
    dict(tag="pairs_reversed", expect="detect", kw={T: [(
        "            for p, (a, b) in enumerate(pairs):\n",
        "            for p, (a, b) in reversed(list(enumerate(pairs))):\n")]},
        desc="가속도 합산에서 쌍(pair) 순서를 뒤집음 (3체부터 반올림이 달라짐)"),
    dict(tag="div_int_recip", expect="detect", kw={T: [(
        "    return a / n\n", "    return a * DI(Fraction(1, n))\n")]},
        desc="a/(k+1) 를 a * DI(1/(k+1)) 로(역수를 먼저 반올림)"),
    dict(tag="k_short", expect="detect", kw={T: [(
        "        for k in range(K):\n            acc = [DI(0)] * (2 * n)",
        "        for k in range(max(K - 1, 0)):\n            acc = [DI(0)] * (2 * n)")]},
        desc="재귀를 한 차수 덜 돌려 계수 목록이 짧음"),
    dict(tag="consts_at_make", expect="detect", kw={T: [(
        "    def coeffs(X, K):\n        M = [DI(m) for m in masses]",
        "    _Mpre = [DI(m) for m in masses]\n\n    def coeffs(X, K):\n        M = _Mpre")]},
        desc="질량 상수를 system 을 만들 때의 정밀도로 고정(결정 A 위반: 호출 시점 정밀도여야 함)"),
    dict(tag="consts_cached", expect="detect", kw={T: [(
        "    def coeffs(X, K):\n        M = [DI(m) for m in masses]",
        "    _cache = []\n\n    def coeffs(X, K):\n        if not _cache:\n"
        "            _cache.append([DI(m) for m in masses])\n        M = _cache[0]")]},
        desc="질량 상수를 첫 호출 때의 정밀도로 캐시하고 이후 호출에 재사용"),
    # ---- lohner_step (_GD: 값 + 기울기 쌍) ----------------------------------------------
    dict(tag="gd_div_split", expect="detect", kw={L: [(
        "                    g.append((a - q * b) / bv)\n", "                    g.append(a / bv - (q * b) / bv)\n")]},
        desc="_GD 나눗셈 미분 (da - q*db)/b 를 da/b - (q*db)/b 로 쪼갬"),
    dict(tag="gd_rsub_sign", expect="detect", kw={L: [(
        "        return _GD(o - self.v, [None if a is None else -a for a in self.g])\n",
        "        return _GD(o - self.v, [None if a is None else a for a in self.g])\n")]},
        desc="_GD 의 (상수 - _GD) 에서 기울기 부호를 뒤집지 않음"),
    dict(tag="gd_rtruediv_sign", expect="detect", kw={L: [(
        "        return _GD(q, [None if a is None else -(qu * a) for a in self.g])\n",
        "        return _GD(q, [None if a is None else (qu * a) for a in self.g])\n")]},
        desc="_GD 의 (상수 / _GD) 미분에서 음수 부호를 빠뜨림"),
    dict(tag="gd_sqr_factor", expect="detect", kw={L: [(
        "        return _GD(av.sqr(), [None if a is None else (av * a) * 2 for a in self.g])\n",
        "        return _GD(av.sqr(), [None if a is None else (av * a) for a in self.g])\n")]},
        desc="_GD.sqr() 의 미분에서 2배를 빠뜨림"),
    dict(tag="gd_sqrt_factor", expect="detect", kw={L: [(
        "        r2 = r * 2\n", "        r2 = r\n")]},
        desc="_GD.sqrt() 의 미분 ds/(2 sqrt s) 에서 2를 빠뜨림"),
    dict(tag="gd_mul_none", expect="detect", kw={L: [(
        "                    g.append(None if b is None else av * b)\n", "                    g.append(None)\n")]},
        desc="_GD 곱셈에서 한쪽 기울기가 None 일 때 다른 쪽 항을 버림"),
    # ---- 정밀도·상태·예외·모양 (소스 패치 없이 감싸기) --------------------------------------
    dict(tag="prec_plus_one", expect="detect", kw={}, extra="prec_delta=1",
        desc="내부 계산을 요구 정밀도보다 1비트 높게"),
    dict(tag="prec_at_make", expect="detect", kw={}, extra="prec_at_make=True",
        desc="모든 계산을 system 을 만든 때의 정밀도로(호출 시점 정밀도 무시)"),
    dict(tag="prec_clamp128", expect="detect", kw={}, extra="prec_fn=lambda p: min(p, 128)",
        desc="정밀도가 128비트를 넘으면 128로 잘라서 계산(고정 폭 한계 이식 실수)"),
    dict(tag="exc_wide", expect="detect", kw={},
        prelude="from dyad import DI as _DI\n\n\ndef _wide(system, K):\n"
                "    big = _DI._new((-1, 8192), (1, 8192))\n"
                "    return [[big] * (K + 1) for _ in range(system.dim)]\n",
        extra="on_singular=_wide",
        desc="특이점에서 예외 대신 거대한 구간을 돌려줌"),
    dict(tag="exc_type", expect="detect", kw={},
        prelude="def _zde(system, K):\n    raise ZeroDivisionError('singular')\n",
        extra="on_singular=_zde",
        desc="특이점에서 _Singular 대신 ZeroDivisionError"),
    dict(tag="none_to_di0", expect="detect", kw={},
        prelude="from dyad import DI as _DI\nfrom lohner_step import _GD as _G\n\n\ndef _dense(X):\n"
                "    return [(_G(x.v, [_DI(0) if g is None else g for g in x.g]) if type(x) is _G else x)\n"
                "            for x in X]\n",
        extra="pre=_dense",
        desc="_GD 입력의 None(정확히 0)을 DI(0) 으로 바꿔 계산 — 값은 같지만 None 위치가 달라짐"),
    dict(tag="zero_to_none", expect="detect", kw={},
        prelude="from lohner_step import _GD as _G\n\n\ndef _sparsify(res):\n"
                "    out = []\n    for comp in res:\n        row = []\n        for e in comp:\n"
                "            if type(e) is _G:\n"
                "                e = _G(e.v, [None if (g is not None and g._l == (0, 0) and g._h == (0, 0)) else g\n"
                "                             for g in e.g])\n"
                "            row.append(e)\n        out.append(row)\n    return out\n",
        extra="post=_sparsify",
        desc="결과 기울기의 정확한 0 인 DI(0) 을 None 으로 바꿈"),
    dict(tag="mutates_input", expect="detect", kw={},
        prelude="def _mutate(X):\n    X.append(X[0])\n    return X\n",
        extra="pre=_mutate",
        desc="입력 리스트 X 를 제자리에서 바꿈(끝에 원소 추가)"),
    dict(tag="prec_side_effect", expect="detect", kw={},
        prelude="import dyad as _dy\n\n\ndef _bump(res):\n    _dy.set_prec(_dy.get_prec() + 1)\n    return res\n",
        extra="post=_bump",
        desc="계산 뒤에 전역 정밀도를 바꿔 놓음"),
    dict(tag="container_tuple", expect="detect", kw={},
        prelude="def _tup(res):\n    return tuple(tuple(c) for c in res)\n",
        extra="post=_tup",
        desc="결과를 리스트가 아닌 튜플로 돌려줌(모양 차이만; 값은 같음)"),
    # ---- 동치 변형: 비트 단위로 같아서 못 잡는 것이 정상 ------------------------------------------
    dict(tag="eq_acc_plus_acc", expect="equivalent", kw={T: [(
        "        acc = acc * 2", "        acc = acc + acc")]},
        desc="[동치] 2배를 acc*2 대신 acc+acc 로 (2의 거듭제곱 곱은 반올림과 교환되므로 비트 동일)"),
    dict(tag="eq_acc_times_di2", expect="equivalent", kw={T: [(
        "        acc = acc * 2", "        acc = acc * DI(2)")]},
        desc="[동치] 2배를 정수 2 대신 DI(2) 로"),
    dict(tag="eq_add_commuted", expect="equivalent", kw={T: [(
        "                acc[2 * a] = acc[2 * a] + M[b] * ux\n",
        "                acc[2 * a] = M[b] * ux + acc[2 * a]\n")]},
        desc="[동치] 덧셈의 피연산자 순서 교환(구간 덧셈은 교환법칙이 비트 단위로 성립)"),
]

HEADER = ("# -*- coding: utf-8 -*-\n# 자동 생성(fakes/_gen_fakes.py) — 검사기 자기시험용 가짜 대상. 직접 고치지 말고 생성기를 고칠 것.\n"
          "from _mutant import build\n\n")


def render(f):
    kw = f["kw"]
    parts = []
    for k, v in kw.items():
        parts.append("%s=%s" % (k, repr(v)))
    if f.get("extra"):
        parts.append(f["extra"])
    body = HEADER
    body += "EXPECT = %r\nDESCRIPTION = %r\n\n" % (f["expect"], f["desc"])
    if f.get("prelude"):
        body += f["prelude"] + "\n\n"
    body += "_mutant, make_nbody_fast = build(%r%s)\n" % (f["tag"], "".join(", " + p for p in parts))
    return body


PROCESS_EXIT = '''# -*- coding: utf-8 -*-
# 수동 작성 — 충돌(프로세스 사망) 처리 시험용 가짜. 특정 입력에서 프로세스를 즉시 끝낸다(세그폴트 모사).
import os

import taylor_step

EXPECT = "detect"
DESCRIPTION = "특정 입력에서 프로세스가 죽음(세그폴트 모사: os._exit)"


class _Sys(object):
    def __init__(self, masses):
        self._inner = taylor_step.make_nbody(masses)
        self.name = self._inner.name
        self.dim = self._inner.dim

    def coeffs(self, X, K):
        x0 = X[0]
        v = x0.v if hasattr(x0, "v") else x0
        if K == 20 and v._h[0] % 16 == 5:
            os._exit(139)
        return self._inner.coeffs(X, K)


def make_nbody_fast(masses):
    return _Sys(masses)
'''


def main():
    names = []
    for f in FAKES:
        path = os.path.join(HERE, "fake_%s.py" % f["tag"])
        with open(path, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(render(f))
        names.append(f["tag"])
    with open(os.path.join(HERE, "fake_process_exit.py"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write(PROCESS_EXIT)
    print("wrote %d fakes + process_exit" % len(names))


if __name__ == "__main__":
    main()
