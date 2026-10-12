"""자체 시험이 '일부러 틀린 dyad.py'를 잡아내는지 확인한다 (변이 시험).

dyad.py를 임시 폴더로 복사해 한 곳만 틀리게 바꾸고, test_dyad_impl.py를 그 복사본에
대해 실행한다. 시험이 실패(= 변이를 잡음)해야 정상이다. 원본 dyad.py는 건드리지 않는다.

실행:  python mutation_check_impl.py
"""

import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))

# (이름, 원본 문자열, 바꿀 문자열)  - 원본 문자열은 dyad.py에서 정확히 1번 나와야 한다.
MUTATIONS = [
    ("dn_truncate_toward_zero",
     "    if s > 0:\n        m >>= s\n        e += s\n    return _canon(m, e)",
     "    if s > 0:\n        m = m >> s if m >= 0 else -((-m) >> s)\n        e += s\n    return _canon(m, e)"),
    ("up_floor_instead_of_ceil",
     "        m = -((-m) >> s)\n        e += s",
     "        m = m >> s\n        e += s"),
    ("qup_ignore_remainder",
     "    q, r = divmod(ma << k, mb)\n    if r:\n        q += 1\n    return _up",
     "    q, r = divmod(ma << k, mb)\n    return _up"),
    ("sqrt_up_no_plus_one",
     "    if not exact:\n        s += 1\n    return _up(s, e2, p)",
     "    return _up(s, e2, p)"),
    ("sqrt_dn_uses_up",
     "    s, _exact, e2 = _sqrt_parts(a, p)\n    return _dn(s, e2, p)",
     "    s, _exact, e2 = _sqrt_parts(a, p)\n    return _up(s + 1, e2, p)"),
    ("pow_even_straddle_lo_not_zero",
     "        return DI._new((0, 0), _up(*top, p))",
     "        return DI._new(_dn(*a, p), _up(*top, p))"),
    ("div_min_picks_wrong_divisor_end",
     "        lo = _qdn(xl, yh if xl[0] >= 0 else yl, p)\n        hi = _qup(xh, yl if xh[0] >= 0 else yh, p)\n        return DI._new(lo, hi)",
     "        lo = _qdn(xl, yl if xl[0] >= 0 else yh, p)\n        hi = _qup(xh, yl if xh[0] >= 0 else yh, p)\n        return DI._new(lo, hi)"),
    ("div_negative_divisor_no_flip",
     "            return DI._new(_dn(*_neg(hi), p), _up(*_neg(lo), p))",
     "            return DI._new(_dn(*lo, p), _up(*hi, p))"),
    ("mul_only_two_corners",
     "        t1 = _mulx(a, d)\n        t2 = _mulx(b, c)",
     "        t1 = t0\n        t2 = t3"),
    ("sub_wrong_pairing",
     "        return DI._new(_dn(*_addx(self._l, _neg(o._h)), p),\n                       _up(*_addx(self._h, _neg(o._l)), p))",
     "        return DI._new(_dn(*_addx(self._l, _neg(o._l)), p),\n                       _up(*_addx(self._h, _neg(o._h)), p))"),
    ("frac_ceil_ignores_remainder",
     "        return _dn(q, k, p), _up(q + 1 if r else q, k, p)",
     "        return _dn(q, k, p), _up(q, k, p)"),
    ("frac_neg_lower_uses_floor_mag",
     "    return _dn(-(q + 1 if r else q), k, p), _up(-q, k, p)",
     "    return _dn(-q, k, p), _up(-q, k, p)"),
    ("canon_drops_exponent_shift",
     "        return (m >> t, e + t)",
     "        return (m >> t, e)"),
    ("intersect_never_empty",
     "        if _lt(hi, lo):\n            return None\n        return DI._new(_dn(*lo, p), _up(*hi, p))",
     "        if _lt(hi, lo):\n            lo, hi = hi, lo\n        return DI._new(_dn(*lo, p), _up(*hi, p))"),
    ("div_zero_check_missing",
     "        if yl[0] <= 0 <= yh[0]:\n            raise ZeroDivisionError(\"divisor interval contains 0\")",
     "        if yl[0] == 0 and yh[0] == 0:\n            raise ZeroDivisionError(\"divisor interval contains 0\")"),
]


def main():
    src = open(os.path.join(HERE, "dyad.py"), encoding="utf-8").read()
    survivors = []
    for name, old, new in MUTATIONS:
        if src.count(old) != 1:
            print("BAD MUTATION SPEC (match count %d): %s" % (src.count(old), name))
            survivors.append(name + " (spec error)")
            continue
        with tempfile.TemporaryDirectory() as td:
            with open(os.path.join(td, "dyad.py"), "w", encoding="utf-8") as f:
                f.write(src.replace(old, new))
            shutil.copy(os.path.join(HERE, "test_dyad_impl.py"), td)
            r = subprocess.run([sys.executable, "-I", os.path.join(td, "test_dyad_impl.py")],
                               capture_output=True, text=True, cwd=td)
            caught = r.returncode != 0
            fails = [ln for ln in r.stderr.splitlines() if ln.startswith(("FAIL:", "ERROR:"))]
            print("%-34s %s  (%d failing tests)" % (name, "CAUGHT" if caught else "SURVIVED", len(fails)))
            if not caught:
                survivors.append(name)
    print("\nmutations: %d, caught: %d, survived: %d" % (
        len(MUTATIONS), len(MUTATIONS) - len(survivors), len(survivors)))
    for s in survivors:
        print("  SURVIVED:", s)
    return 1 if survivors else 0


if __name__ == "__main__":
    sys.exit(main())
