// fastcore.cpp - Gate FC1: 계수 재귀 핵심부의 C++17 이식 (원본과 비트 단위 동일을 목표).
//
// 원본(오라클, 수정 금지, import만):
//   rigorous/dyad.py        DI 끝점 연산 _canon _dn _up _lt _addx _mulx _powx _qdn _qup _sqrt_parts
//                           DI.__add__ __sub__ __neg__ __mul__ __pow__(2)=sqr _div sqrt, int 강제변환 _coerce
//   rigorous/lohner_step.py _GD 의 + - 단항- * / 역나눗셈 sqr sqrt (None = 정확한 0 표시)
//   rigorous/taylor_step.py make_nbody.coeffs, _conv, _sq_conv, _pow_m32_next, _div_int
// 대응표는 README_KO.md. 원칙: 원본의 연산 순서·분기·반올림 위치를 한 줄씩 그대로 옮긴다.
//
// 표현: 끝점 (m, e) = m * 2^e. m 은 GMP mpz, e 는 int64 (넘침 검사; 넘치면 Unsupported).
// 입력 끝점은 받은 그대로 둔다(표준형으로 바꾸지 않음): 원본도 입력 튜플을 그대로 쓰기 때문이다
// (_div_int(a, 1) 은 a 를 그대로 돌려줌, _qdn 의 k 는 가수 비트 길이에 의존).
//
// Unsupported: 아주 큰 지수 차이(왼쪽 밀기 > 2^26 비트), 지수 넘침, 원본 클래스가 아닌 입력 등.
// 이때 C++ 은 계산을 버리고 _fastcore.Unsupported 를 던진다. 어댑터(fast_system.py)는 그 호출을
// 원본 파이썬 경로로 다시 계산한다(결과는 원본 그 자체). 횟수는 통계에 남긴다.
//
// GIL 을 잡은 채로 실행한다(정적 임시 변수 사용). 스레드 안전성은 GIL 에 기댄다.

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <gmp.h>

#include <cstdint>
#include <cstring>
#include <memory>
#include <new>
#include <string>
#include <unordered_map>
#include <vector>

#if PY_VERSION_HEX < 0x030D0000
#error "fastcore needs Python >= 3.13 (PyLong_AsNativeBytes / PyLong_FromNativeBytes)"
#endif

namespace {

// ---------------------------------------------------------------------------
// 오류 종류
// ---------------------------------------------------------------------------
struct Unsupported { std::string why; };
struct SingularEx {};   // taylor_step._Singular 에 대응 (어댑터가 원본 클래스로 다시 던짐)
struct ZeroDivEx {};    // dyad: ZeroDivisionError("divisor interval contains 0")
struct SqrtNegEx {};    // dyad: ValueError("sqrt of interval with negative lower bound")
struct PyErrSet {};     // 파이썬 예외가 이미 설정됨

const int64_t EXP_LIMIT = int64_t(1) << 40;     // 입력 지수 |e| 상한 (넘으면 Unsupported)
const int64_t SHIFT_LIMIT = int64_t(1) << 26;   // 왼쪽 밀기 비트 수 상한 (넘으면 Unsupported)
const int64_t PREC_LIMIT = int64_t(1) << 24;    // 정밀도 상한
const long long K_LIMIT = 1000000;              // 차수 상한

inline int64_t eadd(int64_t a, int64_t b) {
  int64_t r;
  if (__builtin_add_overflow(a, b, &r)) throw Unsupported{"exponent overflow"};
  return r;
}
inline int64_t esub(int64_t a, int64_t b) {
  int64_t r;
  if (__builtin_sub_overflow(a, b, &r)) throw Unsupported{"exponent overflow"};
  return r;
}
inline int64_t emul(int64_t a, int64_t b) {
  int64_t r;
  if (__builtin_mul_overflow(a, b, &r)) throw Unsupported{"exponent overflow"};
  return r;
}
inline void shift_guard(int64_t d) {
  if (d < 0 || d > SHIFT_LIMIT) throw Unsupported{"shift too large"};
}

// RAII mpz
struct Z {
  mpz_t z;
  Z() { mpz_init(z); }
  Z(const Z& o) { mpz_init_set(z, o.z); }
  Z(Z&& o) noexcept { mpz_init(z); mpz_swap(z, o.z); }
  Z& operator=(const Z& o) {
    if (this != &o) mpz_set(z, o.z);
    return *this;
  }
  Z& operator=(Z&& o) noexcept {
    mpz_swap(z, o.z);
    return *this;
  }
  ~Z() { mpz_clear(z); }
};

struct EP {  // 끝점 (m, e)
  Z m;
  int64_t e = 0;
};
struct DIv {  // 구간 [l, h]
  EP l, h;
};

// int.bit_length() (부호 무시, 0 -> 0)
inline int64_t bitlen(const mpz_t m) {
  return mpz_sgn(m) == 0 ? 0 : (int64_t)mpz_sizeinbase(m, 2);
}

// dyad._canon: 꼬리 0 비트를 지수로 (값 불변)
inline void canon(EP& x) {
  if (mpz_sgn(x.m.z) == 0) {
    x.e = 0;
    return;
  }
  mp_bitcnt_t t = mpz_scan1(x.m.z, 0);  // 음수도 2의 보수 기준 최하위 1 비트 = |m| 과 같은 위치
  if (t) {
    mpz_tdiv_q_2exp(x.m.z, x.m.z, t);  // 나누어떨어짐 -> 정확
    x.e = eadd(x.e, (int64_t)t);
  }
}

// dyad._dn: 내림 (-무한대 쪽). out.m 이 m 과 같은 객체여도 됨.
inline void dn(EP& out, const mpz_t m, int64_t e, int64_t p) {
  if (mpz_sgn(m) == 0) {  // 원본: s<0 이라 밀기 없음 -> _canon(0,e) = (0,0)
    mpz_set_ui(out.m.z, 0);
    out.e = 0;
    return;
  }
  int64_t s = bitlen(m) - p;
  if (s > 0) {
    mpz_fdiv_q_2exp(out.m.z, m, (mp_bitcnt_t)s);  // 파이썬 >> (floor)
    out.e = eadd(e, s);
  } else {
    if (out.m.z != m) mpz_set(out.m.z, m);
    out.e = e;
  }
  canon(out);
}

// dyad._up: 올림 (+무한대 쪽). -((-m) >> s) = ceil(m / 2^s)
inline void up(EP& out, const mpz_t m, int64_t e, int64_t p) {
  if (mpz_sgn(m) == 0) {
    mpz_set_ui(out.m.z, 0);
    out.e = 0;
    return;
  }
  int64_t s = bitlen(m) - p;
  if (s > 0) {
    mpz_cdiv_q_2exp(out.m.z, m, (mp_bitcnt_t)s);
    out.e = eadd(e, s);
  } else {
    if (out.m.z != m) mpz_set(out.m.z, m);
    out.e = e;
  }
  canon(out);
}

// 정확한 값 비교 sign(a - b) (dyad._lt 는 vcmp(a,b) < 0). 값만으로 결정되므로 원본의 지수 맞춤 비교와 같다.
int vcmp(const mpz_t am, int64_t ae, const mpz_t bm, int64_t be) {
  static Z t;
  int sa = mpz_sgn(am), sb = mpz_sgn(bm);
  if (sa != sb) return sa < sb ? -1 : 1;
  if (sa == 0) return 0;
  int64_t ta = eadd(bitlen(am), ae);  // 최상위 비트 위치 + 1
  int64_t tb = eadd(bitlen(bm), be);
  int mag;
  if (ta != tb) {
    mag = ta < tb ? -1 : 1;
  } else if (ae >= be) {
    int64_t d = ae - be;  // = bitlen(bm) - bitlen(am), 작다
    mpz_mul_2exp(t.z, am, (mp_bitcnt_t)d);
    mag = mpz_cmpabs(t.z, bm);
  } else {
    int64_t d = be - ae;
    mpz_mul_2exp(t.z, bm, (mp_bitcnt_t)d);
    mag = -mpz_cmpabs(t.z, am);
  }
  if (mag > 0) mag = 1;
  if (mag < 0) mag = -1;
  return sa > 0 ? mag : -mag;
}
inline int vcmp(const EP& a, const EP& b) { return vcmp(a.m.z, a.e, b.m.z, b.e); }

// dyad._addx(a, b) 또는 _addx(a, _neg(b)) (negb). 결과 (rm, re) 는 반올림 전 정확한 합.
// rm 은 a.m, b.m 과 다른 객체여야 한다.
void addx(mpz_t rm, int64_t& re, const EP& a, const EP& b, bool negb) {
  int sa = mpz_sgn(a.m.z), sbv = mpz_sgn(b.m.z);
  if (sbv == 0) {  // if mb == 0: return a
    mpz_set(rm, a.m.z);
    re = a.e;
    return;
  }
  if (sa == 0) {  // if ma == 0: return b
    if (negb)
      mpz_neg(rm, b.m.z);
    else
      mpz_set(rm, b.m.z);
    re = b.e;
    return;
  }
  if (a.e >= b.e) {  // ((ma << (ea - eb)) + mb, eb)
    int64_t d = esub(a.e, b.e);
    shift_guard(d);
    mpz_mul_2exp(rm, a.m.z, (mp_bitcnt_t)d);
    if (negb)
      mpz_sub(rm, rm, b.m.z);
    else
      mpz_add(rm, rm, b.m.z);
    re = b.e;
  } else {  // (ma + (mb << (eb - ea)), ea)
    int64_t d = esub(b.e, a.e);
    shift_guard(d);
    mpz_mul_2exp(rm, b.m.z, (mp_bitcnt_t)d);
    if (negb)
      mpz_sub(rm, a.m.z, rm);
    else
      mpz_add(rm, a.m.z, rm);
    re = a.e;
  }
}

// dyad._mulx
inline void mulx(mpz_t rm, int64_t& re, const EP& a, const EP& b) {
  mpz_mul(rm, a.m.z, b.m.z);
  if (mpz_sgn(rm) == 0)
    re = 0;
  else
    re = eadd(a.e, b.e);
}

// dyad._powx(a, 2)
inline void powx2(mpz_t rm, int64_t& re, const EP& a) {
  if (mpz_sgn(a.m.z) == 0) {
    mpz_set_ui(rm, 0);
    re = 0;
    return;
  }
  mpz_mul(rm, a.m.z, a.m.z);
  re = emul(a.e, 2);
}

// ---------------------------------------------------------------------------
// DI 연산 (out 은 입력과 같은 객체여도 된다: 입력을 모두 읽은 뒤 쓴다)
// ---------------------------------------------------------------------------

// DI.__add__: [dn(l1+l2), up(h1+h2)]
void di_add(DIv& out, const DIv& x, const DIv& y, int64_t p) {
  static Z sl, sh;
  int64_t el, eh;
  addx(sl.z, el, x.l, y.l, false);
  addx(sh.z, eh, x.h, y.h, false);
  dn(out.l, sl.z, el, p);
  up(out.h, sh.z, eh, p);
}

// DI.__sub__: [dn(l1 + (-h2)), up(h1 + (-l2))]
void di_sub(DIv& out, const DIv& x, const DIv& y, int64_t p) {
  static Z sl, sh;
  int64_t el, eh;
  addx(sl.z, el, x.l, y.h, true);
  addx(sh.z, eh, x.h, y.l, true);
  dn(out.l, sl.z, el, p);
  up(out.h, sh.z, eh, p);
}

// DI.__neg__: [dn(-h), up(-l)]
void di_neg(DIv& out, const DIv& x, int64_t p) {
  static Z nl, nh;
  int64_t el = x.h.e, eh = x.l.e;
  mpz_neg(nl.z, x.h.m.z);
  mpz_neg(nh.z, x.l.m.z);
  dn(out.l, nl.z, el, p);
  up(out.h, nh.z, eh, p);
}

// DI.__mul__: 네 끝점 곱 t0..t3 중 (원본의 순서대로) 최소·최대 선택, dn/up 한 번씩
void di_mul(DIv& out, const DIv& x, const DIv& y, int64_t p) {
  static Z T[4];
  int64_t te[4];
  mulx(T[0].z, te[0], x.l, y.l);  // t0 = a*c
  mulx(T[1].z, te[1], x.l, y.h);  // t1 = a*d
  mulx(T[2].z, te[2], x.h, y.l);  // t2 = b*c
  mulx(T[3].z, te[3], x.h, y.h);  // t3 = b*d
  int lo = 0, hi = 0;
  for (int i = 1; i < 4; ++i) {
    if (vcmp(T[i].z, te[i], T[lo].z, te[lo]) < 0) lo = i;  // if _lt(t, lo): lo = t
    if (vcmp(T[hi].z, te[hi], T[i].z, te[i]) < 0) hi = i;  // if _lt(hi, t): hi = t
  }
  dn(out.l, T[lo].z, te[lo], p);
  up(out.h, T[hi].z, te[hi], p);
}

// DI.__pow__(2) (= sqr)
void di_sqr(DIv& out, const DIv& x, int64_t p) {
  static Z a, b;
  int64_t ea, eb;
  if (mpz_sgn(x.l.m.z) >= 0) {  // (n & 1) or lo[0] >= 0  (n = 2)
    powx2(a.z, ea, x.l);
    powx2(b.z, eb, x.h);
    dn(out.l, a.z, ea, p);
    up(out.h, b.z, eb, p);
    return;
  }
  if (mpz_sgn(x.h.m.z) <= 0) {  // hi[0] <= 0
    powx2(a.z, ea, x.h);
    powx2(b.z, eb, x.l);
    dn(out.l, a.z, ea, p);
    up(out.h, b.z, eb, p);
    return;
  }
  powx2(a.z, ea, x.l);  // a = lo^2
  powx2(b.z, eb, x.h);  // b = hi^2
  bool b_is_top = vcmp(a.z, ea, b.z, eb) < 0;  // top = b if _lt(a, b) else a
  if (b_is_top)
    up(out.h, b.z, eb, p);
  else
    up(out.h, a.z, ea, p);
  mpz_set_ui(out.l.m.z, 0);  // (0, 0)
  out.l.e = 0;
}

// dyad._qdn(a, b, p): b 가수 > 0
void qdn(EP& out, const EP& a, const EP& b, int64_t p) {
  static Z t, q;
  if (mpz_sgn(a.m.z) == 0) {
    mpz_set_ui(out.m.z, 0);
    out.e = 0;
    return;
  }
  int64_t k = p + 2 + bitlen(b.m.z) - bitlen(a.m.z);
  if (k < 0) k = 0;
  shift_guard(k);
  mpz_mul_2exp(t.z, a.m.z, (mp_bitcnt_t)k);
  mpz_fdiv_q(q.z, t.z, b.m.z);  // (ma << k) // mb
  dn(out, q.z, esub(esub(a.e, b.e), k), p);
}

// dyad._qup(a, b, p): q, r = divmod(ma << k, mb); if r: q += 1  (mb > 0 -> ceil)
void qup(EP& out, const EP& a, const EP& b, int64_t p) {
  static Z t, q;
  if (mpz_sgn(a.m.z) == 0) {
    mpz_set_ui(out.m.z, 0);
    out.e = 0;
    return;
  }
  int64_t k = p + 2 + bitlen(b.m.z) - bitlen(a.m.z);
  if (k < 0) k = 0;
  shift_guard(k);
  mpz_mul_2exp(t.z, a.m.z, (mp_bitcnt_t)k);
  mpz_cdiv_q(q.z, t.z, b.m.z);
  up(out, q.z, esub(esub(a.e, b.e), k), p);
}

inline void ep_neg(EP& out, const EP& a) {
  mpz_neg(out.m.z, a.m.z);
  out.e = a.e;
}

// DI._div (self / o)
void di_div(DIv& out, const DIv& x, const DIv& y, int64_t p) {
  static EP lo, hi, nyl, nyh, nhi, nlo;
  int syl = mpz_sgn(y.l.m.z), syh = mpz_sgn(y.h.m.z);
  if (syl <= 0 && 0 <= syh) throw ZeroDivEx{};  // yl[0] <= 0 <= yh[0]
  if (syh < 0) {
    ep_neg(nyl, y.h);  // yl, yh = _neg(yh), _neg(yl)
    ep_neg(nyh, y.l);
    qdn(lo, x.l, mpz_sgn(x.l.m.z) >= 0 ? nyh : nyl, p);
    qup(hi, x.h, mpz_sgn(x.h.m.z) >= 0 ? nyl : nyh, p);
    ep_neg(nhi, hi);  // [-hi, -lo]
    ep_neg(nlo, lo);
    dn(out.l, nhi.m.z, nhi.e, p);
    up(out.h, nlo.m.z, nlo.e, p);
    return;
  }
  qdn(lo, x.l, mpz_sgn(x.l.m.z) >= 0 ? y.h : y.l, p);
  qup(hi, x.h, mpz_sgn(x.h.m.z) >= 0 ? y.l : y.h, p);
  out.l = lo;
  out.h = hi;
}

// dyad._sqrt_parts -> (s, exact, e2)
void sqrt_parts(Z& s, bool& exact, int64_t& e2, const EP& a, int64_t p) {
  static Z m, n, r;
  mpz_set(m.z, a.m.z);
  int64_t e = a.e;
  if (e & 1) {  // 파이썬 e & 1 (음수도 2의 보수) 과 같다
    mpz_mul_2exp(m.z, m.z, 1);
    e = esub(e, 1);
  }
  int64_t num = 2 * (p + 2) - bitlen(m.z) + 1;
  int64_t k = num >= 0 ? num / 2 : -((-num + 1) / 2);  // 파이썬 // (floor)
  if (k < 0) k = 0;
  shift_guard(2 * k);
  mpz_mul_2exp(n.z, m.z, (mp_bitcnt_t)(2 * k));
  mpz_sqrtrem(s.z, r.z, n.z);  // isqrt
  exact = mpz_sgn(r.z) == 0;   // s*s == n
  e2 = esub(e / 2, k);         // e 는 짝수 -> 정확
}

void sqrt_dn(EP& out, const EP& a, int64_t p) {
  static Z s;
  if (mpz_sgn(a.m.z) == 0) {
    mpz_set_ui(out.m.z, 0);
    out.e = 0;
    return;
  }
  bool exact;
  int64_t e2;
  sqrt_parts(s, exact, e2, a, p);
  dn(out, s.z, e2, p);
}

void sqrt_up(EP& out, const EP& a, int64_t p) {
  static Z s;
  if (mpz_sgn(a.m.z) == 0) {
    mpz_set_ui(out.m.z, 0);
    out.e = 0;
    return;
  }
  bool exact;
  int64_t e2;
  sqrt_parts(s, exact, e2, a, p);
  if (!exact) mpz_add_ui(s.z, s.z, 1);
  up(out, s.z, e2, p);
}

// DI.sqrt
void di_sqrt(DIv& out, const DIv& x, int64_t p) {
  static EP lo, hi;
  if (mpz_sgn(x.l.m.z) < 0) throw SqrtNegEx{};
  sqrt_dn(lo, x.l, p);
  sqrt_up(hi, x.h, p);
  out.l = lo;
  out.h = hi;
}

// DI(int) = (_dn(n, 0, p), _up(n, 0, p))
void di_from_int(DIv& out, long long v, int64_t p) {
  static Z t;
  mpz_set_si(t.z, (long)v);
  dn(out.l, t.z, 0, p);
  up(out.h, t.z, 0, p);
}

// ---------------------------------------------------------------------------
// Val: DI 상수 또는 _GD (값 v + 기울기 목록 g, 항목 some=false 는 None)
// ---------------------------------------------------------------------------
struct GE {
  bool some = false;
  DIv d;
};
using GVec = std::vector<GE>;

struct Val {
  DIv v;
  std::shared_ptr<const GVec> g;  // null -> DI (상수)
  PyObject* src = nullptr;        // 바뀌지 않은 입력 객체 (출력 때 그대로 돌려줌; 참조는 KeepRefs 가 붙잡음)
  bool gd() const { return (bool)g; }
};

inline Val vconst(const DIv& d) {
  Val r;
  r.v = d;
  return r;
}

struct Engine {
  int64_t p;
  std::unordered_map<long long, DIv> icache;  // DI(n) 캐시 (노드 기반: 참조가 무효화되지 않음)

  explicit Engine(int64_t prec) : p(prec) {}

  const DIv& cint(long long n) {  // DI(n) at precision p (n >= 0)
    if (n < 0) throw Unsupported{"negative int constant"};
    auto it = icache.find(n);
    if (it != icache.end()) return it->second;
    DIv& d = icache[n];
    di_from_int(d, n, p);
    return d;
  }

  // ---- _GD / DI 연산 (파이썬 연산자 분배 그대로) ----

  // x + y
  Val add(const Val& x, const Val& y) {
    Val r;
    if (x.gd() && y.gd()) {  // _GD.__add__(_GD)
      di_add(r.v, x.v, y.v, p);
      const GVec& a = *x.g;
      const GVec& b = *y.g;
      auto g = std::make_shared<GVec>(a.size());
      for (size_t j = 0; j < a.size(); ++j) {
        if (!a[j].some)
          (*g)[j] = b[j];  // a None -> b (None 이든 DI 이든 그대로)
        else if (!b[j].some)
          (*g)[j] = a[j];
        else {
          (*g)[j].some = true;
          di_add((*g)[j].d, a[j].d, b[j].d, p);
        }
      }
      r.g = std::move(g);
    } else if (x.gd()) {  // _GD.__add__(DI): _GD(self.v + o, self.g)
      di_add(r.v, x.v, y.v, p);
      r.g = x.g;
    } else if (y.gd()) {  // DI + _GD -> _GD.__radd__(DI) = y.__add__(x)
      di_add(r.v, y.v, x.v, p);
      r.g = y.g;
    } else {
      di_add(r.v, x.v, y.v, p);
    }
    return r;
  }

  // x - y
  Val sub(const Val& x, const Val& y) {
    Val r;
    if (x.gd() && y.gd()) {  // _GD.__sub__(_GD)
      di_sub(r.v, x.v, y.v, p);
      const GVec& a = *x.g;
      const GVec& b = *y.g;
      auto g = std::make_shared<GVec>(a.size());
      for (size_t j = 0; j < a.size(); ++j) {
        if (!b[j].some)
          (*g)[j] = a[j];  // b None -> a
        else if (!a[j].some) {
          (*g)[j].some = true;  // a None -> -b
          di_neg((*g)[j].d, b[j].d, p);
        } else {
          (*g)[j].some = true;
          di_sub((*g)[j].d, a[j].d, b[j].d, p);
        }
      }
      r.g = std::move(g);
    } else if (x.gd()) {  // _GD.__sub__(DI): _GD(self.v - o, self.g)
      di_sub(r.v, x.v, y.v, p);
      r.g = x.g;
    } else if (y.gd()) {  // DI - _GD -> y.__rsub__(x): _GD(o - self.v, [None or -a])
      di_sub(r.v, x.v, y.v, p);
      r.g = neg_list(*y.g);
    } else {
      di_sub(r.v, x.v, y.v, p);
    }
    return r;
  }

  std::shared_ptr<const GVec> neg_list(const GVec& a) {
    auto g = std::make_shared<GVec>(a.size());
    for (size_t j = 0; j < a.size(); ++j) {
      if (a[j].some) {
        (*g)[j].some = true;
        di_neg((*g)[j].d, a[j].d, p);
      }
    }
    return g;
  }

  // -x
  Val neg(const Val& x) {
    Val r;
    di_neg(r.v, x.v, p);
    if (x.gd()) r.g = neg_list(*x.g);  // _GD(-self.v, [None or -a])
    return r;
  }

  // 목록의 각 항목에 DI 상수 c 를 곱함: [None or a * c]
  std::shared_ptr<const GVec> mul_list(const GVec& a, const DIv& c) {
    auto g = std::make_shared<GVec>(a.size());
    for (size_t j = 0; j < a.size(); ++j) {
      if (a[j].some) {
        (*g)[j].some = true;
        di_mul((*g)[j].d, a[j].d, c, p);
      }
    }
    return g;
  }

  // [None or a / c]
  std::shared_ptr<const GVec> div_list(const GVec& a, const DIv& c) {
    auto g = std::make_shared<GVec>(a.size());
    for (size_t j = 0; j < a.size(); ++j) {
      if (a[j].some) {
        (*g)[j].some = true;
        di_div((*g)[j].d, a[j].d, c, p);
      }
    }
    return g;
  }

  // x * y
  Val mul(const Val& x, const Val& y) {
    Val r;
    if (x.gd() && y.gd()) {  // _GD.__mul__(_GD)
      const DIv& av = x.v;
      const DIv& bv = y.v;
      const GVec& a = *x.g;
      const GVec& b = *y.g;
      auto g = std::make_shared<GVec>(a.size());
      static DIv t1, t2;
      for (size_t j = 0; j < a.size(); ++j) {
        if (!a[j].some) {
          if (b[j].some) {  // None if b is None else av * b
            (*g)[j].some = true;
            di_mul((*g)[j].d, av, b[j].d, p);
          }
        } else if (!b[j].some) {  // bv * a
          (*g)[j].some = true;
          di_mul((*g)[j].d, bv, a[j].d, p);
        } else {  // av * b + bv * a
          di_mul(t1, av, b[j].d, p);
          di_mul(t2, bv, a[j].d, p);
          (*g)[j].some = true;
          di_add((*g)[j].d, t1, t2, p);
        }
      }
      di_mul(r.v, av, bv, p);
      r.g = std::move(g);
    } else if (x.gd()) {  // _GD.__mul__(DI/int): _GD(self.v * o, [None or a * o])
      di_mul(r.v, x.v, y.v, p);
      r.g = mul_list(*x.g, y.v);
    } else if (y.gd()) {  // DI * _GD -> y.__rmul__(x) = y.__mul__(x)
      di_mul(r.v, y.v, x.v, p);
      r.g = mul_list(*y.g, x.v);
    } else {
      di_mul(r.v, x.v, y.v, p);
    }
    return r;
  }

  // x / y
  Val div(const Val& x, const Val& y) {
    Val r;
    if (x.gd() && y.gd()) {  // _GD.__truediv__(_GD)
      const DIv& bv = y.v;
      DIv q;
      di_div(q, x.v, bv, p);  // q = self.v / bv
      const GVec& a = *x.g;
      const GVec& b = *y.g;
      auto g = std::make_shared<GVec>(a.size());
      static DIv t1, t2;
      for (size_t j = 0; j < a.size(); ++j) {
        if (!b[j].some) {
          if (a[j].some) {  // None if a is None else a / bv
            (*g)[j].some = true;
            di_div((*g)[j].d, a[j].d, bv, p);
          }
        } else if (!a[j].some) {  // -((q * b) / bv)
          di_mul(t1, q, b[j].d, p);
          di_div(t2, t1, bv, p);
          (*g)[j].some = true;
          di_neg((*g)[j].d, t2, p);
        } else {  // (a - q * b) / bv
          di_mul(t1, q, b[j].d, p);
          di_sub(t2, a[j].d, t1, p);
          (*g)[j].some = true;
          di_div((*g)[j].d, t2, bv, p);
        }
      }
      r.v = q;
      r.g = std::move(g);
    } else if (x.gd()) {  // _GD.__truediv__(DI/int): _GD(self.v / o, [None or a / o])
      di_div(r.v, x.v, y.v, p);
      r.g = div_list(*x.g, y.v);
    } else if (y.gd()) {  // 상수 / _GD -> y.__rtruediv__(x)
      const DIv& u = y.v;
      DIv q, qu;
      di_div(q, x.v, u, p);   // q = o / u
      di_div(qu, q, u, p);    // qu = q / u
      const GVec& a = *y.g;
      auto g = std::make_shared<GVec>(a.size());
      static DIv t1;
      for (size_t j = 0; j < a.size(); ++j) {
        if (a[j].some) {  // -(qu * a)
          di_mul(t1, qu, a[j].d, p);
          (*g)[j].some = true;
          di_neg((*g)[j].d, t1, p);
        }
      }
      r.v = q;
      r.g = std::move(g);
    } else {
      di_div(r.v, x.v, y.v, p);
    }
    return r;
  }

  // x.sqr()
  Val sqr(const Val& x) {
    Val r;
    if (x.gd()) {  // _GD(av.sqr(), [None or (av * a) * 2])
      const DIv& av = x.v;
      const DIv& two = cint(2);
      const GVec& a = *x.g;
      auto g = std::make_shared<GVec>(a.size());
      static DIv t1;
      for (size_t j = 0; j < a.size(); ++j) {
        if (a[j].some) {
          di_mul(t1, av, a[j].d, p);
          (*g)[j].some = true;
          di_mul((*g)[j].d, t1, two, p);
        }
      }
      di_sqr(r.v, av, p);
      r.g = std::move(g);
    } else {
      di_sqr(r.v, x.v, p);
    }
    return r;
  }

  // x.sqrt()
  Val sqrt(const Val& x) {
    Val r;
    if (x.gd()) {  // r = v.sqrt(); r2 = r * 2; [None or a / r2]
      DIv rv, r2;
      di_sqrt(rv, x.v, p);
      di_mul(r2, rv, cint(2), p);
      r.g = div_list(*x.g, r2);
      r.v = rv;
    } else {
      di_sqrt(r.v, x.v, p);
    }
    return r;
  }

  // ---- taylor_step 도우미 ----

  // _div_int(a, n): n == 1 이면 a 그대로 (반올림 없음, 같은 객체)
  Val div_int(const Val& a, long long n) {
    if (n == 1) return a;
    return div(a, vconst(cint(n)));
  }

  // _conv(a, b, k)
  Val conv(const std::vector<Val>& a, const std::vector<Val>& b, long long k) {
    Val acc = mul(a[0], b[k]);
    for (long long j = 1; j <= k; ++j) acc = add(acc, mul(a[j], b[k - j]));
    return acc;
  }

  // _sq_conv(a, k)
  Val sq_conv(const std::vector<Val>& a, long long k) {
    bool have = false;
    Val acc;
    for (long long j = 0; j < (k + 1) / 2; ++j) {
      Val t = mul(a[j], a[k - j]);
      if (!have) {
        acc = std::move(t);
        have = true;
      } else {
        acc = add(acc, t);
      }
    }
    if (have) acc = mul(acc, vconst(cint(2)));
    if (k % 2 == 0) {
      Val m = sqr(a[k / 2]);
      if (!have) {
        acc = std::move(m);
        have = true;
      } else {
        acc = add(acc, m);
      }
    }
    return acc;
  }

  // _pow_m32_next(s, w, k)
  Val pow_m32_next(const std::vector<Val>& s, const std::vector<Val>& w, long long k) {
    const Val& s0 = s[0];
    if (k == 0) {
      if (mpz_sgn(s0.v.l.m.z) <= 0) throw SingularEx{};  // s0.lo_q <= 0
      Val sq = sqrt(s0);
      Val prod = mul(s0, sq);
      return div(vconst(cint(1)), prod);  // 1 / (s0 * s0.sqrt())
    }
    Val acc = mul(mul(s[1], w[k - 1]), vconst(cint(1 + 2 * k)));
    for (long long i = 2; i <= k; ++i) acc = add(acc, mul(mul(s[i], w[k - i]), vconst(cint(i + 2 * k))));
    return neg(div(acc, mul(s0, vconst(cint(2 * k)))));
  }

  // make_nbody(masses).coeffs(X, K) 본체. X 길이 4n, M 길이 n.
  std::vector<std::vector<Val>> nbody(const std::vector<Val>& X, const std::vector<DIv>& M, long long K) {
    size_t n = M.size();
    std::vector<std::pair<size_t, size_t>> pairs;
    for (size_t a = 0; a < n; ++a)
      for (size_t b = a + 1; b < n; ++b) pairs.emplace_back(a, b);
    std::vector<Val> Mv;
    for (auto& m : M) Mv.push_back(vconst(m));
    std::vector<std::vector<Val>> pos(2 * n), vel(2 * n);
    for (size_t i = 0; i < 2 * n; ++i) {
      pos[i].push_back(X[i]);
      vel[i].push_back(X[2 * n + i]);
    }
    size_t np = pairs.size();
    std::vector<std::vector<Val>> dx(np), dy(np), s(np), w(np);
    DIv zero;
    di_from_int(zero, 0, p);  // DI(0)
    for (long long k = 0; k < K; ++k) {
      std::vector<Val> acc(2 * n, vconst(zero));
      for (size_t pi = 0; pi < np; ++pi) {
        size_t a = pairs[pi].first, b = pairs[pi].second;
        dx[pi].push_back(sub(pos[2 * b][k], pos[2 * a][k]));
        dy[pi].push_back(sub(pos[2 * b + 1][k], pos[2 * a + 1][k]));
        {
          Val sx = sq_conv(dx[pi], k);
          Val sy = sq_conv(dy[pi], k);
          s[pi].push_back(add(sx, sy));
        }
        w[pi].push_back(pow_m32_next(s[pi], w[pi], k));
        Val ux = conv(dx[pi], w[pi], k);
        Val uy = conv(dy[pi], w[pi], k);
        acc[2 * a] = add(acc[2 * a], mul(Mv[b], ux));
        acc[2 * a + 1] = add(acc[2 * a + 1], mul(Mv[b], uy));
        acc[2 * b] = sub(acc[2 * b], mul(Mv[a], ux));
        acc[2 * b + 1] = sub(acc[2 * b + 1], mul(Mv[a], uy));
      }
      for (size_t i = 0; i < 2 * n; ++i) {
        pos[i].push_back(div_int(vel[i][k], k + 1));
        vel[i].push_back(div_int(acc[i], k + 1));
      }
    }
    std::vector<std::vector<Val>> out;
    out.reserve(4 * n);
    for (auto& v : pos) out.push_back(std::move(v));
    for (auto& v : vel) out.push_back(std::move(v));
    return out;
  }
};

// ---------------------------------------------------------------------------
// 파이썬 연결
// ---------------------------------------------------------------------------
PyObject* g_DI_type = nullptr;  // dyad.DI
PyObject* g_DI_new = nullptr;   // dyad.DI._new
PyObject* g_GD_type = nullptr;  // lohner_step._GD
PyObject* g_Unsupported = nullptr;
PyObject* g_Singular = nullptr;
PyObject* s_l = nullptr;
PyObject* s_h = nullptr;
PyObject* s_v = nullptr;
PyObject* s_g = nullptr;

struct KeepRefs {  // 계산·출력 변환 동안 입력 객체를 강한 참조로 붙잡아 둔다(빌린 참조만 쓰지 않음)
  std::vector<PyObject*> v;
  ~KeepRefs() {
    for (PyObject* o : v) Py_DECREF(o);
  }
};

struct Ref {  // 소유 참조
  PyObject* o;
  explicit Ref(PyObject* x = nullptr) : o(x) {}
  ~Ref() { Py_XDECREF(o); }
  Ref(const Ref&) = delete;
  Ref& operator=(const Ref&) = delete;
  PyObject* release() {
    PyObject* t = o;
    o = nullptr;
    return t;
  }
};

// 파이썬 int -> mpz (정확)
void py_to_mpz(PyObject* o, mpz_t out) {
  int overflow = 0;
  long v = PyLong_AsLongAndOverflow(o, &overflow);
  if (!overflow) {
    if (v == -1 && PyErr_Occurred()) throw PyErrSet{};
    mpz_set_si(out, v);
    return;
  }
  Py_ssize_t need = PyLong_AsNativeBytes(o, nullptr, 0, Py_ASNATIVEBYTES_LITTLE_ENDIAN);
  if (need < 0) throw PyErrSet{};
  if (need == 0) need = 1;
  std::vector<unsigned char> buf((size_t)need);
  Py_ssize_t got = PyLong_AsNativeBytes(o, buf.data(), need, Py_ASNATIVEBYTES_LITTLE_ENDIAN);
  if (got < 0) throw PyErrSet{};
  if (got > need) throw Unsupported{"int conversion size"};
  mpz_import(out, (size_t)need, -1, 1, -1, 0, buf.data());  // 크기 그대로(부호 없는 해석)
  if (buf[(size_t)need - 1] & 0x80) {                        // 2의 보수 음수 -> 2^(8n) 빼기
    static Z t;
    mpz_set_ui(t.z, 0);
    mpz_setbit(t.z, (mp_bitcnt_t)(8 * need));
    mpz_sub(out, out, t.z);
  }
}

// mpz -> 파이썬 int (정확)
PyObject* mpz_to_py(const mpz_t m) {
  if (mpz_fits_slong_p(m)) return PyLong_FromLong(mpz_get_si(m));
  size_t nbits = mpz_sizeinbase(m, 2);
  size_t nbytes = nbits / 8 + 1;  // 부호 비트 자리 포함
  std::vector<unsigned char> buf(nbytes, 0);
  size_t count = 0;
  mpz_export(buf.data(), &count, -1, 1, -1, 0, m);  // |m|, 작은 쪽부터
  if (count > nbytes) throw Unsupported{"mpz export size"};
  if (mpz_sgn(m) < 0) {  // 2의 보수
    unsigned carry = 1;
    for (size_t i = 0; i < nbytes; ++i) {
      unsigned v = (unsigned)(unsigned char)(~buf[i]) + carry;
      buf[i] = (unsigned char)(v & 0xff);
      carry = v >> 8;
    }
  }
  PyObject* r = PyLong_FromNativeBytes(buf.data(), nbytes, Py_ASNATIVEBYTES_LITTLE_ENDIAN);
  if (!r) throw PyErrSet{};
  return r;
}

void parse_ep(PyObject* t, EP& out) {
  if (!PyTuple_CheckExact(t) || PyTuple_GET_SIZE(t) != 2) throw Unsupported{"endpoint is not a 2-tuple"};
  PyObject* m = PyTuple_GET_ITEM(t, 0);
  PyObject* e = PyTuple_GET_ITEM(t, 1);
  if (!PyLong_CheckExact(m) || !PyLong_CheckExact(e)) throw Unsupported{"endpoint entries are not exact int"};
  int overflow = 0;
  long long ev = PyLong_AsLongLongAndOverflow(e, &overflow);
  if (overflow) throw Unsupported{"exponent too large"};
  if (ev == -1 && PyErr_Occurred()) throw PyErrSet{};
  if (ev > EXP_LIMIT || ev < -EXP_LIMIT) throw Unsupported{"exponent too large"};
  py_to_mpz(m, out.m.z);
  out.e = ev;
}

void parse_di(PyObject* o, DIv& out) {
  if (!Py_IS_TYPE(o, (PyTypeObject*)g_DI_type)) throw Unsupported{"not an exact dyad.DI"};
  Ref l(PyObject_GetAttr(o, s_l));
  if (!l.o) throw PyErrSet{};
  Ref h(PyObject_GetAttr(o, s_h));
  if (!h.o) throw PyErrSet{};
  parse_ep(l.o, out.l);
  parse_ep(h.o, out.h);
}

// X 의 원소 하나: DI 또는 _GD. L 은 처음 만난 _GD 의 g 길이(-1 = 아직 없음).
void parse_val(PyObject* o, Val& out, Py_ssize_t& L) {
  out.src = o;
  if (Py_IS_TYPE(o, (PyTypeObject*)g_DI_type)) {
    parse_di(o, out.v);
    return;
  }
  if (!Py_IS_TYPE(o, (PyTypeObject*)g_GD_type)) throw Unsupported{"element is neither dyad.DI nor lohner_step._GD"};
  Ref v(PyObject_GetAttr(o, s_v));
  if (!v.o) throw PyErrSet{};
  Ref g(PyObject_GetAttr(o, s_g));
  if (!g.o) throw PyErrSet{};
  parse_di(v.o, out.v);
  if (!PyList_CheckExact(g.o) && !PyTuple_CheckExact(g.o)) throw Unsupported{"_GD.g is not a list/tuple"};
  Py_ssize_t n = PySequence_Fast_GET_SIZE(g.o);
  if (L < 0) L = n;
  if (n != L) throw Unsupported{"_GD.g lengths differ"};
  PyObject** items = PySequence_Fast_ITEMS(g.o);
  auto gv = std::make_shared<GVec>((size_t)n);
  for (Py_ssize_t j = 0; j < n; ++j) {
    if (items[j] == Py_None) continue;
    (*gv)[(size_t)j].some = true;
    parse_di(items[j], (*gv)[(size_t)j].d);
  }
  out.g = std::move(gv);
}

PyObject* ep_to_py(const EP& x) {
  Ref m(mpz_to_py(x.m.z));
  Ref e(PyLong_FromLongLong((long long)x.e));
  if (!e.o) throw PyErrSet{};
  PyObject* t = PyTuple_New(2);
  if (!t) throw PyErrSet{};
  PyTuple_SET_ITEM(t, 0, m.release());
  PyTuple_SET_ITEM(t, 1, e.release());
  return t;
}

PyObject* di_to_py(const DIv& d) {
  Ref l(ep_to_py(d.l));
  Ref h(ep_to_py(d.h));
  PyObject* r = PyObject_CallFunctionObjArgs(g_DI_new, l.o, h.o, nullptr);  // DI._new(l, h)
  if (!r) throw PyErrSet{};
  return r;
}

PyObject* val_to_py(const Val& v) {
  if (v.src) {
    Py_INCREF(v.src);
    return v.src;
  }
  Ref d(di_to_py(v.v));
  if (!v.gd()) return d.release();
  const GVec& g = *v.g;
  Ref lst(PyList_New((Py_ssize_t)g.size()));
  if (!lst.o) throw PyErrSet{};
  for (size_t j = 0; j < g.size(); ++j) {
    PyObject* it;
    if (g[j].some) {
      it = di_to_py(g[j].d);
    } else {
      Py_INCREF(Py_None);
      it = Py_None;
    }
    PyList_SET_ITEM(lst.o, (Py_ssize_t)j, it);
  }
  PyObject* r = PyObject_CallFunctionObjArgs(g_GD_type, d.o, lst.o, nullptr);  // _GD(v, g)
  if (!r) throw PyErrSet{};
  return r;
}

bool ready() {
  if (!g_DI_type || !g_GD_type || !g_DI_new) {
    PyErr_SetString(PyExc_RuntimeError, "_fastcore.setup(DI, _GD) was not called");
    return false;
  }
  return true;
}

// C++ 예외 -> 파이썬 예외
void translate(const Unsupported& u) { PyErr_SetString(g_Unsupported, u.why.c_str()); }

#define FC_CATCH_ALL                                                                         \
  catch (const Unsupported& u) {                                                            \
    translate(u);                                                                           \
    return nullptr;                                                                         \
  }                                                                                         \
  catch (const SingularEx&) {                                                               \
    PyErr_SetNone(g_Singular);                                                              \
    return nullptr;                                                                         \
  }                                                                                         \
  catch (const ZeroDivEx&) {                                                                \
    PyErr_SetString(PyExc_ZeroDivisionError, "divisor interval contains 0");                \
    return nullptr;                                                                         \
  }                                                                                         \
  catch (const SqrtNegEx&) {                                                                \
    PyErr_SetString(PyExc_ValueError, "sqrt of interval with negative lower bound");        \
    return nullptr;                                                                         \
  }                                                                                         \
  catch (const PyErrSet&) {                                                                 \
    if (!PyErr_Occurred()) PyErr_SetString(PyExc_SystemError, "fastcore: missing error");   \
    return nullptr;                                                                         \
  }                                                                                         \
  catch (const std::bad_alloc&) {                                                           \
    PyErr_NoMemory();                                                                       \
    return nullptr;                                                                         \
  }                                                                                         \
  catch (const std::exception& ex) {                                                        \
    PyErr_SetString(PyExc_RuntimeError, ex.what());                                         \
    return nullptr;                                                                         \
  }

int64_t parse_prec(PyObject* o) {
  if (!PyLong_CheckExact(o)) throw Unsupported{"prec is not an exact int"};
  int overflow = 0;
  long long p = PyLong_AsLongLongAndOverflow(o, &overflow);
  if (overflow || p < 1 || p > PREC_LIMIT) throw Unsupported{"prec out of range"};
  return (int64_t)p;
}

// setup(DI, _GD)
PyObject* fc_setup(PyObject*, PyObject* args) {
  PyObject *di, *gd;
  if (!PyArg_ParseTuple(args, "OO", &di, &gd)) return nullptr;
  if (!PyType_Check(di) || !PyType_Check(gd)) {
    PyErr_SetString(PyExc_TypeError, "setup expects (DI class, _GD class)");
    return nullptr;
  }
  PyObject* dn_ = PyObject_GetAttrString(di, "_new");
  if (!dn_) return nullptr;
  Py_XDECREF(g_DI_type);
  Py_XDECREF(g_GD_type);
  Py_XDECREF(g_DI_new);
  Py_INCREF(di);
  Py_INCREF(gd);
  g_DI_type = di;
  g_GD_type = gd;
  g_DI_new = dn_;
  Py_RETURN_NONE;
}

// nbody_coeffs(X: list, K: int, masses: list[DI], prec: int) -> list[list]
PyObject* fc_nbody_coeffs(PyObject*, PyObject* args) {
  PyObject *X, *Ko, *Mo, *Po;
  if (!PyArg_ParseTuple(args, "OOOO", &X, &Ko, &Mo, &Po)) return nullptr;
  if (!ready()) return nullptr;
  try {
    if (!PyList_CheckExact(X) && !PyTuple_CheckExact(X)) throw Unsupported{"X is not a list/tuple"};
    if (!PyList_CheckExact(Mo)) throw Unsupported{"masses is not a list"};
    if (!PyLong_Check(Ko)) throw Unsupported{"K is not an int"};
    int overflow = 0;
    long long K = PyLong_AsLongLongAndOverflow(Ko, &overflow);
    if (K == -1 && PyErr_Occurred()) throw PyErrSet{};
    if (overflow || K > K_LIMIT) throw Unsupported{"K out of range"};
    int64_t p = parse_prec(Po);
    Py_ssize_t n = PyList_GET_SIZE(Mo);
    if (n < 1) throw Unsupported{"no bodies"};
    if (PySequence_Fast_GET_SIZE(X) != 4 * n) throw Unsupported{"len(X) != 4 * n"};
    std::vector<DIv> M((size_t)n);
    for (Py_ssize_t i = 0; i < n; ++i) parse_di(PyList_GET_ITEM(Mo, i), M[(size_t)i]);
    Py_ssize_t L = -1;
    std::vector<Val> Xv((size_t)(4 * n));
    PyObject** items = PySequence_Fast_ITEMS(X);
    KeepRefs keep;  // Val.src 가 가리키는 입력 객체를 끝까지 살려 둔다
    keep.v.reserve((size_t)(4 * n));
    for (Py_ssize_t i = 0; i < 4 * n; ++i) {
      Py_INCREF(items[i]);
      keep.v.push_back(items[i]);
    }
    for (Py_ssize_t i = 0; i < 4 * n; ++i) parse_val(keep.v[(size_t)i], Xv[(size_t)i], L);
    Engine eng(p);
    std::vector<std::vector<Val>> res = eng.nbody(Xv, M, K);
    Ref out(PyList_New((Py_ssize_t)res.size()));
    if (!out.o) throw PyErrSet{};
    for (size_t i = 0; i < res.size(); ++i) {
      Ref row(PyList_New((Py_ssize_t)res[i].size()));
      if (!row.o) throw PyErrSet{};
      for (size_t k = 0; k < res[i].size(); ++k) PyList_SET_ITEM(row.o, (Py_ssize_t)k, val_to_py(res[i][k]));
      PyList_SET_ITEM(out.o, (Py_ssize_t)i, row.release());
    }
    return out.release();
  }
  FC_CATCH_ALL
}

// 시험용: val_op(op, x, y, prec). x, y 는 DI, _GD, int(>= 0, 상수 DI(n)) 또는 None(단항).
// 파이썬 연산자 x OP y 와 같은 결과를 내야 한다(test_fc_ops.py 가 비교).
PyObject* fc_val_op(PyObject*, PyObject* args) {
  const char* op;
  PyObject *xo, *yo, *Po;
  if (!PyArg_ParseTuple(args, "sOOO", &op, &xo, &yo, &Po)) return nullptr;
  if (!ready()) return nullptr;
  try {
    int64_t p = parse_prec(Po);
    Engine eng(p);
    Py_ssize_t L = -1;
    auto conv_arg = [&](PyObject* o, Val& v) {
      if (PyLong_CheckExact(o)) {
        long long n = PyLong_AsLongLong(o);
        if (n == -1 && PyErr_Occurred()) throw PyErrSet{};
        v = vconst(eng.cint(n));
      } else {
        parse_val(o, v, L);
        v.src = nullptr;
      }
    };
    Val x, y, r;
    conv_arg(xo, x);
    std::string s(op);
    if (s == "neg" || s == "sqr" || s == "sqrt") {
      if (s == "neg") r = eng.neg(x);
      else if (s == "sqr") r = eng.sqr(x);
      else r = eng.sqrt(x);
    } else {
      conv_arg(yo, y);
      if (s == "add") r = eng.add(x, y);
      else if (s == "sub") r = eng.sub(x, y);
      else if (s == "mul") r = eng.mul(x, y);
      else if (s == "div") r = eng.div(x, y);
      else throw Unsupported{"unknown op"};
    }
    return val_to_py(r);
  }
  FC_CATCH_ALL
}

// 시험용: int 왕복 변환 (py -> mpz -> py)
PyObject* fc_int_roundtrip(PyObject*, PyObject* arg) {
  try {
    if (!PyLong_CheckExact(arg)) throw Unsupported{"not an exact int"};
    Z z;
    py_to_mpz(arg, z.z);
    return mpz_to_py(z.z);
  }
  FC_CATCH_ALL
}

PyObject* fc_gmp_version(PyObject*, PyObject*) { return PyUnicode_FromString(gmp_version); }

PyMethodDef fc_methods[] = {
    {"setup", fc_setup, METH_VARARGS, "setup(DI, _GD): register the original classes"},
    {"nbody_coeffs", fc_nbody_coeffs, METH_VARARGS, "nbody_coeffs(X, K, masses_DI, prec)"},
    {"val_op", fc_val_op, METH_VARARGS, "val_op(op, x, y, prec) (test helper)"},
    {"int_roundtrip", fc_int_roundtrip, METH_O, "int_roundtrip(n) (test helper)"},
    {"gmp_version", fc_gmp_version, METH_NOARGS, "GMP version string"},
    {nullptr, nullptr, 0, nullptr}};

PyModuleDef fc_module = {PyModuleDef_HEAD_INIT, "_fastcore", "Gate FC1 C++ core", -1, fc_methods,
                         nullptr, nullptr, nullptr, nullptr};

}  // namespace

PyMODINIT_FUNC PyInit__fastcore(void) {
  PyObject* m = PyModule_Create(&fc_module);
  if (!m) return nullptr;
  g_Unsupported = PyErr_NewException("_fastcore.Unsupported", nullptr, nullptr);
  g_Singular = PyErr_NewException("_fastcore.Singular", nullptr, nullptr);
  s_l = PyUnicode_InternFromString("_l");
  s_h = PyUnicode_InternFromString("_h");
  s_v = PyUnicode_InternFromString("v");
  s_g = PyUnicode_InternFromString("g");
  if (!g_Unsupported || !g_Singular || !s_l || !s_h || !s_v || !s_g) {
    Py_DECREF(m);
    return nullptr;
  }
  Py_INCREF(g_Unsupported);
  Py_INCREF(g_Singular);
  if (PyModule_AddObject(m, "Unsupported", g_Unsupported) < 0 ||
      PyModule_AddObject(m, "Singular", g_Singular) < 0) {
    Py_DECREF(m);
    return nullptr;
  }
  return m;
}
