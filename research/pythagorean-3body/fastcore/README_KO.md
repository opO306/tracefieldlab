# fastcore — Gate FC1: 테일러 계수 계산을 C++로 (원본과 비트 단위 동일)

계획: `../GATE_FC1_PLAN_KO.md` (SHA-256 `37276492…552ed`, 봉인 `../GATE_FC1_SEAL.json`).
결과 보고: `FC1_IMPL_RESULT_KO.md`.

## 0. 쉬운 설명

- **비유**: 원본 계산은 부품(128비트 수) 하나 옮길 때마다 서류 작업(파이썬 처리)이 붙는 공장이다.
  가장 바쁜 공정(테일러 계수 계산)만 **같은 작업 순서를 그대로 따르는 기계(C++)** 로 바꿨다.
- **실제 의미**: `taylor_step.make_nbody(masses).coeffs(X, K)` 와 같은 일을 하는 함수를 C++ 로 다시 썼다.
  구간 끝점 하나하나를 원본과 **같은 정수 연산 순서, 같은 방향 반올림**으로 계산하므로, 결과 끝점의
  (가수, 지수) 정수 쌍이 원본과 같아야 한다(같지 않으면 불합격 — 계획 4절).
- 원본(`../rigorous/`)은 한 글자도 바꾸지 않았다. fastcore 는 원본을 **읽기(import)만** 한다.
  원본은 fastcore 를 모른다. 새 실행 스크립트 `r5_run.py` 를 쓸 때만 C++ 경로가 켜진다.

## 1. 파일

| 파일 | 내용 |
|---|---|
| `fastcore.cpp` | C++17 + GMP(mpz) + CPython C API. DI 끝점 연산, DI 연산, `_GD` 연산, N체 계수 재귀, 파이썬 변환 |
| `build.sh` | 빌드 (`g++ -O2 -std=c++17 … -lgmp`), `build_info.txt` 기록 |
| `fast_system.py` | 어댑터 `make_nbody_fast(masses)` (원본 `_System` 객체를 돌려줌) |
| `r5_run.py` | 실행 감싸개: 원본 `r4_run.py → r3c_run.py` 를 그대로 돌리되 system 생성만 교체 |
| `check_seals.py` | 봉인 4개(R3c, R3c′, R4, FC1)의 파일 해시 재계산 |
| `fc_common.py` | 자체 시험 공통(비트 비교) |
| `test_fc_ops.py` | FC-0(계획 밖 추가): 연산 하나씩 비교 |
| `test_fc_equiv.py` | FC-1, FC-2, FC-3, FC-X(경계, 계획 밖 추가) |
| `run_selftests.sh` | FC-0~3, FC-X 를 2 프로세스로 실행 |
| `fc4_compare.py` | FC-4 걸음 기록 비교, FC-5 걸음당 시간 |
| `prof_fc.py` | 측정 보조(cProfile / 걸음 시간) |

## 2. 사용법 (Debian, Python 3.13)

```sh
sudo apt-get install -y --no-install-recommends g++ libgmp-dev python3-dev   # Debian 공식 저장소
sh build.sh                       # → _fastcore.cpython-313-x86_64-linux-gnu.so, build_info.txt
python3 -I -B check_seals.py      # 봉인 해시 확인 (ALL_OK)
sh run_selftests.sh <새 결과 폴더>   # FC-0~3, FC-X (2코어)
python3 -I -B r5_run.py --outdir <새 폴더, 이름에 fc> --prec 128 --pivot width [--max-steps N]
python3 -I -B fc4_compare.py --orig-dir <원본 out> --fc-dir <r5 outdir> --tag 128w [--orig-timing-dir D]
```

파이썬에서 직접:

```python
import fast_system as FS                 # rigorous/ 를 sys.path 에 넣고 원본 모듈을 import
import dyad, lohner_step as ls
dyad.set_prec(128)
sysf = FS.make_nbody_fast([3, 4, 5])     # taylor_step._System 그대로
c = sysf.coeffs(X, 20)                   # X: dyad.DI 목록 → DI 경로, _GD 목록 → GD 경로
FS.fc_stats(sysf)                        # {'fast_calls', 'fallback_calls', 'fallback_reasons', 'singular'}
```

- `fast_system` 을 import 하면 원본 세 파일(`dyad.py`, `taylor_step.py`, `lohner_step.py`)의 SHA-256 이
  FC1 봉인 값과 같은지 먼저 확인한다. 다르면 ImportError (C++ 은 그 판을 옮긴 것이므로).
- 원본 모듈이 다른 폴더에서 이미 import 돼 있으면 ImportError (DI 클래스가 둘이 되는 것을 막음).

## 3. 설계

### 3.1 표현
- 끝점 `(m, e)` = m·2^e. m 은 GMP `mpz`, e 는 int64(덧셈·곱 넘침 검사).
- 입력 끝점은 **받은 그대로** 둔다(표준형으로 바꾸지 않음). 원본도 입력 튜플을 그대로 쓴다
  (`_div_int(a, 1)` 은 a 를 그대로 돌려주고, `_qdn` 의 k 는 가수 비트 길이에 의존).
- `_GD` 는 값 구간 + 기울기 목록. 목록 항목마다 "None(정확한 0)" 표시를 따로 들고 다닌다.
  DI 상수와 `_GD` 를 실행 중에 구별해 원본의 연산자 분배(`__radd__`, `__rsub__`, `__rmul__`,
  `__rtruediv__`)를 그대로 따른다(1체·DI/_GD 혼합 입력도 원본과 같은 결과).

### 3.2 원본과 대응표

| 원본 (파이썬) | C++ (`fastcore.cpp`) | 비고 |
|---|---|---|
| `dyad._canon` | `canon` | mpz_scan1 (음수도 같은 최하위 비트 위치) |
| `dyad._dn` / `_up` | `dn` / `up` | `>>` = mpz_fdiv_q_2exp, `-((-m)>>s)` = mpz_cdiv_q_2exp |
| `dyad._lt` | `vcmp(...) < 0` | 정확한 값 비교(최상위 비트 위치 비교 후 짧은 정렬). 원본의 지수 맞춤 비교와 같은 값 |
| `dyad._addx` (+`_neg`) | `addx(…, negb)` | 0 지름길 두 개 그대로 |
| `dyad._mulx`, `_powx(·,2)` | `mulx`, `powx2` | 곱이 0 이면 (0,0) |
| `dyad._qdn` / `_qup` | `qdn` / `qup` | k = p+2+len(mb)−len(ma), 음수면 0. floor / ceil 정수 나눗셈 |
| `dyad._sqrt_parts`, `_sqrt_dn`, `_sqrt_up` | 같은 이름 | isqrt = mpz_sqrtrem, 홀수 지수 처리·`//` 내림 그대로 |
| `DI(int)` (`_q_floor_ceil` int 갈래) | `di_from_int` | 정수 상수 강제 변환(`_coerce`)에 씀 |
| `DI.__add__`, `__sub__`, `__neg__` | `di_add`, `di_sub`, `di_neg` | |
| `DI.__mul__` | `di_mul` | t0..t3 순서, `if _lt(t, lo)` / `if _lt(hi, t)` 두 비교 그대로 |
| `DI.__pow__(2)` = `sqr` | `di_sqr` | 세 갈래(lo≥0 / hi≤0 / 0 걸침) |
| `DI._div` | `di_div` | 제수 0 포함 → ZeroDivisionError(같은 메시지), 음 제수 갈래 |
| `DI.sqrt` | `di_sqrt` | 하한 음수 → ValueError(같은 메시지) |
| `_GD.__add__/__radd__` | `Engine::add` | None 규칙: a None→b, b None→a, 둘 다→a+b; GD+DI 는 g 공유 |
| `_GD.__sub__/__rsub__` | `Engine::sub` | b None→a, a None→−b, 둘 다→a−b; DI−GD → (o−v, [−a]) |
| `_GD.__neg__` | `Engine::neg` | |
| `_GD.__mul__/__rmul__` | `Engine::mul` | a None→(b None? None : av·b), b None→bv·a, 둘 다→av·b+bv·a; 상수는 [a·o] |
| `_GD.__truediv__` | `Engine::div` | q=v/bv 먼저; b None→a/bv, a None→−((q·b)/bv), 둘 다→(a−q·b)/bv |
| `_GD.__rtruediv__` | `Engine::div`(상수/GD) | q=o/u, qu=q/u, [−(qu·a)] |
| `_GD.sqr` / `_GD.sqrt` | `Engine::sqr` / `Engine::sqrt` | (av·a)·2 / a/(r·2) |
| `taylor_step._div_int` | `Engine::div_int` | n==1 이면 그대로(반올림 없음) |
| `taylor_step._conv` / `_sq_conv` | `Engine::conv` / `sq_conv` | 더하는 순서 그대로, 대칭쌍 2배 후 가운데 sqr |
| `taylor_step._pow_m32_next` | `Engine::pow_m32_next` | k=0: s0 하한 ≤ 0 → _Singular, 1/(s0·sqrt(s0)) |
| `make_nbody.coeffs` | `Engine::nbody` | 쌍 순서 (a<b 사전순), acc 갱신 순서, pos/vel 추가 순서 그대로 |
| `M = [DI(m) for m in masses]` | `fast_system.coeffs` | **원본 `dyad.DI` 로 파이썬에서** 호출 시점 정밀도로 만든 뒤 넘김 |

### 3.3 정밀도
- `coeffs` 호출 시점의 `dyad.get_prec()` 하나를 쓴다(R0 결정 A). 원본도 한 호출 안에서 `_PREC` 를 바꾸지 않는다.

### 3.4 예외
- 특이(어느 쌍의 s_[0] 하한 ≤ 0, k=0 에서 쌍 순서대로 검사) → C++ 이 신호 → 어댑터가 원본
  `taylor_step._Singular()` 를 던진다. K ≤ 0 이면 원본처럼 검사 자체를 하지 않는다.
- ZeroDivisionError / ValueError 는 원본과 같은 메시지(정상 입력의 재귀에서는 생기지 않음).

### 3.5 대체 경로(fallback)
- C++ 이 다루지 않는 입력: 원본 클래스가 아닌 원소(하위 클래스 포함), `_GD.g` 길이가 서로 다름,
  `len(X) ≠ 4N`, X 가 list/tuple 이 아님, K 가 int 가 아님, 입력 지수 |e| > 2^40, 왼쪽 밀기 > 2^26 비트,
  지수 넘침, 정밀도 > 2^24, K > 10^6.
- 이때 C++ 은 계산을 버리고 `Unsupported` 를 알리고, 어댑터가 **원본 파이썬 `coeffs` 로 다시 계산**한다
  (결과·예외는 원본 그 자체). 횟수는 `fc_stats(system)["fallback_calls"]` 와 `r5_run` 의 메타 파일에 남는다.
  `make_nbody_fast(masses, allow_fallback=False)` 이면 대신 `FastcoreUnsupported` 를 던진다(시험용).

## 4. 분리 규칙 (계획 3절)
- 소스: `fastcore/` 만. `rigorous/` 는 import(읽기)만.
- 실행: `r5_run.py` 만 C++ 경로를 쓴다. 결과는 새 폴더(이름에 `fc` 필수, `rigorous/` 안 거부),
  파일 이름 앞에 `fc_` (예: `fc_r3c_steps_128w.jsonl`). 덮어쓰기 거부는 원본 r3c_run/r4_run 검사를 그대로 씀.
- `r5_run.py` 가 바꾸는 것은 r3c_run 모듈의 `ts`(→ make_nbody 만 교체)와 r3c_run·r4_run 모듈의 `os`
  (→ 결과 파일 이름에 `fc_`)뿐이다. 보폭·QR·정확한 역행렬·포함 점검·격자 판정은 원본 그대로.

## 5. 한계·주의
- GIL 을 잡은 채 실행하고 정적 임시 변수를 쓴다: 한 스레드에서 부르는 것을 전제로 한다.
- Python ≥ 3.13 필요(`PyLong_AsNativeBytes` / `PyLong_FromNativeBytes`). 서버: Debian 13, Python 3.13.5,
  g++ 14.2.0, GMP 6.3.0 에서 빌드·시험.
- 남은 파이썬 부분(QR, 정확한 역행렬, 구간 행렬곱, 사전 상자 포함 검사)은 이번 범위가 아니다(FC2 후보).
- 증거 수준: 자체 시험 범위 안에서 "원본과 비트 단위 동일"(tested). 형식 증명이 아니며, 시험하지 않은
  입력에서 다를 가능성은 남는다(계획 7절). 최종 판정은 독립 검사(FC-6)와 조율 세션이 한다.
