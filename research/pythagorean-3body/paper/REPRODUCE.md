# 재현 꾸러미 (Reproduction package) — 피타고라스 삼체 문제 컴퓨터 보조 증명

> 작성 2026-10-12. 내부용. **외부 공개 전에 영어판이 필요하다 [TODO].**
> **공개 전제:** 저자(Jeomseon, 필명)는 현재 arXiv 추천인이 없다. 이 꾸러미는 **DOI가 붙는 공개 보관소(예: Zenodo)** 또는 **GitHub 공개 저장소**로 내는 것을 전제로 쓴다. 두 서비스의 업로드 크기 한도·라이선스 요구·DOI 발급 조건·GitHub–Zenodo 연동 방식은 **확인하지 않았다 [CHECK]**. 7절 참고.
> 경로는 모두 `research/pythagorean-3body/` 기준이다.
> 이 문서를 만들면서 코드는 실행하지 않았다. 한 일은 (1) 봉인 파일 `GATE_*_SEAL.json`과 보고서에서 해시를 옮겨 적기, (2) 노트북에서 실제 파일의 SHA-256을 다시 계산해 대조하기, (3) 서버 보관본 tgz의 목록 보기와 작은 요약 파일 몇 개 풀어 읽기(스크래치 폴더)뿐이다.

> 📦 **쉬운 설명**
> - **비유:** 요리책의 "재료 목록 + 순서 + 오븐 온도"다. 다른 사람이 같은 재료(파일)를 같은 순서로 넣으면 같은 결과가 나와야 한다. 재료가 바꿔치기되지 않았는지는 **지문(SHA-256 해시)**으로 확인한다.
> - **실제 의미:** SHA-256은 파일 내용에서 계산한 64자리 지문이다. 한 글자만 바뀌어도 완전히 달라진다. 실행 전에 봉인한 지문과 지금 파일의 지문이 같으면, 실행에 쓴 파일과 지금 파일이 같다는 뜻이다.

---

## 1. 무엇을 재현하나

| 단계 | 결과 | 사람이 읽을 보고서 |
|---|---|---|
| R4 | t ∈ [0, 61.0055] 엄밀 감싸기(128비트), t = 59.50~61.00 상자 302개 저장 | `GATE_R4_PLAN_KO.md`, `results_r4/r3c_final_128w.json` |
| Y1 | 상자 302개에 (★) 판정. t = 59.73~61.00 전부 PASS | `yoshida_check/Y1_RESULT_KO.md`, `Y1_FINAL_JUDGMENT_KO.md` |
| 검사기 | 구간 층(R0), 테일러 걸음(R1), Lohner 걸음(R3c v2)의 독립 검사 | `rigorous/R0_FINAL_RUN_KO.md`, `rigorous/R1_FINAL_RUN_KO.md`, `rigorous/R3C_CHECK_V2_FINAL_KO.md` |
| 증명 점검 | (★)·보조정리의 기호·수치 점검(증명을 대신하지 않음) | `proof_review/PROOF_REVIEW_KO.md` |

---

## 2. 파일과 SHA-256

"봉인값" = 실행 전에 기록한 해시(출처 열). "재계산" = 2026-10-12 노트북(Git Bash `sha256sum`)에서 다시 계산한 값. **아래 모든 행에서 둘이 같았다.**

### 2.1 적분기 (R4가 쓴 코드와 입력)

| 파일 | SHA-256 | 봉인 출처 | 재계산 |
|---|---|---|---|
| `rigorous/dyad.py` | `4d49cf8f707061d128fce3ed44b59ccb3a6bda14937e32dcd469cd201fbabaab` | `GATE_R4_SEAL.json` | 일치 |
| `rigorous/taylor_step.py` | `0183f0f8dc76dcf9c1a9933d8f16a636ac011bc15e9a95e0749fb2459be01c31` | `GATE_R4_SEAL.json` | 일치 |
| `rigorous/lohner_step.py` | `936138ea73b878c811df8db7b582ddb6883590be0f47cc7ece76d5ffc2db6b62` | `GATE_R4_SEAL.json` | 일치 |
| `rigorous/r3c_run.py` | `ce113414f17053d41d03a147a8298b3813e11f07cb8a5b1b1a711306f86dbe5b` | `GATE_R4_SEAL.json` | 일치 |
| `rigorous/r4_run.py` | `9b10b434e9b12e040091dcb9923e67c7e400f711259579b5855ae20ef2de93d9` | `GATE_R4_SEAL.json` | 일치 |
| `rigorous/lohner_step_width.py` (검사기용 감싸개, pivot=width 고정) | `8992d8b6c83dc3aea8c2b86413b06b914d5d8c13e177bc0add58386b349cc693` | `GATE_R3CP_SEAL.json` | 일치 |
| `results/ref_B_grid.json` (비엄밀 60자리 기준 궤적, 대조용) | `3bfe3e44024d90bfe495c365bebb6c3c2ea7a92c18869d177bd5e6beaf8d3478` | `GATE_R4_SEAL.json` | 일치 |

### 2.2 R4 결과 (Y1 입력)

| 파일 | SHA-256 | 봉인 출처 | 재계산 |
|---|---|---|---|
| `results_r4/r4_boxes_128w.jsonl` (302줄, 정확한 분수) | `f65101c194b6f9d4a67573773a9c517c49ccde8f59924954bd6598f6301d59f5` | `GATE_Y1_SEAL.json` | 일치 |
| `results_r4/r3c_final_128w.json` (실행 요약) | `7241e8c9788617ed882d0ec6ab23b173e4b2a1e93c2fa18a5e1965883827c1f4` | `GATE_Y1_SEAL.json` | 일치 |

### 2.3 계획과 봉인 (실행 전 고정 문서)

| 파일 | SHA-256 | 봉인 출처 | 재계산 |
|---|---|---|---|
| `GATE_R3CP_PLAN_KO.md` | `5b83160b45c75c6f1d46ee01b4f2bd4f35af45278ba2a0f458da5fb577b4ba4e` | `GATE_R3CP_SEAL.json` | 일치 |
| `GATE_R4_PLAN_KO.md` | `e59ce0ee1ce1a73eb22c5dc5b62f4878f522b54282d6edb37f059b64ac43acce` | `GATE_R4_SEAL.json` | 일치 |
| `GATE_R4_ADDENDUM_1_KO.md` (112·96비트 측정 실행 추가) | `75e167d26a15578294a64747231e73095deeda3e7763b3e4bea6057ff2be9a4b` | `GATE_R4_ADDENDUM_1_SEAL.json` | 일치 |
| `GATE_Y1_PLAN_KO.md` | `416b408ee41c63044b90122b1fee790ed14d826a01435c4870e85a706195cab5` | `GATE_Y1_SEAL.json` | 일치 |
| `YOSHIDA_REDERIVATION_KO.md` ((★) 원 진술) | `fab037d1c5ad124402fea6c9f6a7f733818de63d4af7b181625cc0f83c7f04a5` | `GATE_Y1_SEAL.json` | 일치 |
| `GATE_R3C_SEAL.json` | `2315f76ece4df82ad6227175c2a4af7df1bea5db9b2855d918a76c41ddcee6c1` | `GATE_R4_SEAL.json` | 일치 |
| `GATE_R3CP_SEAL.json` | `e4a12399c1ded092129f2b6818c44cda2fb5191bac70ad0ebb84576409cb8c1c` | `GATE_R4_SEAL.json` | 일치 |
| `GATE_R4_SEAL.json` | `bbfd40f4c163ad2684a14422b73c7491a6380c16891526f691978596fa672960` | `GATE_Y1_SEAL.json`, `GATE_FC1_SEAL.json` | 일치 |

### 2.4 Y1 판정기

| 파일 | SHA-256 | 기록 출처 | 재계산 |
|---|---|---|---|
| `yoshida_check/y1_exact.py` | `3ade72c85c0a681e0f28af5b94b3a9c36cb522be9aaf79c5f510eb8fbd262c63` | `Y1_RESULT_KO.md` 표(243행) | 일치 |
| `yoshida_check/y1_mpiv.py` | `a546ebb395b234b58ce14f3b6eb6be524726edbb39c5b77dcc5e381fa5008f8b` | 같음 | 일치 |
| `yoshida_check/y1_selftest.py` | `38ac3919d44ab39cc2f851535b02cc53e682b739da9df7032aec020dafb868d0` | 같음 | 일치 |
| `yoshida_check/y1_compare.py` | `06a54b5035046962d56087470138d5e33703557387a372ee19aa000107b44988` | 같음 | 일치 |
| `yoshida_check/y1_results_exact.json` | `a664bc87252cb1da4eccf73910d6488ef1f62c2289c1953a1c07ce86606fb2c7` | 같음 | 일치 |
| `yoshida_check/y1_results_mpiv.json` | `121f9e65205ce39ae3b8061e4a9b2f97945e0f2706abc1172af9bd73e1523ce7` | 같음 | 일치 |

주의: Y1 판정기 해시는 **실행 전 봉인이 아니라** 판정 세션이 실행 뒤 보고서에 적은 값이다(판정기는 Y1 계획 봉인 뒤에 만들어졌다).

### 2.5 독립 검사기

| 파일 | SHA-256 | 기록 출처 | 재계산 |
|---|---|---|---|
| `rigorous/check_dyad.py` (R0) | `7323e236a062420b48108216708c76c40bb4e266d5f46a138a59b4ab01d21f91` | `rigorous/R0_FINAL_RUN_KO.md` | 일치 |
| `rigorous/check_r1.py` (R1) | `13fe422df8427f58b34cdcf18fe3b6700796236370ce915610b5bce8d5b0592d` | `rigorous/R1_FINAL_RUN_KO.md` | 일치 |
| `rigorous/check_r3c.py` (R3c v1, 수정 안 함) | `031f11e28f63a2767b01c6e8a532983ee9882e8b82b7e02192eb8bbc7cbe45b3` | `GATE_R3CP_SEAL.json` | 일치 |
| `rigorous/check_r3c_v2.py` (R3c v2) | `5f13766f2a461f49c8b36e86456574bd38b5f432cd0c384ea0d2d6f4ab078e0c` | 공유 노트 r3c-check-v2, `R3C_CHECK_V2_FINAL_KO.md`(앞 8자리) | 일치 |
| `rigorous/R3C_CHECK_V2_lohner_step_width.json` (v2 결과, width) | `001ec9b4377ffce328f39e55159cfd684dfa2b66997012320508b35e0e38ba38` | `R3C_CHECK_V2_FINAL_KO.md`(앞 8자리) | 앞 8자리 일치 |
| `rigorous/R3C_CHECK_V2_lohner_step.json` (v2 결과, colnorm) | `014d7c5c0d5f8212b943f9c1818968edb5e7ad9655d265313f0078b535ea7c5c` | 같음 | 앞 8자리 일치 |
| `rigorous/R0_CHECK_RESULT_dyad.json` | `b8064bb50dff859b00ad58a3ab0b629ebc57e9e0cd5d66f11fcf3620f49274a7` | 기록 없음 | 이번에 처음 계산 |
| `rigorous/R1_CHECK_RESULT_taylor_step.json` | `7280e009f2c74b3fbd6fffc84f08ced5db2ed9f9ea0b4441a6bcdab70daf8231` | 기록 없음 | 이번에 처음 계산 |

### 2.6 증명 점검 스크립트

| 파일 | SHA-256 | 기록 출처 | 재계산 |
|---|---|---|---|
| `proof_review/symbolic_checks.py` | `0275b32a2ed24d57baebd3321bcd6bbacd3e4078ce3fa322e06095d73c749938` | 기록 없음 | 이번에 처음 계산 |
| `proof_review/lemma_adversarial.py` | `be533c373df4764a8168d4feac886ca6bebea4bc3570d3e172c4372661cdb5a8` | 기록 없음 | 이번에 처음 계산 |
| `proof_review/star_dynamic_falsify.py` | `dee5494d70ea533f881f425816cc2ba448e51663dc21ba47c072425fbb9c025d` | `PROOF_REVIEW_KO.md` §6.3(앞 8자리 `dee5494d`) | 앞 8자리 일치 |

---

## 3. 환경

| 항목 | 값 | 출처 |
|---|---|---|
| R4 실행 | Linux, Python 3.13.5 (GCC 14.2.0 빌드), 순수 파이썬(외부 라이브러리 없음), 실행당 코어 1개 | `results_r4/r3c_final_128w.json` "python" |
| 서버 | Google Cloud n2 인스턴스, CPU 12개. 기계 형식 n2-custom-12는 메모리 기록이며 **이번에 다시 확인하지 않음 [CHECK]**. 인스턴스 이름 `instance-20260925-221723` | `GATE_R3CP_PLAN_KO.md` §4; `Y1_RESULT_KO.md` §2.4 |
| Y1 | Python 3.13.5, mpmath 1.3.0(순수 파이썬, 서버에 설치되어 있지 않아 폴더째 복사해 `PYTHONPATH`로 사용) | `Y1_RESULT_KO.md` §2.4 |
| 검사기 R3c v2 | mpmath 필요. 시스템에 없으면 `mpmath/` 폴더를 `check_r3c_v2.py`·`dyad.py`와 같은 폴더에 둔다(`-I`에서는 `PYTHONPATH`가 무시됨) | 공유 노트 r3c-check-v2 |
| 증명 점검(로컬) | Python 3.12.7, sympy 1.14.0, mpmath 1.3.0 | `PROOF_REVIEW_KO.md` §10 |
| 이 문서 작성 노트북 | Windows 11, Python 3.12.7 (해시 계산만) | 이번 세션 |

---

## 4. 실행 순서

**순서: R4 → Y1 → 검사기 → 증명 점검.** R4가 가장 오래 걸린다. Y1은 R4 결과 파일만 있으면 되므로, 저장된 `results_r4/r4_boxes_128w.jsonl`로 바로 할 수도 있다.

### 4.1 R4 (적분)

`rigorous/` 폴더의 다섯 파일과 `results/ref_B_grid.json`을 같은 상대 경로로 두고:

```bash
python3 -I -B rigorous/r4_run.py --prec 128 --pivot width --wall-limit 172800 --outdir OUT_R4
```

- 결과: `OUT_R4/r4_boxes_128w.jsonl`, `r3c_final_128w.json`, `r3c_grid_128w.jsonl`, `r3c_steps_128w.jsonl`, `r3c_jac_128w.jsonl.gz`, `r3c_progress_128w.log`.
- 기대값: `t_final = 2047006249/33554432`, `steps_accepted = 31811`, `stop_reason = C_T_REACHED`, 격자 위반 0. 결정적 실행이므로(같은 입력 → 같은 걸음) `r4_boxes_128w.jsonl`의 SHA-256이 2.2절과 같아야 한다고 **예상**한다. 다른 기계·다른 파이썬 판에서 비트 단위로 같은지는 **확인하지 않았다.** (확인된 결정성: 서버와 로컬 진단 실행의 첫 4,494걸음이 같았다. `GATE_R4_PLAN_KO.md` §2.)
- 걸린 시간(실측): 128비트 37,015초(약 10.3시간). 144비트 55,633초(약 15.5시간). 112비트는 21,613초에 t ≈ 54.498에서 터짐, 96비트는 9,110초에 t ≈ 27.789에서 멈춤. 네 실행은 같은 서버에서 동시에 돌았다(시작 2026-10-11T14:23:54Z).
- 기록된 `argv`에는 `r3c_run.py` 경로가 찍힌다. `r4_run.py`가 `r3c_run`을 import해 그 `main`을 부르기 때문으로 보인다(추론, 코드 머리말 기준).
- 다른 정밀도는 `--prec 144`, `112`, `96`만 바꾼다(`GATE_R4_ADDENDUM_1_KO.md`).

<!-- src: GATE_R4_PLAN_KO.md §4; results_r4/r3c_final_128w.json; r4_out.tgz → r3c_final_{144,112,96}w.json, started_utc.txt -->

### 4.2 Y1 (탈출 판정)

`yoshida_check/` 폴더에서(입력은 `../results_r4/r4_boxes_128w.jsonl`):

```bash
python3 y1_selftest.py > ../y1_selftest_log.txt
```

```bash
python3 y1_exact.py ../results_r4/r4_boxes_128w.jsonl ../y1_results_exact.json
```

```bash
python3 y1_mpiv.py ../results_r4/r4_boxes_128w.jsonl ../y1_results_mpiv.json
```

```bash
python3 y1_compare.py ../y1_results_exact.json ../y1_results_mpiv.json ../y1_compare_summary.json ../y1_tables.md
```

- 기대값: 자기시험 253개 검사 실패 0(종료 코드 0). 두 판정기 모두 PASS 256개, 처음 PASS k = 5973, 불일치 0.
- 걸린 시간: 판정기 각 1초 미만, 자기시험 몇 초.
- 위 명령은 `Y1_RESULT_KO.md`에 적힌 그대로이며 출력 경로가 `yoshida_check/`의 상위 폴더를 가리킨다. 저장소에 있는 결과 파일은 `yoshida_check/` 안에 있으므로, 실제 실행 위치는 서버 작업 폴더 구조(`~/y1_dev/`)에 맞춰져 있었던 것으로 보인다(추론). 재현할 때는 출력 경로를 원하는 곳으로 바꿔도 된다.

<!-- src: yoshida_check/Y1_RESULT_KO.md 235–238행, §1, §2.4 -->

### 4.3 독립 검사기

```bash
python -I rigorous/check_dyad.py rigorous/dyad.py --out R0_CHECK_RESULT_dyad.json
```

```bash
python -I rigorous/check_r1.py rigorous/taylor_step.py --out R1_CHECK_RESULT_taylor_step.json
```

```bash
python -I rigorous/check_r3c_v2.py rigorous/lohner_step_width.py --out R3C_CHECK_V2_lohner_step_width.json --ref-cache rigorous/r3c_v2_selftest_runs/ref_cache.pkl
```

- 기대값: R0 포함 위반 0(S1 477,198, S2 434,956, S3 17,000건). 검사기의 보조 항목 R0-API는 혼합 정밀도 `intersect` 126건 때문에 FAIL로 찍힌다. 이것은 "정밀도 하나" 규칙(`rigorous/R0_DECISION_KO.md`)으로 범위 밖으로 분류했다. R1 위반 0. R3c v2(width) PASS, 위반 0, 거부 0.
- 걸린 시간: R0 52.7초, R1 83.4초(실행 기록). R3c v2 본 실행 시간은 보고서에 없음 **[UNKNOWN]**. v2 자기시험(참조 구현, 캐시 없이) 849초.
- 저장소 경로 기준 명령이다(`rigorous/R0_FINAL_RUN_KO.md`, `R1_FINAL_RUN_KO.md`, 공유 노트 r3c-check-v2에 적힌 형태를 이 폴더 기준으로 옮김).

### 4.4 증명 점검 (선택, 증명을 대신하지 않음)

`proof_review/` 폴더에서:

```bash
python -I -X utf8 symbolic_checks.py
```

```bash
python -I -X utf8 lemma_adversarial.py 20261012
```

- 동적 반증 시도(서버, 2코어, 약 11분): `python3 -I star_dynamic_falsify.py {G1|G2|M1|M2|M3} N SEED OUT`. 사용한 값은 `dynamic_raw/run_all.sh`.
- 기대값: sympy 33/33, 음성 대조 4/4. 적대 탐색 위반 0. 동적 반증 (★) 표본 800개 위반 0.

<!-- src: proof_review/PROOF_REVIEW_KO.md §6, §10 -->

---

## 5. 서버 결과 보관본 `server_archive_20261012/`

| 파일 | 크기(바이트) | SHA-256 (`SHA256SUMS_part1`) | 재계산 |
|---|---|---|---|
| `r4_out.tgz` | 119,743,991 | `9592c1ba3621567425e2ca54f0a541142ce0b4b9eeac4f3660cacc755144dea5` | 일치 |
| `batch1.tgz` | 39,386,825 | `61a8f0d6a53ea0ab0bae86139a59a85c8b945c5d6b2d883c045dff6b12835a16` | 일치 |

- 파일 이름이 `SHA256SUMS_part1`이다. **part2 이후가 있는지, 아직 내려받지 않은 보관본이 있는지는 모른다 [UNKNOWN].**

### 5.1 `r4_out.tgz` — 서버 `~/p3b_r3cp/out_r4/` 전체 (R4 네 실행)

| 정밀도 | 들어 있는 것 | 비고 |
|---|---|---|
| 128 (R4-A) | `r3c_final_128w.json`, `r4_boxes_128w.jsonl`, `r3c_grid_128w.jsonl`, `r3c_steps_128w.jsonl`(약 39 MB), `r3c_jac_128w.jsonl.gz`, `r3c_progress_128w.log` | 주 결과. 로컬 `results_r4/`에는 요약과 상자 파일 두 개만 있다 |
| 144 (R4-B) | 같은 종류 + `r4_boxes_144w.jsonl`(684,674 바이트) | t = 61.000086 도달. 상자는 Y1에서 **판정하지 않음** |
| 112 (R4-C) | 같은 종류, `r4_boxes_112w.jsonl`은 **0 바이트** | t ≈ 54.498에서 터져 저장 구간에 못 닿음 |
| 96 (R4-D) | 같은 종류, `r4_boxes_96w.jsonl`은 **0 바이트** | t ≈ 27.789에서 멈춤 |
| 공통 | `started_utc.txt`(2026-10-11T14:23:54Z), `stdout_*.txt`·`stderr_*.txt`(모두 0 바이트) | |

- 적분기 코드 자체(`~/p3b_r3cp/rigorous/`)는 이 보관본에 **들어 있지 않다.** 실행에 쓴 코드가 로컬 봉인 사본과 같다는 근거는 실행 기록 `r3c_final_128w.json`의 `sha256` 항목(lohner_step, r3c_run, dyad, taylor_step 해시)이 봉인값과 같다는 것이다.

### 5.2 `batch1.tgz` — 서버 작업 폴더 여러 개 (672개 항목)

| 서버 폴더 | 내용(목록으로 확인한 범위) |
|---|---|
| `p3b_an/` | R3c 분석 실행(`r3c_analyze.py`, `r3c_results/`의 128w 기록·분석, 당시 `lohner_step.py`·`taylor_step.py`·`dyad.py` 사본, `ref_B_grid.json`) |
| `p3b_r3cp/` | R3c′·R4 계획·봉인 사본, `out_r4/` 일부(4개 항목) |
| `r4smoke/` | R4 감싸개 스모크 시험 |
| `r3c_check_v2_dev/`, `r3c_v2_final/` | 검사기 v2 개발·자기시험과 실제 구현 최종 실행(`V2_lohner_step_width.json` 등) |
| `y1_dev/`, `y1_verify/` | Y1 판정 세션 실행과 조율 쪽 재실행 |
| `proof_review_dev/` | (★) 동적 반증 시도 결과(`out_G1`~`out_M3.jsonl`) |
| `fastcore_dev/`, `fastcore_check_dev/`, `fc1_final/`, `fc4v/` | C++ 이식판(FC1) 개발·검사·최종 실행. 증명 기록에는 쓰지 않음 |
| `archive_20261012/tmp_ord` | **측정 부록(증명에 불필요).** 조율 세션이 서버 /tmp에서 돌린 측정 실험 기록(26개 항목). `ord20`·`ord30`·`ord40`: 테일러 차수 비교(128비트 width, t=0→1.8): 20차 387초·216걸음, 30차 375초·77걸음, 40차 342초·45걸음. `pp`: PyPy 7.3.19로 300걸음, CPython 기록과 걸음마다 일치, 1.35배. `prof.out`: CPython cProfile 120걸음, 테일러 계수 계산 85%. 출처: 조율 세션 설명(2026-10-12) + `GATE_FC1_PLAN_KO.md` §1(숫자 대조함), 폴더 목록은 tgz에서 확인 |

- 코드 검토 세션의 서버 폴더(`~/integrator_review_dev/`)와 R3c′ 검사기 v1 실행 결과(`~/p3b_r3cp/out/C0a.json`, `C0b.json`)는 이 보관본 목록에서 보이지 않았다 **[확인 필요]**.

---

## 6. 재현에서 확인하지 않은 것

- **R4를 처음부터 다시 돌려 같은 해시가 나오는지 확인한 적이 없다.** 같은 서버에서 한 번 돌렸다. 다른 파이썬 판·다른 기계의 비트 동일성은 모른다.
- Y1은 조율 세션이 서버에 새로 올린 파일로 다시 실행해 같은 결과를 얻었다(`Y1_FINAL_JUDGMENT_KO.md`).
- 144비트 상자로 Y1을 돌리지 않았다.
- `ref_cache.pkl`(검사기 v2의 기준값 캐시)은 이 문서에서 해시를 기록하지 않았다.
- `mpmath` 폴더 해시는 Y1 보고서에 디렉터리 해시(`46baa23e…`)로만 있다. 재현 꾸러미에 mpmath 사본을 넣을지 정해야 한다 **[사용자 결정]**.
- 저장소는 Git에 아직 커밋되지 않은 파일이 많다(`research/`가 미추적). 어느 커밋을 공개 기준으로 할지 **[사용자 결정]**.

---

## 7. 공개 경로 (Zenodo DOI 또는 GitHub)

> 📦 **쉬운 설명**
> - **비유:** 연구 결과물을 "번호표가 붙은 공공 창고"에 맡기는 것이다. 한 번 맡긴 상자는 번호표(DOI)로 언제든 찾을 수 있고, 나중에 고친 상자는 새 번호표를 받는다.
> - **실제 의미:** DOI는 영구 식별 번호다. 논문에서 "재현 꾸러미: DOI xxx"로 인용하면 검토자가 정확히 같은 파일을 받는다. arXiv 추천인이 없어도 쓸 수 있는 경로로 검토 중이다.

**공개 단위 제안 (설계 초안, 확정 아님):**

| 묶음 | 내용 | 크기(대략) |
|---|---|---|
| A. 코드·계획·봉인 | 2.1·2.3·2.4·2.5·2.6절 파일, `GATE_*_SEAL.json`, 이 문서, 논문 초안 | 수 MB (`ref_B_grid.json` 5.6 MB 포함) |
| B. 핵심 결과 | `results_r4/`(상자 302개·요약), Y1 결과 JSON, 검사기 결과 JSON | 1 MB 미만(추정) |
| C. 전체 실행 기록 | `server_archive_20261012/r4_out.tgz`(약 120 MB, 네 실행의 걸음 기록 전부) | 약 120 MB |
| D. 개발 보관본 | `batch1.tgz`(약 39 MB) — FC1 개발물·임시 폴더 포함 | 약 39 MB |

- A·B만으로 Y1 판정과 검사기 실행은 재현된다. R4 재실행에는 A만 있으면 된다(약 10시간). C는 걸음 단위 대조용이다.
- D에는 측정 부록(`archive_20261012/tmp_ord`: 차수 비교·PyPy·프로파일, 증명에 불필요)과 FC1 개발물 등 다른 서버 작업 사본이 섞여 있다. 증명 재현에는 필요 없다. 공개할 경우 "측정 부록"으로 표시한다. 공개 여부 **[사용자 결정]**.
- 서버 경로·인스턴스 이름·사용자 이름이 실행 기록(`r3c_final_*.json`의 `argv`, 보고서)과 tar 머리(소유자 이름)에 들어 있다. **결정(2026-10-12, 사용자): 경로를 가린 사본을 따로 만들어 공개한다.** 방법은 8절.
- 공개 전 할 일:
  1. 위 파일들의 해시를 공개본에서 다시 계산해 이 문서와 대조한다.
  2. 라이선스를 정한다 **[사용자 결정]**.
  3. 서비스별 크기 한도와 DOI 발급 방식을 공식 문서로 확인한다 **[CHECK]**.
  4. 논문 초안의 재현성 절(9절)에 DOI 또는 저장소 주소를 넣는다.
  5. mpmath 사본 포함 여부를 정한다(6절).


---

## 8. 경로를 가린 공개용 사본 (계획 — 아직 만들지 않음)

> 결정: 사용자, 2026-10-12 ("경로를 가린 사본을 따로 만들어 공개 하자"). 사본 생성은 공개 직전에 조율 쪽이 한다.
> 이 절과 스크립트 `paper/redact_verify.py`, 규칙 틀 `paper/redact_rules.example.json`은 2026-10-12에 **작성만** 했다. 스크립트는 실행하지 않았다(구문 검사만 함). 공개 전에 작은 가짜 폴더로 먼저 시험해야 한다.

> 📦 **쉬운 설명**
> - **비유:** 원본 장부는 금고에 그대로 두고, 복사본에서 이름 칸만 먹칠한다. 그리고 "먹칠한 칸에 원래 이름을 다시 써 넣으면 원본과 글자 하나까지 똑같다"를 기계로 확인한다. 그러면 이름 말고는 숫자 하나도 안 바뀌었다는 것이 보장된다.
> - **실제 의미:** 정확한 문자열 치환 → 역치환 → 바이트 비교. 봉인 해시는 원본 기준 그대로이고, 대조표가 원본 해시와 사본 해시를 이어 준다.

### 8.1 원칙
1. **원본은 바꾸지 않는다.** 봉인(`GATE_*_SEAL.json`)과 이 문서 2절의 해시는 원본 기준으로 유지한다.
2. 공개용 사본에서만 다음 문자열을 **정확한 문자열 치환**으로 바꾼다. 숫자·분수·계산값은 바꾸지 않는다.

| 원본에 나오는 것 | 사본 | 규칙 id |
|---|---|---|
| GCP 인스턴스 이름(`instance-…`) | `INSTANCE` | `gcp_instance` |
| GCP 프로젝트 ID(`project-…`, 전체 문자열) | `PROJECT` | `gcp_project` |
| 이메일 주소 | `EMAIL_REMOVED` (제거 표시) | `email` |
| 사용자명: `/home/<사용자명>` → `/home/USER`, `C:\Users\<사용자명>` → `C:\Users\USER`, 단독 등장 → `USER` | `USER` | `username` |

- 경로 세 형태는 사용자명 하나를 `USER`로 바꾸는 규칙 하나로 모두 처리된다(경로의 나머지 글자는 그대로).
- 실제 문자열은 비공개 규칙 파일에만 적는다. `redact_rules.example.json`을 공개 폴더 **밖**으로 복사해 `<<…>>`를 채운다. 채운 파일은 공개하지 않는다. 스크립트는 `<<`가 남아 있으면 실행을 거부한다.
- 자리표시가 원본에 원래 들어 있으면(예: 셸 스크립트의 `USER`) 역치환이 모호해진다. 스크립트는 원본 전체에 없는 첫 후보를 고른다(`USER` → `REDACTED_USER` → `XUSERX`). 실제로 고른 자리표시는 대조표에 기록된다.
- 대소문자가 다른 변형(예: 대문자로 쓴 사용자명)은 정확한 규칙에 걸리지 않으므로 `scan`이 FLAG로 알린다. 사람이 처리한다.

### 8.2 파일 종류별 처리
| 종류 | 처리 |
|---|---|
| 텍스트(.md .py .json .jsonl .log .txt .sh .csv .svg 등) | 바이트 단위 치환 |
| 단일 파일 `.gz`(`r3c_jac_*.jsonl.gz` 등) | 풀어서 내용에 치환, gzip(mtime=0)으로 다시 압축. 검증은 풀린 내용 기준 |
| 그 밖의 이진 파일(.so .pkl .tar .tgz .zip 등) | 바꾸지 않고 그대로 복사. secret이 들어 있으면(압축된 tgz는 풀어서 훑음) FLAG로 보고하고 `redact`를 거부 → 사람이 처리 |
| 파일·폴더 **이름** | 바꾸지 않음. 이름에 secret이 있으면 FLAG |

- `r4_out.tgz`, `batch1.tgz`는 먼저 원본 폴더에 풀고, 사본 폴더를 다시 묶는다. **tar 머리의 소유자 이름**(현재 목록에 사용자명이 소유자로 보임)은 다시 묶을 때 숫자 0 등으로 바꿔야 한다. GNU tar라면 `--owner=0 --group=0 --numeric-owner --sort=name`류 옵션을 쓴다 **[CHECK: 옵션 이름은 공식 문서로 확인]**.
- `batch1.tgz` 안의 중첩 보관본(`y1_out.tar`, `fc1.tar`, `fakes.tar`, `pack*.tgz`, `all_out.tgz` 등)은 `scan`이 FLAG를 낼 가능성이 크다(추정). 풀어서 따로 처리하거나 공개에서 뺀다 **[사용자 결정]**.

### 8.3 대조표 `REDACTION_MAP.json` (공개함, secret은 적지 않음)
- 규칙마다: id, 설명, 실제로 쓴 자리표시, 적용 순서(긴 secret 먼저).
- 파일마다: 상대 경로, 종류, **원본 SHA-256**, **사본 SHA-256**, 규칙별 치환 횟수. `.gz`는 풀린 내용의 원본·사본 SHA-256도.
- 봉인 해시(원본) → 대조표 → 사본 해시로 이어진다. 외부 사람은 사본 해시가 대조표와 같은지 확인할 수 있다. 역치환 바이트 비교는 원본을 가진 저자 쪽이 하고, 그 결과(`verify` 출력)를 함께 공개한다.

### 8.4 절차 (공개 직전, 조율 쪽)

```bash
python -I redact_verify.py scan --rules PRIVATE_RULES.json --src ORIG
```

FLAG·PROBLEM이 0이 될 때까지 사람이 처리한다(이진 파일, 대소문자 변형, 이름 속 secret, 자리표시 충돌).

```bash
python -I redact_verify.py redact --rules PRIVATE_RULES.json --src ORIG --dst COPY --map REDACTION_MAP.json
```

```bash
python -I redact_verify.py verify --rules PRIVATE_RULES.json --src ORIG --dst COPY --map REDACTION_MAP.json
```

`verify`가 확인하는 것(하나라도 어기면 FAIL, 종료 코드 1):
1. 원본·사본의 파일 목록이 대조표와 같다.
2. 원본·사본 해시가 대조표와 같다.
3. 고른 자리표시가 원본에 원래 없다(역치환이 모호하지 않다).
4. **사본에 역치환을 하면 원본과 바이트 단위로 같다**(`.gz`는 풀린 내용). → 바뀐 것은 규칙 문자열뿐이다.
5. 치환 횟수가 대조표와 같고, 사본 속 자리표시 개수가 치환 횟수와 같다.
6. 사본에 secret이 대소문자 무시로도 남아 있지 않다.
7. 이진 파일은 원본과 바이트 단위로 같다.

### 8.5 한계
- 스크립트는 작성만 했고 실행하지 않았다. 실제 원본에서 FLAG가 얼마나 나올지 모른다 **[UNKNOWN]**.
- 바이트 치환은 UTF-8 텍스트를 가정한다. UTF-16 등 다른 인코딩 파일에 secret이 있으면 정확 규칙과 `scan`이 놓칠 수 있다.
- 압축 tgz 안쪽은 훑기만 하고 바꾸지 않는다. zip 안쪽은 훑지 않는다(FLAG만).
- 공개 대조표만으로 외부 사람이 역치환 검사를 직접 할 수는 없다(원본이 비공개이므로). 외부 사람이 확인할 수 있는 것은 사본 해시와 대조표의 일치, 그리고 공개된 `verify` 출력이다.
- 이 문서 자체(REPRODUCE.md)와 논문 초안·보고서에도 사용자명이 든 경로가 있다. 공개할 때는 이 문서들도 같은 방식으로 가린 사본을 쓴다.
