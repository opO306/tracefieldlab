# R3c 검사기 v2 × 실제 구현 최종 실행 (조율 세션, 2026-10-12)

- 검사기 `check_r3c_v2.py` SHA-256 `5f13766f…` (v1 `031f11e2…`는 수정 없음). 작성: 검사기 v2 세션(`R3C_CHECK_V2_NOTE_KO.md`).
- 실행: GCP `~/r3c_v2_final/`. 서버 원본 dyad.py·taylor_step.py·lohner_step.py·lohner_step_width.py가 봉인 사본(`~/p3b_r3cp/rigorous/`)과 바이트 단위로 같음을 `cmp`로 확인한 뒤 실행. 기준값 캐시 `r3c_v2_selftest_runs/ref_cache.pkl` 사용.

| 대상 | 봉인 기준 판정 | 위반 | 비고 |
|---|---|---|---|
| `lohner_step_width.py` (pivot=width, **R4가 쓴 규칙**) | **PASS** (보충 포함 PASS) | 0 | R3-S1 160, S2 2,356, C 4,098, X4 732 등 전부 PASS. 거부(REFUSED) 0 — v1의 미시험 10건 해소 |
| `lohner_step.py` (pivot=colnorm, 기본값) | UNRESOLVED | 0 | S1에서 REFUSED 2, NOT_RUN 50. 원인: `LohnerStepError: singular … collision configuration in box` — colnorm 규칙에서 상자가 커져 원점(충돌 배치)을 포함하자 구현이 정상적으로 거부. 엄밀성 위반 아님. colnorm은 R4·증명에 쓰지 않았다 |

결과 파일: `R3C_CHECK_V2_lohner_step_width.json` (`001ec9b4…`), `R3C_CHECK_V2_lohner_step.json` (`014d7c5c…`).

결론: 증명에 쓴 적분 규칙(width)에 대해 독립 검사기의 봉인 기준이 미시험 없이 모두 통과했다. 시험 범위 안의 결론이며 형식 증명은 아니다. 검사기가 못 잡는 반올림 한 단위 수준 오류는 코드 검토(`integrator_review/INTEGRATOR_REVIEW_KO.md`)에 기댄다.
