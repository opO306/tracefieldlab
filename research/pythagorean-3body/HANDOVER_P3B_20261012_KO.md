# 인수인계서 — 피타고라스 삼체 컴퓨터 보조 증명 (2026-10-12 밤)

새 세션은 이 문서 → `STATUS_20261012_KO.md` → 필요한 판정 문서 순서로 읽는다. 사용자는 비전공 최종 설계자(한국어만), 필명 **Jeomseon**.

## 1. 한 줄 상태
"내부 검토를 통과한 컴퓨터 보조 증명 초안 (외부 검토 대기)". **Zenodo 게시 직전.** 신규성("세계 최초") 주장 금지, "찾지 못함(부재 증명 아님)"까지만.

## 2. 지금 진행 중인 것 (이 세션이 끊긴 지점)
- 사용자가 Zenodo 초안(로그인 계정 jeomseon@proton.me)에 파일 3개를 올리는 중. C 묶음(119MB) 업로드가 느림.
- 예약 DOI: ~~10.5281/zenodo.23315691~~ → **10.5281/zenodo.23317644** (2026-10-13 변경: 23315691 초안 업로드가 계속 실패해 23317644 초안을 쓰기로 사용자 결정. DRAFT_EN 450행·PDF·A.zip 안 사본 갱신). 메타데이터 칸은 사용자가 채우고 Save draft 함(Preprint, CC BY 4.0(+MIT), 키워드, v1, 요약문).
- 올릴 파일 (폴더 `research/pythagorean-3body/zenodo_upload/`, Git 제외, `SHA256SUMS.txt` 있음):

| 파일 | 바이트 | SHA-256 |
|---|---|---|
| `Jeomseon_2026_Pythagorean_three-body_v1.pdf` | 642,961 | `d7b6b7e528ee5f2a26c95c2f5087b17a1488582ba92aa5f7dc1dbbd7501cf7fa` (DOI 갱신판) |
| `jeomseon_p3b_release_A.zip` (86개) | 3,751,768 | `04aab00a827e87634e8df32c559695c8a90f50fa7d1150ee0e805becb5996f5b` (paper/DRAFT_EN.md·PDF 두 항목만 교체) |
| `jeomseon_p3b_release_C.tar.gz` (Zenodo에는 10MB 조각 12개 `part0{0,1,2}_{0..3}`, md5 `zenodo_upload/C_PARTS_MD5.txt`) | 119,476,199 | `41f338db76c737af7413a6bca81e7c6c82ab42db4938069fc3d3239dcc0cdcb6` |

- 세 파일 모두 개인 문자열(zun24, GCP 인스턴스·프로젝트 ID, lamjwid94@gmail.com, opO306, tracefieldlab) 0건 확인(조율 세션 직접 검사). 가리기 대조표 `paper/release/REDACTION_MAP.json`, 목록 `paper/release/RELEASE_FILE_LIST.md`(묶음 갱신 전 판일 수 있음 — 해시는 위 표가 최신).

## 3. 다음 단계 (순서대로)
1. 사용자가 업로드 완료 → **Preview 캡처 확인**(제목, 저자 Jeomseon, DOI, 라이선스, 파일 3개) → 사용자가 **Publish**. 게시는 되돌릴 수 없음. 조율 쪽은 누르지 않는다.
2. 게시 후: DOI 링크가 열리는지 확인, 내려받은 파일 해시를 위 표와 대조.
3. Git: 커밋 대상 — 커밋 `5f7e563` 이후 바뀐 문서(DRAFT_EN·KO, REPRODUCE, PROOF_STAR_FINAL, LITERATURE_RECHECK, STATUS, build_pdf.py, R4R 판정, release 목록 등). **커밋은 사용자 승인 후**, 브랜치 `p3b-proof-record-20261012`. 큰 파일(*.tgz, zenodo_upload/)은 Git 밖.
4. 푸시는 사용자가 직접: `git push private p3b-proof-record-20261012` (remote `private` = 비공개 저장소 tracefieldlab). **공개 저장소 origin(numerical-audit-lab)에는 절대 푸시하지 않는다**(경로·필명 노출).
5. 전문가 문의: 초안 `paper/EXPERT_EMAIL_DRAFT.md`(영어+한국어, 발신 jeomseon@proton.me). 받는 사람 이름·주소는 공식 페이지에서 확인 후 기입, [LINK]에 DOI. **발송은 사용자 결정.**
6. 남은 약점(논문 8절): 1 ulp 오류는 코드 읽기만, 충돌 정칙화는 인용, 형식 검증 없음, 외부 심사 없음, MathSciNet 미검색. 사용자가 받은 외부 AI 검토(별도 세션 EXTERNAL_REVIEW_KO.md, 프로젝트 밖)가 정칙화를 인용 없이 닫는 스케치를 제안함 — 미반영.

## 4. 증명 구조와 근거 파일
| 고리 | 결과 | 파일 |
|---|---|---|
| R4: t=0~61 엄밀 상자 (128비트, width 피벗) | t=61.0055, 위반 0, 상자 302개 | `GATE_R4_PLAN_KO.md`, `results_r4/` |
| R4R: C++ 경로 전 구간 재실행 | 31,811걸음 비트 동일, 2.2h | `R4R_FINAL_JUDGMENT_KO.md` |
| Y1: 요시다 수정 조건 (★) 판정 | k=5973~6100 PASS, 이전 0, 판정기 2개 | `Y1_FINAL_JUDGMENT_KO.md` |
| (★) 증명 | 독립 검토, 보조정리 A·B 손 증명 | `proof_review/PROOF_STAR_FINAL_EN.md` |
| 적분기 코드 검토 / 검사기 v2 | 결함 0 / 미시험 0·위반 0 | `integrator_review/`, `rigorous/R3C_CHECK_V2_FINAL_KO.md` |
| FC1 C++ 핵심부 | 원본과 비트 동일 23,452건, 4.5배 | `FC1_FINAL_JUDGMENT_KO.md` |
| 문헌 | 선행 엄밀 증명 찾지 못함. Yoshida Corollary 빈틈은 저자가 **1973 Erratum**(PASJ 25, 285)에서 정정. Brutus 136비트·t≈63.4는 **2018** 논문(CNSNS 61) | `LITERATURE_RECHECK_KO.md`(+부록), `LITERATURE_RECHECK_ZBMATH_KO.md` |

## 5. 환경·규칙
- GCP: 인스턴스 **중지됨(TERMINATED, 2026-10-12 03:24 PDT)**, 서버 자료 보존. 크레딧 약 5일 남음. 다시 켤 때 매번 사용자 승인. `CLOUDSDK_CONFIG='C:\Users\zun24\gcloud-task9'`.
- OCI(ARM 24코어): `ssh -i C:\Users\zun24\.ssh\oci_p3b.key ubuntu@163.192.123.153`, 아직 사용 안 함.
- 서버 결과 보관본: `server_archive_20261012/*.tgz` 4개(Git 밖, `MANIFEST_SHA256.txt`). 사용자에게 클라우드 백업 권고함(완료 여부 미확인).
- 노트북 CPU는 사용자가 다른 일에 씀 → 무거운 계산은 서버로.
- 금지: 승인 없는 커밋·푸시, 공개 저장소 푸시, 서버 자료 삭제, 계정 생성·약관 동의·게시 버튼(사용자 몫), 봉인 계획 사후 수정.
- 메모리(`~/.claude/projects/.../memory/`): p3b-proof-status, pen-name-jeomseon, tools-allowed-for-proof, heavy-compute-on-gcp, oci-server 등.

## 6. 그 뒤 후보 (사용자와 논의만 함)
- A: 피타고라스 주변 초기조건 상자 전체 증명. B: 새로 보고된 삼체 주기 궤도 엄밀 검증(위상적 방법 필요). 둘 다 선행 연구 조사 먼저.
- 속도: FC2(나머지 15% C++), 릴레이(구간 병렬), 차수 재조정(B 이후).
- 부수입 목표(크몽)는 별도로 남아 있음.
