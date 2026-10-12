# 전문가 문의 이메일 초안 (보내지 않음 — 발송은 사용자 결정)

작성: 조율 세션, 2026-10-12. 상태: **초안**. 받는 사람·주소·첨부 링크는 비어 있다.

## 0. 보내기 전에 정할 것 (사용자)
| 항목 | 상태 |
|---|---|
| 받는 사람 | 아래 2절 후보 중 선택. **이름·소속·주소는 반드시 공식 학교/연구소 페이지에서 직접 확인** (아래 후보는 기억에 의존한 목록이며 확인하지 않았다) |
| 보내는 주소 | jeomseon@proton.me (2026-10-12 사용자 생성) |
| 첨부·링크 | 비공개 저장소는 남이 못 본다. 보낼 때는 (a) 논문 PDF 첨부, (b) 경로를 가린 공개 사본을 Zenodo에 올린 링크, (c) 비공개 저장소에 상대를 초대 중 하나가 필요 |
| 시점 | [CHECK] 해소와 R4R 재실행 대조가 끝난 뒤를 권함 |

## 1. 원칙
- 짧게(읽는 데 2분). 질문은 두 개만.
- 신규성 주장 금지("first" 금지). "찾지 못했다, 알려 주시면 감사"로.
- AI 사용을 숨기지 않는다.
- 한 번에 여러 명에게 같은 메일을 묶어 보내지 않는다(각자 따로, 개인화).

## 2. 받는 사람 후보 (기억에 의존, 발송 전 확인 필요)
| 후보 | 왜 | 확인할 것 |
|---|---|---|
| CAPD 그룹(야기엘론스키 대학, 폴란드) 연구자들 — 예: Daniel Wilczak, Piotr Zgliczyński, Tomasz Kapela | 구간 산술·Lohner 방법 기반 엄밀 적분의 대표 그룹. 선행 엄밀 증명이 있다면 가장 잘 알 사람들 | 현재 소속·연락처, 이 주제 관련 논문 |
| Carles Simó (바르셀로나 대학) | 삼체 주기 궤도의 컴퓨터 보조 증명 공동 연구로 알려짐 | 현재 활동 여부 |
| Hiroshi Hirayama | 2019년 피타고라스 문제 4배 정밀도 테일러 계산 논문 저자(zbMATH 확인) | 연락처, 논문 본문 |

## 3. 영어 본문 (초안)

**Subject:** Question about prior rigorous proofs of the Pythagorean three-body outcome

Dear Dr. [NAME],

I am an independent researcher writing under the name Jeomseon. I am contacting you because of your work on [THEIR RELEVANT WORK — fill in after checking].

I have written a draft computer-assisted proof of the final evolution of the Pythagorean (Burrau) three-body problem. The claim is that mass 3 escapes for t ≥ 59.73 and masses 4 and 5 stay bound. The approach has two parts:

1. A rigorous enclosure for t ∈ [0, 61]. It uses interval arithmetic on dyadic numbers with outward rounding, Taylor steps with remainder bounds, and a Lohner–QR enclosure with a width-weighted pivot, at 128 bits.
2. A check that the enclosures for t ∈ [59.73, 61] satisfy an escape criterion. The criterion is Yoshida's (1972), with the r(0) ≤ d condition from his 1974 correction. Two independently written interval checkers did this check.

I could not find a published rigorous proof of this outcome. I searched ADS, zbMATH Open, scite and general web search; I did not search MathSciNet.

I would be grateful for answers to two questions:
1. Do you know of an existing rigorous or computer-assisted proof of this result?
2. If you have a few minutes, do you see an obvious flaw in this approach?

For transparency: the code, independent checkers and reviews were produced with an AI assistant (Anthropic's Claude), working in separate sessions under sealed test plans. I made the design decisions and reviewed the results. Draft, code and data: [LINK].

Thank you for your time.

Best regards,
Jeomseon
jeomseon@proton.me

## 4. 한국어 대조 (사용자 확인용, 보내는 글 아님)

**제목:** 피타고라스 삼체 문제 결말의 기존 엄밀 증명에 대한 질문

[이름] 박사님께,

저는 Jeomseon이라는 이름으로 활동하는 독립 연구자입니다. 박사님의 [관련 연구 — 확인 후 기입]을 보고 연락드립니다.

피타고라스(Burrau) 삼체 문제의 결말에 대한 컴퓨터 보조 증명 초안을 만들었습니다. 주장은 "질량 3은 t ≥ 59.73부터 탈출하고, 질량 4·5는 묶여 있다"입니다. 방법은 두 부분입니다.

1. t ∈ [0, 61] 구간의 엄밀한 포함입니다. 이진 분수 구간 산술(바깥 반올림), 나머지 상한을 포함한 테일러 단계, 폭 가중 피벗을 쓴 Lohner–QR 포함을 128비트로 사용했습니다.
2. t ∈ [59.73, 61]의 상자들이 탈출 조건을 만족하는지 확인했습니다. 조건은 요시다(1972)의 기준에 그의 1974년 정정(r(0) ≤ d)을 더한 것입니다. 이 확인은 서로 독립적으로 작성한 구간 판정기 두 개로 했습니다.

이 결과의 엄밀 증명이 출판된 것은 찾지 못했습니다. ADS, zbMATH Open, scite, 일반 웹 검색으로 찾았고, MathSciNet은 보지 못했습니다.

두 가지를 여쭙고 싶습니다.
1. 이 결과의 기존 엄밀 증명이나 컴퓨터 보조 증명을 알고 계신가요?
2. 시간이 되신다면, 이 방법에서 명백한 결함이 보이시나요?

투명성을 위해 말씀드립니다. 코드, 독립 검사기, 검토는 AI 도구(Anthropic Claude)를 사용해 봉인된 시험 계획 아래에서 서로 분리된 세션으로 만들었습니다. 설계 결정과 결과 검토는 제가 했습니다. 초안, 코드, 데이터: [링크]

시간 내 주셔서 감사합니다.

Jeomseon 드림
jeomseon@proton.me

## 5. 메모
- 답장이 영어로 오면 번역과 답장 초안은 조율 세션이 돕는다.
- 답이 없을 수도 있다(흔한 일). 2~3주 뒤 한 번만 짧게 다시 묻는 정도가 예의 범위다(일반적 관행, 추정).
