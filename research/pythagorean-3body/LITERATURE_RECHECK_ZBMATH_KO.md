# 문헌 재확인 보충 — zbMATH Open 검색 (조율 세션, 2026-10-12)

`LITERATURE_RECHECK_KO.md`(ADS·scite·Consensus·웹)에서 빠진 **수학 문헌 데이터베이스**를 보충한다.

## 방법
- zbMATH Open 웹 검색 페이지는 자동 접근을 거부했다(HTTP 403). 대신 공개 API `https://api.zbmath.org/v1/document/_search`로 검색했다.
- MathSciNet은 유료 구독이 필요해 **검색하지 못했다(UNKNOWN)**.

## 검색어와 결과
| 검색어 | 건수 | 관련 문헌 |
|---|---|---|
| Pythagorean three-body | 6 | Aarseth–Anosova–Orlov–Szebehely 1994 (두 편: "Global chaoticity in the Pythagorean three-body problem", "Close triple approaches and escape…"), Janssens 2011 ("On Free Fall in the Three Body Problem"), Rose–Dullin 2013 (영 각운동량 정칙화 심플렉틱 적분기), **Hirayama 2019** ("Numerical calculation by quadruple precision higher order Taylor series method of the Pythagorean problem…") |
| Pythagorean three body problem | 8 | 위와 같은 문헌 + 무관 2편 |
| Pythagorean problem three bodies | 5 | 위 1994, 2019 + 무관 |
| Burrau | 44 | Szebehely 1967 "Burrau's problem of three bodies", Burrau 1906·1913 계열. 나머지는 동명이인·무관 |
| three-body interval arithmetic | 6 | Lohner 1988(학위논문, 우리가 쓴 기울인 상자 방법의 원전), Adams 1990, Shibayama–Yagasaki 2009(이등변 삼체 이종 연결, 컴퓨터 보조), Baldomá–Capiński 외 2025. **피타고라스 결말을 다룬 것 없음** |
| validated integration three-body | 7 | 무관 |
| three-body computer-assisted proof escape / rigorous numerics three-body hyperbolic-elliptic | API 오류(결과 없음) | — |

## 판단
- zbMATH Open에서 피타고라스(Burrau) 삼체 문제의 결말을 **엄밀하게(구간 산술·컴퓨터 보조 증명)** 증명했다는 제목의 문헌은 **찾지 못했다**. 부재 증명은 아니다.
- 가장 가까운 최근 문헌은 Hirayama 2019(4배 정밀도 고차 테일러 **수치** 계산)다. 제목상 수치 계산이며, 본문·초록은 확인하지 않았다(UNKNOWN). 논문 초안 Introduction에서 선행 고정밀 수치 연구로 인용 후보([CHECK]).
- Lohner 1988 학위논문은 방법론 원전으로 인용 후보([CHECK]: 피타고라스 예제를 포함하는지 본문 미확인).
- MathSciNet 미검색은 한계로 남긴다.
