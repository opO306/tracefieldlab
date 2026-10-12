# 공개 꾸러미 파일 목록 (Zenodo 업로드용, 2026-10-12 작성)

업로드·게시는 하지 않았다. 사용자가 직접 한다. 원본은 수정하지 않았다.

## 묶음

| 묶음 | 크기(바이트) | SHA-256 |
|---|---|---|
| `jeomseon_p3b_release_A.zip` | 3751819 | `ca3668180ef44bc066d5d9155d1332566266285ebcb1e35ecf306c3fa67822c2` |
| `jeomseon_p3b_release_C.tar.gz` | 119476199 | `41f338db76c737af7413a6bca81e7c6c82ab42db4938069fc3d3239dcc0cdcb6` |

- A 묶음(A+B 범위) 안 경로 = `jeomseon_p3b_release_A/` + 아래 표의 `A/` 뒤 경로(= `research/pythagorean-3body/` 기준 상대 경로).
- C 묶음 안 경로 = `jeomseon_p3b_release_C/` + 아래 표의 `C/` 뒤 경로(= `server_archive_20261012/r4_out.tgz` 안 경로).
- 묶음 안 날짜는 모두 2026-10-12 00:00 UTC로 고정, 소유자 0/0·이름 없음, 파일 권한 644. 이름 정렬. (C는 Python tarfile+gzip mtime=0로 만들었다. GNU tar 옵션과 같은 목적.)
- `.gz` 4개(C의 `r3c_jac_*w.jsonl.gz`)는 redact_verify가 다시 압축했으므로 압축 바이트는 원본과 다르다. 대조표의 해시는 압축 파일 기준이고, 검증(역치환 바이트 비교)은 풀린 내용 기준이다. 이 4개는 치환 0회.
- `orig_sha256`가 봉인(`GATE_*_SEAL.json`, REPRODUCE §2)의 기준이다. `copy_sha256`가 공개 사본 해시다. 둘이 다른 파일만 '치환' 열에 횟수가 있다.

## 파일 (총 119개)

| 경로 | 크기(사본) | 사본 SHA-256 | 원본 SHA-256(다를 때만) | 치환 |
|---|---|---|---|---|
| `A/FC1_FINAL_JUDGMENT_KO.md` | 2221 | `83e3c8c2f951df2f48005b73ecf6d0fb66547eaa444990dcb26689be80347a48` |  | - |
| `A/GATE_FC1_PLAN_KO.md` | 8219 | `37276492d6855c2e9598d3ad69869e4650b8db2802db0b9c224f3c26670552ed` |  | - |
| `A/GATE_FC1_SEAL.json` | 913 | `f7006be638f8b11fc1ffb1044215fbb103ce610dfc8460ac4952115f71ba02f5` |  | - |
| `A/GATE_P3B0_PLAN_KO.md` | 9723 | `416b51196fedfbdfc84e190d56d4ca28a83d291099e430eef1ffcd90fbb84551` |  | - |
| `A/GATE_P3B0_SEAL.json` | 315 | `9cee6fe973a97e416ebcacfd2a50ffd8dbeb08568c85f43eb3b4f5312194e156` |  | - |
| `A/GATE_R0_PLAN_KO.md` | 4789 | `f8537baa4bde7ada8a851e7c4fbc1390e13da528266a34b9e88cee5d3a6d5a97` |  | - |
| `A/GATE_R0_SEAL.json` | 320 | `ff3457d7d2d530f2b2841eeb4256ff35f1ab8cc18dfb29b24cd00a55783af369` |  | - |
| `A/GATE_R1_PLAN_KO.md` | 6153 | `edc5d648531a57d7165870d4ca9b4ce397703a3255ca042cd15f0fc8fc5e611f` |  | - |
| `A/GATE_R1_SEAL.json` | 637 | `20c874169ebb39df12619e3caeb17030d8a181b4e5ed056816d53f03f95ac86e` |  | - |
| `A/GATE_R2_PLAN_KO.md` | 4451 | `79ed3d4ebfdf8f503cea4b215c1fa447de478282fea09f5a71a5ba0ce76b02c5` |  | - |
| `A/GATE_R2_SEAL.json` | 746 | `1b14504315601a7e4e5f81175cf781e3d56d81642f9303c8fda5278edc5afbfd` |  | - |
| `A/GATE_R3CP_ADDENDUM_1_KO.md` | 1638 | `da79ed0f9de6c722fb0f7de4893ed9f5682729406f27f5abcacb606be73da8ab` |  | - |
| `A/GATE_R3CP_ADDENDUM_1_SEAL.json` | 170 | `c1ac13dd530d77c0fd290a7818a06f53efca328d11f984c34e9283fe6835310d` |  | - |
| `A/GATE_R3CP_PLAN_KO.md` | 10391 | `5b83160b45c75c6f1d46ee01b4f2bd4f35af45278ba2a0f458da5fb577b4ba4e` |  | - |
| `A/GATE_R3CP_SEAL.json` | 1126 | `e4a12399c1ded092129f2b6818c44cda2fb5191bac70ad0ebb84576409cb8c1c` |  | - |
| `A/GATE_R3C_PLAN_KO.md` | 9620 | `26a6d98727320d718b646f6693c95aea3c4c950a4080a2a61d5774bcb83d560a` |  | - |
| `A/GATE_R3C_SEAL.json` | 718 | `2315f76ece4df82ad6227175c2a4af7df1bea5db9b2855d918a76c41ddcee6c1` |  | - |
| `A/GATE_R4R_PLAN_KO.md` | 2182 | `d7e14e80b16a856aaaa15df3f5cecfe1e2b499fc8ddf7363dec983aaf4e9cfd0` |  | - |
| `A/GATE_R4R_SEAL.json` | 638 | `076b27125eceb3a3489c720cea1764dabedd3f688ce1d1203a6c88900bb49fdd` |  | - |
| `A/GATE_R4_ADDENDUM_1_KO.md` | 1485 | `75e167d26a15578294a64747231e73095deeda3e7763b3e4bea6057ff2be9a4b` |  | - |
| `A/GATE_R4_ADDENDUM_1_SEAL.json` | 166 | `46229771eb403d832736ded59f35b5c8772841a402d2d289724685eb7e65fcbd` |  | - |
| `A/GATE_R4_PLAN_KO.md` | 6111 | `e59ce0ee1ce1a73eb22c5dc5b62f4878f522b54282d6edb37f059b64ac43acce` |  | - |
| `A/GATE_R4_SEAL.json` | 917 | `bbfd40f4c163ad2684a14422b73c7491a6380c16891526f691978596fa672960` |  | - |
| `A/GATE_Y1_PLAN_KO.md` | 7238 | `416b408ee41c63044b90122b1fee790ed14d826a01435c4870e85a706195cab5` |  | - |
| `A/GATE_Y1_SEAL.json` | 561 | `3a61d102210dd094da85f4b6ce63bea40a6fd7391ce1db40ccc6982c5a4fc4fd` |  | - |
| `A/LITERATURE_RECHECK_KO.md` | 24513 | `1b672d01977bdccafdd23773907bb8a36adec9d89bc27311de0afc777b9159ca` |  | - |
| `A/LITERATURE_RECHECK_ZBMATH_KO.md` | 2400 | `538d859fb06add57bf0b0bbe77756b143c7881c5a9a0610c1e34dbf375e11430` |  | - |
| `A/PRIOR_ART_R3_METHODS_KO.md` | 38804 | `09babb30c1c49cede86fad12994e53b5fb752df931dc7e8f2d96abe6057f6fb4` |  | - |
| `A/PRIOR_ART_RIGOROUS_KO.md` | 26595 | `5e200a3e1e31fde085a7a7642568df8dd7faf8f78c4c012c2fc6f4d11e51c825` |  | - |
| `A/R4R_FINAL_JUDGMENT_KO.md` | 1349 | `a72172232c173ed63ce570c0978763b29bad3c611c9425cb74a5584910a34e08` |  | - |
| `A/RESULT_ADDENDUM_1_KO.md` | 1277 | `15e49b60c66e222447e3aaaef15f12bf26fb1972ba1a8c0bedfeed2e68d66f09` |  | - |
| `A/RESULT_KO.md` | 22943 | `f75d3378abf43d3c853127e65362077600737fe7f99a6550e7aff5ca4668ea4e` |  | - |
| `A/ROADMAP_AMENDMENT_1_KO.md` | 1994 | `44e98d3c4badb65e2b0478284b5a365b96f342bc92d19940958f167faf0290ce` |  | - |
| `A/ROADMAP_AMENDMENT_2_KO.md` | 2737 | `6a83bb20d7cd5fe2f9885fec1b8359dae277b069e89ec177d4e5a42bd547e83b` |  | - |
| `A/ROADMAP_AMENDMENT_3_KO.md` | 2132 | `9e080c99858635c2012c5112a0b9f0bda72fb8dcfeca68fd48eaa1ae80e5bef3` |  | - |
| `A/ROADMAP_AMENDMENT_4_KO.md` | 5188 | `33df32cfdb79fdd3029b97982de669f2e5597e8f22f222d551c782d4db41b09c` |  | - |
| `A/ROADMAP_RIGOROUS_KO.md` | 3122 | `a0e861ce686f0d9a406c226840b4f67a555a75af3a50162440051b33995d69d5` |  | - |
| `A/STATUS_20261012_KO.md` | 4135 | `0d6559382bd9ba93f235730c43f95a072763bb171cd126225b16ceeb34757abd` | `9bfed718a8f8bfc8f855eab9f98f4bdb6d4175f57e199667c8cf9485794fe2cb` | github_org×1 |
| `A/SURVEY_KO.md` | 9864 | `f268d105f49e2de0f4f295a4e552389d4bc880bc56ba5bd3aa58911c7b5dd91e` |  | - |
| `A/Y1_FINAL_JUDGMENT_KO.md` | 1955 | `0e6fc2b1392f9d7678d2e0a517a9efaf18378dab283e74635496dd8e2b9b727f` |  | - |
| `A/YOSHIDA_REDERIVATION_KO.md` | 13594 | `fab037d1c5ad124402fea6c9f6a7f733818de63d4af7b181625cc0f83c7f04a5` |  | - |
| `A/fastcore/FC1_IMPL_RESULT_KO.md` | 20334 | `14fd65c11b07111643b2ac239c002df4bceba182cf9a656b9e8c6c1d4a802188` | `d86f17a8b05cd45ac6c248b87704e8c6cbb7442a003b165f1569c71bcb899fbc` | username×1 |
| `A/fastcore_check/FC1_CHECK_REPORT_KO.md` | 42884 | `392e95b09859508f158f390d271267b48918cf14038d07636e6440c1ab547a66` |  | - |
| `A/integrator_review/INTEGRATOR_REVIEW_KO.md` | 38715 | `282bea5cd652d4d7b469ae39cb2f6c0e5a2bcb65ee9025126ef143a1397e4748` |  | - |
| `A/paper/DRAFT_EN.md` | 67741 | `95a77f0b727a4862dec30e0f03782175a7ba681c4dcf73cb000b7f8332b9183d` |  | - |
| `A/paper/Jeomseon_2026_Pythagorean_three-body_v1.pdf` | 642958 | `1e54766dcc09881e2757c5349fc5cc4cf608bfb09246fe7174911d4bfe913e9e` |  | - |
| `A/paper/REPRODUCE.md` | 31368 | `4acf6e807f4b5758d0de74a006c466328a1124959df5c6b9acd9b9a2adb86ce6` | `6e8d732583ff27fd3ee102caa4b232c5c18d24dd5e0caa5177599c94fccbf4c8` | gcp_instance×1 |
| `A/paper/build_pdf.py` | 3390 | `1f35e06f244deff0083cd2b356d29362637fbb11d43a6e85103dea7bf4d8c808` |  | - |
| `A/proof_review/PROOF_REVIEW_KO.md` | 34582 | `9c8857af850c70efe4a3e88ed7beaa0704b5179d1ab469f2f3953ccbb4229b6a` |  | - |
| `A/proof_review/PROOF_STAR_FINAL_EN.md` | 24437 | `4dfaea735cf3a90fa07d610be0c46c57dbc95150400630feb215acfac411ce4f` |  | - |
| `A/proof_review/PROOF_STAR_FINAL_KO.md` | 30248 | `99b81391d97283b5dc924eceeecb502f176e0389e4e403aae23e51657f05f2bf` |  | - |
| `A/proof_review/lemma_adversarial.py` | 5279 | `be533c373df4764a8168d4feac886ca6bebea4bc3570d3e172c4372661cdb5a8` |  | - |
| `A/proof_review/star_dynamic_falsify.py` | 12599 | `dee5494d70ea533f881f425816cc2ba448e51663dc21ba47c072425fbb9c025d` |  | - |
| `A/proof_review/symbolic_checks.py` | 13290 | `0275b32a2ed24d57baebd3321bcd6bbacd3e4078ce3fa322e06095d73c749938` |  | - |
| `A/results/ref_B_grid.json` | 5566246 | `3bfe3e44024d90bfe495c365bebb6c3c2ea7a92c18869d177bd5e6beaf8d3478` |  | - |
| `A/results_r4/r3c_final_128w.json` | 2593 | `49ad3805eb927262d0a05059d7ce944a86b5cd78e0876d78259ca3122e449ce7` | `7241e8c9788617ed882d0ec6ab23b173e4b2a1e93c2fa18a5e1965883827c1f4` | username×3 |
| `A/results_r4/r4_boxes_128w.jsonl` | 613937 | `f65101c194b6f9d4a67573773a9c517c49ccde8f59924954bd6598f6301d59f5` |  | - |
| `A/rigorous/R0_CHECK_RESULT_dyad.json` | 19128 | `b8064bb50dff859b00ad58a3ab0b629ebc57e9e0cd5d66f11fcf3620f49274a7` |  | - |
| `A/rigorous/R0_FINAL_RUN_KO.md` | 3757 | `847efae35914ac375ee4d616469736c252a7af0340c0f219e6e328d3b109e6a2` |  | - |
| `A/rigorous/R1_CHECK_RESULT_taylor_step.json` | 56636 | `7280e009f2c74b3fbd6fffc84f08ced5db2ed9f9ea0b4441a6bcdab70daf8231` |  | - |
| `A/rigorous/R1_FINAL_RUN_KO.md` | 2048 | `7a0732b21287ecdba5090457fa200aedcbaeecc5d8373511b431569e95a1ddcc` |  | - |
| `A/rigorous/R3C_CHECK_V2_FINAL_KO.md` | 1574 | `6d279082043a79d1d2baa327157b27a7ee8df01ea69c1dca889bd22eb25b4f22` |  | - |
| `A/rigorous/R3C_CHECK_V2_NOTE_KO.md` | 20048 | `f3270d030c41ee5834a23478eb6a0f3093778282f275637a37ce6df5ad08c8c7` | `fb6604a6fa4e14858189dc4fefadcef72295e514cd09c53ef2e47c1d82694597` | username×1 |
| `A/rigorous/R3C_CHECK_V2_lohner_step.json` | 20272 | `c8af9ef3e6ed36c310514c137bd260b396031285945e540ba955f63dde4c7f5b` | `014d7c5c0d5f8212b943f9c1818968edb5e7ad9655d265313f0078b535ea7c5c` | username×8 |
| `A/rigorous/R3C_CHECK_V2_lohner_step_width.json` | 19366 | `22bfb59bedee16eea938ca1b06afcc47ed80e2cc6269ae953a1ec8143f05451b` | `001ec9b4377ffce328f39e55159cfd684dfa2b66997012320508b35e0e38ba38` | username×2 |
| `A/rigorous/check_dyad.py` | 48441 | `7323e236a062420b48108216708c76c40bb4e266d5f46a138a59b4ab01d21f91` |  | - |
| `A/rigorous/check_r1.py` | 66390 | `13fe422df8427f58b34cdcf18fe3b6700796236370ce915610b5bce8d5b0592d` |  | - |
| `A/rigorous/check_r3c.py` | 62198 | `031f11e28f63a2767b01c6e8a532983ee9882e8b82b7e02192eb8bbc7cbe45b3` |  | - |
| `A/rigorous/check_r3c_v2.py` | 64399 | `5f13766f2a461f49c8b36e86456574bd38b5f432cd0c384ea0d2d6f4ab078e0c` |  | - |
| `A/rigorous/dyad.py` | 15484 | `4d49cf8f707061d128fce3ed44b59ccb3a6bda14937e32dcd469cd201fbabaab` |  | - |
| `A/rigorous/lohner_step.py` | 34173 | `936138ea73b878c811df8db7b582ddb6883590be0f47cc7ece76d5ffc2db6b62` |  | - |
| `A/rigorous/lohner_step_width.py` | 640 | `8992d8b6c83dc3aea8c2b86413b06b914d5d8c13e177bc0add58386b349cc693` |  | - |
| `A/rigorous/r3c_run.py` | 21055 | `ce113414f17053d41d03a147a8298b3813e11f07cb8a5b1b1a711306f86dbe5b` |  | - |
| `A/rigorous/r4_run.py` | 2999 | `9b10b434e9b12e040091dcb9923e67c7e400f711259579b5855ae20ef2de93d9` |  | - |
| `A/rigorous/taylor_step.py` | 17740 | `0183f0f8dc76dcf9c1a9933d8f16a636ac011bc15e9a95e0749fb2459be01c31` |  | - |
| `A/yoshida_check/Y1_RESULT_KO.md` | 59709 | `0af6222317541e72a022a249e41f9123e926b7f18df0a7a4a358c4592576b4ca` | `753699ca18dbb2fe93bdea869d1863861c5060f396d4f374ea11d13aa2811e73` | gcp_instance×1, username×1 |
| `A/yoshida_check/y1_compare.py` | 12145 | `06a54b5035046962d56087470138d5e33703557387a372ee19aa000107b44988` |  | - |
| `A/yoshida_check/y1_compare_log.txt` | 2437 | `10cd4e57b502f9386de693249308110d69a4a2d0a090229eebd8ab36d043e454` |  | - |
| `A/yoshida_check/y1_compare_summary.json` | 6526 | `16cd67f24ecf8e838e4b8c8bc7a2e54c8f45cd69e2eb89eefff083396de920a8` |  | - |
| `A/yoshida_check/y1_exact.py` | 12719 | `3ade72c85c0a681e0f28af5b94b3a9c36cb522be9aaf79c5f510eb8fbd262c63` |  | - |
| `A/yoshida_check/y1_mpiv.py` | 11366 | `a546ebb395b234b58ce14f3b6eb6be524726edbb39c5b77dcc5e381fa5008f8b` |  | - |
| `A/yoshida_check/y1_results_exact.json` | 332663 | `a664bc87252cb1da4eccf73910d6488ef1f62c2289c1953a1c07ce86606fb2c7` |  | - |
| `A/yoshida_check/y1_results_mpiv.json` | 332761 | `397a7158b4bc8e6265cecc2ad9dd43c9af78c1d99a891f66c76e0a68614470a4` | `121f9e65205ce39ae3b8061e4a9b2f97945e0f2706abc1172af9bd73e1523ce7` | username×1 |
| `A/yoshida_check/y1_selftest.py` | 20070 | `38ac3919d44ab39cc2f851535b02cc53e682b739da9df7032aec020dafb868d0` |  | - |
| `A/yoshida_check/y1_selftest_log.txt` | 33459 | `0185be81b9bd0cb40bb5515ed772d0f17dfe3b53fbe418dd827a33fd5a1888fa` |  | - |
| `A/yoshida_check/y1_tables.md` | 29289 | `020166dfc37a561e5e142768ff7549c681bdab67774af69fe3cf124dd802c66c` |  | - |
| `C/p3b_r3cp/out_r4/r3c_final_112w.json` | 2590 | `082a74fdd7ffddee50b5c15db94609e6c72e89282155ef8e572176dcb1f7cb99` | `e091e3ac5fbbba12e4fb955b6d82b8734a371a78526674e0049a25e11016ecd7` | username×3 |
| `C/p3b_r3cp/out_r4/r3c_final_128w.json` | 2593 | `49ad3805eb927262d0a05059d7ce944a86b5cd78e0876d78259ca3122e449ce7` | `7241e8c9788617ed882d0ec6ab23b173e4b2a1e93c2fa18a5e1965883827c1f4` | username×3 |
| `C/p3b_r3cp/out_r4/r3c_final_144w.json` | 2581 | `0f5797b8d75a7afa50a7b40fb8bf9a536c4bf833a73989dea0707e38f3802f44` | `1c904b259c24ff4f7626b37cf95a63fcf33240976ab57d39bdad7221a323da66` | username×3 |
| `C/p3b_r3cp/out_r4/r3c_final_96w.json` | 2626 | `210e4ea3b583d184b783b120dbb3fef3e6eb5ba8d4db8f7edefc851b6631e88a` | `83f9c91bd0746d643f728ce1918491e4ef5c218ce57a33394a957f8995f63033` | username×3 |
| `C/p3b_r3cp/out_r4/r3c_grid_112w.jsonl` | 1281885 | `7de5426b5273722deb646030f95e2ce57ad8acf426e4454ebb50c2b61868dcc2` |  | - |
| `C/p3b_r3cp/out_r4/r3c_grid_128w.jsonl` | 1447122 | `fe53ded84282033c000d97e4ed541ab68f8c3f04d995b77e0bc7840da4d62c79` |  | - |
| `C/p3b_r3cp/out_r4/r3c_grid_144w.jsonl` | 1455798 | `83fa6cf13afb02961a8fd6727fb667458a46cbd5e5ef16f53bca2835765926b7` |  | - |
| `C/p3b_r3cp/out_r4/r3c_grid_96w.jsonl` | 644863 | `6db2ef7a15451d3fee7deae387c122e8cd961e6b2f317c5c6744f46072f38725` |  | - |
| `C/p3b_r3cp/out_r4/r3c_jac_112w.jsonl.gz` | 13792874 | `67740275c866646e903ef50d9f50ba2971ce20fbcc566e1f3c983d556e18fae5` |  | - |
| `C/p3b_r3cp/out_r4/r3c_jac_128w.jsonl.gz` | 27198424 | `8dec5c761b2d79024245feabafb40b965054e9e78606a7c7cc27c0eddc0a39d2` |  | - |
| `C/p3b_r3cp/out_r4/r3c_jac_144w.jsonl.gz` | 43887794 | `cd5d59eace984dc1cb3ce5e5ee4f3b6919f4364d6477f9a125ce5eb7b4f8ca13` |  | - |
| `C/p3b_r3cp/out_r4/r3c_jac_96w.jsonl.gz` | 4530669 | `846e969c1486393e760235d3c302f6ba80f2c834af1320be1153503d04637f5d` |  | - |
| `C/p3b_r3cp/out_r4/r3c_progress_112w.log` | 108013 | `2cf9565a9d2d3faad4aa4cbc61ffaccfdebd7eea2f6fbecf53d918762a2452ac` | `fb13eb18adcbc426b4234120fff42836573163625f43e6bf641e3b011b62bfba` | username×1 |
| `C/p3b_r3cp/out_r4/r3c_progress_128w.log` | 196194 | `17b28543778458a8f0fb8327ed2922ac02d747bca115e2293207d6f2da743be2` | `9d55265b14f35bdd458622fccf78d060f2cbdffd030e288ec885e3b923c36e58` | username×1 |
| `C/p3b_r3cp/out_r4/r3c_progress_144w.log` | 304333 | `08bcb3f072666cd53eea8877b9cc6cc18e36eddd5787fbb19910c1062db7d1c2` | `698972c3fae55e2e71dbd5f465751a465415e4b19c3f9c04f371a72d2092b542` | username×1 |
| `C/p3b_r3cp/out_r4/r3c_progress_96w.log` | 42141 | `73a25a0db77b1aca7af879f25cc060bf69e6c5195d9df7adf1947d12870dbda7` | `990968b809aebfd75a1d8818d044de7f8738991944da0ce4c17cb2ac112564d5` | username×1 |
| `C/p3b_r3cp/out_r4/r3c_steps_112w.jsonl` | 19560051 | `5856e38ff4380ac774eaa51801cb761da0c377182d31c7300ef7dce000f7d46a` |  | - |
| `C/p3b_r3cp/out_r4/r3c_steps_128w.jsonl` | 39451752 | `6a09a8d511f8a8ed2251cb704d49672482b5a732ee6fadd75b7f0a9bf1476749` |  | - |
| `C/p3b_r3cp/out_r4/r3c_steps_144w.jsonl` | 64873813 | `5bf20e5cd943b7f1fcfafa1c09dbc94d4a2d99fcaac06d57dda0ee71108bcb28` |  | - |
| `C/p3b_r3cp/out_r4/r3c_steps_96w.jsonl` | 6309463 | `a5d959c5dcc7e92ed87e7d3402c2a21b5e3baad4a2f07b699e8235672cdabfb1` |  | - |
| `C/p3b_r3cp/out_r4/r4_boxes_112w.jsonl` | 0 | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |  | - |
| `C/p3b_r3cp/out_r4/r4_boxes_128w.jsonl` | 613937 | `f65101c194b6f9d4a67573773a9c517c49ccde8f59924954bd6598f6301d59f5` |  | - |
| `C/p3b_r3cp/out_r4/r4_boxes_144w.jsonl` | 684674 | `91986c710b0f16fdaa0cbcdbd9f5095ff0c83d16fd1aed96462e7930598329eb` |  | - |
| `C/p3b_r3cp/out_r4/r4_boxes_96w.jsonl` | 0 | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |  | - |
| `C/p3b_r3cp/out_r4/started_utc.txt` | 21 | `b9faea124897c4056bfc314e68cccb9e565a1e8385a551974d72df754a846b99` |  | - |
| `C/p3b_r3cp/out_r4/stderr_112.txt` | 0 | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |  | - |
| `C/p3b_r3cp/out_r4/stderr_128.txt` | 0 | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |  | - |
| `C/p3b_r3cp/out_r4/stderr_144.txt` | 0 | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |  | - |
| `C/p3b_r3cp/out_r4/stderr_96.txt` | 0 | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |  | - |
| `C/p3b_r3cp/out_r4/stdout_112.txt` | 0 | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |  | - |
| `C/p3b_r3cp/out_r4/stdout_128.txt` | 0 | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |  | - |
| `C/p3b_r3cp/out_r4/stdout_144.txt` | 0 | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |  | - |
| `C/p3b_r3cp/out_r4/stdout_96.txt` | 0 | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |  | - |

## 제외

| 대상 | 이유 |
|---|---|
| `server_archive_20261012/batch1.tgz`(D), `batch3.tgz`, `r3cp_out.tgz` | 사용자 결정(D 제외)·작업 범위 밖 |
| `paper/EXPERT_EMAIL_DRAFT.md` | 지시로 제외(개인 연락 초안) |
| `rigorous/r3c_results/` 전체(깨진 `r3c_jac_*.gz` 포함) | 지시로 제외. 포함 목록에 없음 |
| `__pycache__`/`.pyc` | 지시로 제외 |
| 비공개 규칙 파일(secret 원문) | 비공개. 스크래치에만 있음. 묶음·release 폴더에 넣지 않음 |
| `server_archive_20261012/r4_out.tgz` 원본 자체 | 풀어서 가린 사본으로 다시 묶음(C). tar 머리에 소유자 사용자명이 있었음 |

## 추가 (보류였다가 사용자 승인으로 A에 넣음, 2026-10-12)

- Y1_RESULT_KO.md·Y1 로그·표, R3C_CHECK_V2_NOTE_KO.md, R0/R1_FINAL_RUN_KO.md, PROOF_REVIEW_KO.md, 나머지 GATE_*_PLAN/ADDENDUM 문서, ROADMAP_*, RESULT_*, LITERATURE_*, PRIOR_ART_*, SURVEY_KO.md, INTEGRATOR_REVIEW_KO.md, FC1_CHECK_REPORT_KO.md, FC1_IMPL_RESULT_KO.md (31개 추가). 새 문서에서 이메일·전화·실명·주소 형태를 정규식으로 검색: 해당 없음(사람이 전부 읽은 것은 아님).


## 최종판 변경 (2026-10-12, DOI 10.5281/zenodo.23315691 예약 뒤)

- `paper/DRAFT_KO.md` 제외(사용자 확인용, 비공개). `paper/DRAFT_EN_preview.pdf` 제외. `paper/EXPERT_EMAIL_DRAFT.md` 계속 제외.
- 추가: `paper/Jeomseon_2026_Pythagorean_three-body_v1.pdf`(이진, 치환 없음), `paper/build_pdf.py`(PDF 생성 도구, http(s) 외 링크 변환 안 하도록 수정된 판).
- PDF 잔존 검사: 원 바이트 + zlib로 풀리는 내부 스트림 98개에서 secret 0, `file:` 0, `Users` 0. 링크는 dlmf.nist.gov, DOI, ipsj, `mailto:jeomseon@proton.me`(필명 공개 연락처로 보고 가리지 않음 - 조율 세션 확인 대상). zlib로 안 풀리는 스트림(글꼴 등 다른 압축)은 검사하지 않음.
- PDF 첫판(645,165바이트)에는 수식이 링크로 잘못 바뀌어 `file:///C:/Users/<사용자명>/...` 주소 3개가 있었다. 고친 판(642,958바이트)만 넣었다.
## 의존성

- mpmath 1.3.0 필요(사본은 넣지 않음). sympy는 proof_review 점검 스크립트에만 필요.

## 주의

- `paper/REPRODUCE.md`, `paper/DRAFT_EN.md`, `paper/DRAFT_KO.md`는 작성 시점(2026-10-12) Git 작업 트리의 **커밋 안 된 수정본**을 담았다(HEAD `5f7e563`와 다름). 이후 문서를 고치면 묶음을 다시 만들어야 한다.
- 가린 문자열 규칙 6개(규칙 id·자리표시는 `REDACTION_MAP.json`): username, gcp_instance, gcp_project, email, github_account, github_org. 실제로 치환이 일어난 것은 username·gcp_instance·github_org뿐이다.
- 원본 문서에 원래 있던 일반 예시 `/home/USER`(3회), `C:\Users\USER`(1회)는 개인 이름이 아니므로 그대로다.

## 검사 결과

- `redact_verify.py verify`: `VERIFY PASS: 119 files, 0 failures`.
- 별도 잔존 검사(대소문자 무시, 사본 119개 + `.gz` 풀린 내용 + 두 묶음 바이트·zip 항목 내용·이름 + 대조표): secret 0건. `/home/`·`C:\Users\` 뒤 이름은 자리표시 또는 원래 예시 `USER`뿐.
- 확인 안 함: 묶음을 다른 컴퓨터에서 풀어 해시를 다시 계산하는 것(같은 노트북에서 만든 직후 검사만 함).

## DOI 변경 (2026-10-13)
- 게시 초안을 23317644로 바꿈(사용자 결정, 23315691 초안 업로드 실패). DRAFT_EN.md 450행 DOI 교체, PDF 재생성(본문 텍스트는 DOI 외 동일 — pypdf 추출 비교), A.zip은 paper/DRAFT_EN.md·PDF 두 항목만 교체(순서·압축 방식 유지, 86개).
- 새 해시: PDF `d7b6b7e528ee5f2a26c95c2f5087b17a1488582ba92aa5f7dc1dbbd7501cf7fa` (642,961 B), A `04aab00a827e87634e8df32c559695c8a90f50fa7d1150ee0e805becb5996f5b` (3,751,768 B). C 변경 없음. 위 표의 A·PDF 해시는 이전 판.
- 개인 문자열 재검사(A 86개 항목, PDF): 0건.
