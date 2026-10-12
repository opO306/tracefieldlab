# Reviewer quick start (about 10 minutes)

Paper: *A computer-assisted proof of the final evolution of the Pythagorean three-body problem*, Jeomseon, 2026, https://doi.org/10.5281/zenodo.23317644 (preprint; not externally reviewed; not formally verified).

**Claim.** For the Pythagorean (Burrau) problem (masses 3, 4, 5 at rest at the vertices of a 3-4-5 triangle, G = 1), mass 3 escapes for t ≥ 59.73 while masses 4 and 5 remain bound.

**How the claim is split (paper §2, §8):**
1. **Enclosure (R4).** A 128-bit interval Lohner–Taylor integrator encloses the solution on t ∈ [0, 61.0055]. The stored output is `results_r4/r4_boxes_128w.jsonl`: 302 boxes, exact fractions.
2. **Escape criterion (★).** A corrected form of Yoshida's 1972 condition (Yoshida's own 1973 erratum) is proved by hand in `proof_review/PROOF_STAR_FINAL_EN.md`.
3. **Check (Y1).** Two independent judges evaluate (★) on the stored boxes. One uses exact rationals, the other interval arithmetic in mpmath.

## Step 1: run the quick verifier (about 2–3 minutes)

From the root of bundle A:

```
pip install mpmath
python -I -X utf8 verify_quick.py
```

| Step | What it re-checks | Expected |
|---|---|---|
| 1 | SHA-256 of every file listed in the sealed plans (`GATE_*_SEAL.json`) that is shipped in A | All match. One file (`results_r4/r3c_final_128w.json`) differs only by a listed privacy redaction: a server username was replaced in 3 path strings. `REDACTION_MAP.json` gives its original and published hashes. |
| 2 | Self-test of the Y1 judges | 253 checks, 0 failures |
| 3 | Both judges on the 302 stored boxes | 256/256 PASS for both judges. First PASS at k = 5973 (t = 59.73). No PASS at k ≤ 5972. 0 disagreements. Input hash matches the seal. |
| 4 | R0 interval-arithmetic checker (`dyad.py`) | Inclusion checks R0-S1..S4 have 0 violations. The auxiliary item R0-API reports 126 mixed-precision `intersect` cases, so the checker's overall word is "FAIL". These cases are outside the one-precision contract (see `paper/REPRODUCE.md` §4.3) and are the same as the published `rigorous/R0_CHECK_RESULT_dyad.json`. |
| 5 | R1 Taylor-step checker (`taylor_step.py`) | PASS |

Use `--no-checkers` to run steps 1–3 only (a few seconds).

## Step 2: read one thing

`proof_review/PROOF_STAR_FINAL_EN.md` contains the hand proof of (★), Lemmas A and B. This is the part a computer cannot check for you here.

## What this quick path does **not** cover

- **Re-running the R4 integration.** It takes about 10 hours of pure Python, or 2.2 hours with the C++ core. See `paper/REPRODUCE.md` §4.1. Steps 1, 4 and 5 check that the integrator code is the sealed version and that its arithmetic and Taylor layers pass independent tests. They do not re-derive the boxes.
- **The R3c v2 Lohner-step checker**, which is long (`paper/REPRODUCE.md` §4.3).
- **Collision regularization.** The paper cites it rather than implementing it (§8).
- **The one-ulp rounding argument**, which was checked by code reading only (§8).
- **Formal verification.** None has been done.

## The most useful question for a reviewer

If you have limited time: **is the corrected escape criterion (★) and its hand proof correct, and is (★) a sufficient condition for the stated final evolution?** Every computational step above is mechanical. (★) is where human judgment matters most.

AI assistance (Anthropic Claude) was used for code, checking and writing. Design decisions and final judgments are the author's.
