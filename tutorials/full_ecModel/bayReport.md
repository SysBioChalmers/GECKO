# Bayesian Kcat Tuning Diagnostic Report

**Generated:** 2026-06-09 12:05:16

---

## Instructions for Claude

Please analyze this diagnostic report from a Bayesian parameter tuning run and provide specific recommendations for improving the hyperparameter settings.

Focus on:
1. Whether convergence was healthy or problematic
2. Whether source-specific regularization is working correctly
3. Specific hyperparameter adjustments with rationale
4. Any red flags or concerning patterns

---

## Current Hyperparameter Settings

### Initial Uncertainty
```matlab
sigma0logDefault = 0.30  % Default log-space std dev
kcatSources = {'OpenKineticsPredictor', 'brenda', 'custom'}
sigma0logSource = [0.25; 0.20; 0.10]  % Per-source initial uncertainty
```

### Regularization (Source-Specific Constraints)
```matlab
shrinkThrDefault = 3.0  % Unlabelled: σ deviations for full update
shrinkThrSource = [5.0; 10.0; 12.0]  % Per-source shrinkage

varianceCapDefault = 1.50  % Unlabelled: max variance growth
varianceCapSource = [1.50; 1.10; 1.05]  % Per-source caps

forcePriorThrDefault = 3.5  % Unlabelled: snap to prior if dev < this
forcePriorThrSource = [5.5; 11.0; 13.0]  % Per-source force thresholds

sparsityThreshold = 0.50  % Clean up changes < this × sigma0
```

### Sampling & Acceptance
```matlab
scheduleGenerations = [1   2   9  15]
scheduleSamples = [1000   800   600   400]
targetAccept = 10%  % Keep best X% of samples
minKeep = 0.30  % Min fraction to keep
maxKeep = 0.60  % Max fraction to keep
```

### Stopping Criteria
```matlab
rmseThreshold = 0.20  % Target RMSE
maxGenerations = 50  % Hard limit
maxRMSEplateau = 10  % Stop after this many generations without improvement
```

---

## Run Summary

- **Total generations:** 42
- **Converged:** No
- **Final RMSE:** 0.8909
- **Initial RMSE:** 8.6001
- **RMSE reduction:** 89.6%
- **Total parameters:** 4834

### Convergence Pattern

- **Early phase (gen 1-10):** 84.1% RMSE reduction
- **Mid phase (gen 11-20):** 16.0% RMSE reduction
- **Late phase (gen 21+):** 3.7% RMSE reduction
- **Late-stage directedness:** 0.82 (1=converging, 0.5=random, <0.3=oscillating)

---

## Per-Source Analysis

### OpenKineticsPredictor

**Summary:**
- Total parameters: 1175
- Active (shrinkWeight > 0.3): 246 (20.9%)
- Near prior (|σ_post - σ_prior| < 0.1): 929 (79.1%)
- Mean deviation from prior: 1.57 σ
- Variance ratio (σ_post / σ_prior): 1.10

**Change magnitude (|log(final/prior)|):**
- Q25: 0.000
- Median: 0.000
- Q75: 0.000

### brenda

**Summary:**
- Total parameters: 3198
- Active (shrinkWeight > 0.3): 85 (2.7%)
- Near prior (|σ_post - σ_prior| < 0.1): 3198 (100.0%)
- Mean deviation from prior: 0.44 σ
- Variance ratio (σ_post / σ_prior): 1.00

**Change magnitude (|log(final/prior)|):**
- Q25: 0.000
- Median: 0.000
- Q75: 0.000

### custom

**Summary:**
- Total parameters: 207
- Active (shrinkWeight > 0.3): 3 (1.4%)
- Near prior (|σ_post - σ_prior| < 0.1): 207 (100.0%)
- Mean deviation from prior: 0.23 σ
- Variance ratio (σ_post / σ_prior): 1.00

**Change magnitude (|log(final/prior)|):**
- Q25: 0.000
- Median: 0.000
- Q75: 0.000

⚠️ **WARNING:** Variance ratio 1.00 is at cap 1.05 (may be constrained)

### unlabelled

**Summary:**
- Total parameters: 254
- Active (shrinkWeight > 0.3): 115 (45.3%)
- Near prior (|σ_post - σ_prior| < 0.1): 139 (54.7%)
- Mean deviation from prior: 2.68 σ
- Variance ratio (σ_post / σ_prior): 1.23

**Change magnitude (|log(final/prior)|):**
- Q25: 0.000
- Median: 0.000
- Q75: 1.641

---

## Key Diagnostic Traces

### RMSE Evolution (every 5 generations)
```
Gen | RMSE   | Improvement
----|--------|------------
  0 | 8.6001 | -
  5 | 2.5749 | 16.5%
 10 | 1.1191 | 12.4%
 15 | 0.9434 | 0.0%
 20 | 0.9402 | 0.0%
 25 | 0.9250 | 0.0%
 30 | 0.9107 | 0.0%
 35 | 0.8909 | 0.0%
 40 | 0.8909 | 0.0%
 42 | 0.8909 | 0.0%
```

### Acceptance Rates (every 5 generations)
```
Gen | Overall | New Proposals | Boundary Violations
----|---------|---------------|--------------------
  5 |  30.0% |  13.5%       |   2.7%
 10 |  29.9% |  11.3%       |   2.4%
 15 |  30.0% |   5.0%       |   3.0%
 20 |  29.9% |   1.5%       |   3.6%
 25 |  29.9% |   2.0%       |   4.0%
 30 |  29.9% |   2.5%       |   4.3%
 35 |  29.9% |   0.8%       |   4.6%
 40 |  29.9% |   0.2%       |   5.0%
```

### Directedness & Diversity (every 5 generations)
```
Gen | Directedness | Diversity | Sparsity %
----|--------------|-----------|----------
  5 | 1.00        | 3.153     |  99.2%
 10 | 1.00        | 5.740     |  99.2%
 15 | 1.00        | 7.439     |  99.1%
 20 | 0.93        | 8.267     |  98.3%
 25 | 0.86        | 8.623     |  97.2%
 30 | 0.85        | 8.911     |  94.9%
 35 | 0.85        | 9.009     |  93.1%
 40 | 0.80        | 9.208     |  92.2%
```

---

## Potential Issues Detected

### ⚠️ WARNING: Very Low Proposal Acceptance
- Late-stage new proposal acceptance: 1.1% (should be 5-15%)
- **Likely cause:** Proposals too aggressive
- **Recommendation:** Proposal width adaptation may need adjustment

---

## Detailed Diagnostic Data

<details>
<summary>Click to expand full diagnostic traces</summary>

### Complete RMSE Trace
```
0,8.600052
1,7.020520
2,3.842719
3,3.567754
4,3.084181
5,2.574882
6,2.140940
7,1.772273
8,1.528175
9,1.277315
10,1.119115
11,1.119115
12,1.119115
13,1.025652
14,0.943351
15,0.943351
16,0.943351
17,0.943351
18,0.943351
19,0.940233
20,0.940233
21,0.924977
22,0.924977
23,0.924977
24,0.924977
25,0.924977
26,0.910726
27,0.910726
28,0.910726
29,0.910726
30,0.910726
31,0.905443
32,0.905443
33,0.890940
34,0.890940
35,0.890940
36,0.890940
37,0.890940
38,0.890940
39,0.890940
40,0.890940
41,0.890940
42,0.890940
```

### Acceptance Rate Trace
```
Gen,Overall,NewProposals,BoundaryViolations
1,0.2997,0.1000,0.0000
2,0.3000,0.1375,0.0731
3,0.3000,0.1325,0.0372
4,0.2994,0.1325,0.0310
5,0.2997,0.1350,0.0273
6,0.2995,0.1250,0.0253
7,0.2995,0.1187,0.0244
8,0.2995,0.1138,0.0240
9,0.2994,0.1117,0.0239
10,0.2993,0.1133,0.0243
11,0.2998,0.0883,0.0250
12,0.2992,0.0800,0.0262
13,0.2999,0.0733,0.0276
14,0.2999,0.0617,0.0289
15,0.2998,0.0500,0.0302
16,0.2998,0.0500,0.0308
17,0.2988,0.0450,0.0323
18,0.2984,0.0225,0.0337
19,0.2995,0.0400,0.0347
20,0.2995,0.0150,0.0365
21,0.2995,0.0200,0.0371
22,0.2995,0.0225,0.0384
23,0.2995,0.0225,0.0393
24,0.2995,0.0150,0.0399
25,0.2995,0.0200,0.0401
26,0.2995,0.0150,0.0406
27,0.2995,0.0225,0.0408
28,0.2995,0.0250,0.0419
29,0.2995,0.0100,0.0428
30,0.2995,0.0250,0.0434
31,0.2995,0.0100,0.0442
32,0.2995,0.0125,0.0444
33,0.2995,0.0100,0.0451
34,0.2995,0.0100,0.0459
35,0.2995,0.0075,0.0465
36,0.2995,0.0050,0.0473
37,0.2995,0.0150,0.0473
38,0.2995,0.0225,0.0479
39,0.2995,0.0125,0.0487
40,0.2995,0.0025,0.0499
41,0.2995,0.0150,0.0505
42,0.2995,0.0075,0.0503
```

### Directedness Trace
```
2,1.0000
3,1.0000
4,1.0000
5,1.0000
6,1.0000
7,1.0000
8,1.0000
9,1.0000
10,0.9999
11,0.9999
12,0.9994
13,0.9990
14,0.9825
15,0.9967
16,0.9788
17,0.9811
18,0.9751
19,0.9514
20,0.9346
21,0.9476
22,0.9371
23,0.9142
24,0.8891
25,0.8642
26,0.8458
27,0.8557
28,0.8513
29,0.8641
30,0.8537
31,0.8392
32,0.8366
33,0.8484
34,0.8399
35,0.8459
36,0.8313
37,0.8077
38,0.8098
39,0.8023
40,0.7983
41,0.7837
42,0.7978
```

### OpenKineticsPredictor Change Magnitude Quartiles
```
Gen,Q25,Q50,Q75
1,0.0000,0.0000,0.0000
2,0.0000,0.0000,0.0000
3,0.0000,0.0000,0.0000
4,0.0000,0.0000,0.0000
5,0.0000,0.0000,0.0000
6,0.0000,0.0000,0.0000
7,0.0000,0.0000,0.0000
8,0.0000,0.0000,0.0000
9,0.0000,0.0000,0.0000
10,0.0000,0.0000,0.0000
11,0.0000,0.0000,0.0000
12,0.0000,0.0000,0.0000
13,0.0000,0.0000,0.0000
14,0.0000,0.0000,0.0000
15,0.0000,0.0000,0.0000
16,0.0000,0.0000,0.0000
17,0.0000,0.0000,0.0000
18,0.0000,0.0000,0.0000
19,0.0000,0.0000,0.0000
20,0.0000,0.0000,0.0000
21,0.0000,0.0000,0.0000
22,0.0000,0.0000,0.0000
23,0.0000,0.0000,0.0000
24,0.0000,0.0000,0.0000
25,0.0000,0.0000,0.0000
26,0.0000,0.0000,0.0000
27,0.0000,0.0000,0.0000
28,0.0000,0.0000,0.0000
29,0.0000,0.0000,0.0000
30,0.0000,0.0000,0.0000
31,0.0000,0.0000,0.0000
32,0.0000,0.0000,0.0000
33,0.0000,0.0000,0.0000
34,0.0000,0.0000,0.0000
35,0.0000,0.0000,0.0000
36,0.0000,0.0000,0.0000
37,0.0000,0.0000,0.0000
38,0.0000,0.0000,0.0000
39,0.0000,0.0000,0.0000
40,0.0000,0.0000,0.0000
41,0.0000,0.0000,0.0000
42,0.0000,0.0000,0.0000
```

### brenda Change Magnitude Quartiles
```
Gen,Q25,Q50,Q75
1,0.0000,0.0000,0.0000
2,0.0000,0.0000,0.0000
3,0.0000,0.0000,0.0000
4,0.0000,0.0000,0.0000
5,0.0000,0.0000,0.0000
6,0.0000,0.0000,0.0000
7,0.0000,0.0000,0.0000
8,0.0000,0.0000,0.0000
9,0.0000,0.0000,0.0000
10,0.0000,0.0000,0.0000
11,0.0000,0.0000,0.0000
12,0.0000,0.0000,0.0000
13,0.0000,0.0000,0.0000
14,0.0000,0.0000,0.0000
15,0.0000,0.0000,0.0000
16,0.0000,0.0000,0.0000
17,0.0000,0.0000,0.0000
18,0.0000,0.0000,0.0000
19,0.0000,0.0000,0.0000
20,0.0000,0.0000,0.0000
21,0.0000,0.0000,0.0000
22,0.0000,0.0000,0.0000
23,0.0000,0.0000,0.0000
24,0.0000,0.0000,0.0000
25,0.0000,0.0000,0.0000
26,0.0000,0.0000,0.0000
27,0.0000,0.0000,0.0000
28,0.0000,0.0000,0.0000
29,0.0000,0.0000,0.0000
30,0.0000,0.0000,0.0000
31,0.0000,0.0000,0.0000
32,0.0000,0.0000,0.0000
33,0.0000,0.0000,0.0000
34,0.0000,0.0000,0.0000
35,0.0000,0.0000,0.0000
36,0.0000,0.0000,0.0000
37,0.0000,0.0000,0.0000
38,0.0000,0.0000,0.0000
39,0.0000,0.0000,0.0000
40,0.0000,0.0000,0.0000
41,0.0000,0.0000,0.0000
42,0.0000,0.0000,0.0000
```

### custom Change Magnitude Quartiles
```
Gen,Q25,Q50,Q75
1,0.0000,0.0000,0.0000
2,0.0000,0.0000,0.0000
3,0.0000,0.0000,0.0000
4,0.0000,0.0000,0.0000
5,0.0000,0.0000,0.0000
6,0.0000,0.0000,0.0000
7,0.0000,0.0000,0.0000
8,0.0000,0.0000,0.0000
9,0.0000,0.0000,0.0000
10,0.0000,0.0000,0.0000
11,0.0000,0.0000,0.0000
12,0.0000,0.0000,0.0000
13,0.0000,0.0000,0.0000
14,0.0000,0.0000,0.0000
15,0.0000,0.0000,0.0000
16,0.0000,0.0000,0.0000
17,0.0000,0.0000,0.0000
18,0.0000,0.0000,0.0000
19,0.0000,0.0000,0.0000
20,0.0000,0.0000,0.0000
21,0.0000,0.0000,0.0000
22,0.0000,0.0000,0.0000
23,0.0000,0.0000,0.0000
24,0.0000,0.0000,0.0000
25,0.0000,0.0000,0.0000
26,0.0000,0.0000,0.0000
27,0.0000,0.0000,0.0000
28,0.0000,0.0000,0.0000
29,0.0000,0.0000,0.0000
30,0.0000,0.0000,0.0000
31,0.0000,0.0000,0.0000
32,0.0000,0.0000,0.0000
33,0.0000,0.0000,0.0000
34,0.0000,0.0000,0.0000
35,0.0000,0.0000,0.0000
36,0.0000,0.0000,0.0000
37,0.0000,0.0000,0.0000
38,0.0000,0.0000,0.0000
39,0.0000,0.0000,0.0000
40,0.0000,0.0000,0.0000
41,0.0000,0.0000,0.0000
42,0.0000,0.0000,0.0000
```

### unlabelled Change Magnitude Quartiles
```
Gen,Q25,Q50,Q75
1,0.0000,0.0000,0.0000
2,0.0000,0.0000,0.0000
3,0.0000,0.0000,0.0000
4,0.0000,0.0000,0.0000
5,0.0000,0.0000,0.0000
6,0.0000,0.0000,0.0000
7,0.0000,0.0000,0.0000
8,0.0000,0.0000,0.0000
9,0.0000,0.0000,0.0000
10,0.0000,0.0000,0.0000
11,0.0000,0.0000,0.0000
12,0.0000,0.0000,0.0000
13,0.0000,0.0000,0.0000
14,0.0000,0.0000,0.0000
15,0.0000,0.0000,0.0000
16,0.0000,0.0000,0.0000
17,0.0000,0.0000,0.0000
18,0.0000,0.0000,0.0000
19,0.0000,0.0000,0.0000
20,0.0000,0.0000,0.0000
21,0.0000,0.0000,0.0000
22,0.0000,0.0000,0.0000
23,0.0000,0.0000,0.0000
24,0.0000,0.0000,0.0000
25,0.0000,0.0000,0.0000
26,0.0000,0.0000,0.0000
27,0.0000,0.0000,0.0000
28,0.0000,0.0000,0.0000
29,0.0000,0.0000,1.0721
30,0.0000,0.0000,1.1607
31,0.0000,0.0000,1.1381
32,0.0000,0.0000,1.2500
33,0.0000,0.0000,1.2682
34,0.0000,0.0000,1.3038
35,0.0000,0.0000,1.3395
36,0.0000,0.0000,1.3866
37,0.0000,0.0000,1.3560
38,0.0000,0.0000,1.3682
39,0.0000,0.0000,1.3589
40,0.0000,0.0000,1.4518
41,0.0000,0.0000,1.4989
42,0.0000,0.0000,1.6408
```

</details>

---

## Analysis Request

Based on this diagnostic report, please provide:

1. **Overall assessment** of convergence quality
2. **Specific parameter recommendations** with justification:
   - Which shrinkThrSource values to adjust?
   - Which varianceCapSource values to adjust?
   - Should sparsityThreshold change?
   - Should maxRMSEplateau change?
3. **Priority order** for changes (what to fix first)
4. **Expected impact** of recommended changes

Please be specific with numbers (e.g., "increase shrinkThrSource(2) from 10.0 to 12.0") rather than general advice.
