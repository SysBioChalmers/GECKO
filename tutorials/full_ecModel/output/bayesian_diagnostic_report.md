# Bayesian Kcat Tuning Diagnostic Report

**Generated:** 2026-06-10 11:44:47

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
shrinkThrDefault = 5.0  % Unlabelled: σ deviations for full update
shrinkThrSource = [3.0; 10.0; 12.0]  % Per-source shrinkage

varianceCapDefault = 1.50  % Unlabelled: max variance growth
varianceCapSource = [1.50; 1.10; 1.05]  % Per-source caps

forcePriorThrDefault = 3.5  % Unlabelled: snap to prior if dev < this
forcePriorThrSource = [2.5; 11.0; 13.0]  % Per-source force thresholds

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

- **Total generations:** 31
- **Converged:** No
- **Final RMSE:** 0.8715
- **Initial RMSE:** 8.6001
- **RMSE reduction:** 89.9%
- **Total parameters:** 4834

### Convergence Pattern

- **Early phase (gen 1-10):** 81.7% RMSE reduction
- **Mid phase (gen 11-20):** 9.7% RMSE reduction
- **Late phase (gen 21+):** 2.9% RMSE reduction
- **Late-stage directedness:** 0.78 (1=converging, 0.5=random, <0.3=oscillating)

---

## Per-Source Analysis

### OpenKineticsPredictor

**Summary:**
- Total parameters: 1175
- Active (shrinkWeight > 0.3): 478 (40.7%)
- Near prior (|σ_post - σ_prior| < 0.1): 697 (59.3%)
- Mean deviation from prior: 1.72 σ
- Variance ratio (σ_post / σ_prior): 1.20

**Change magnitude (|log(final/prior)|):**
- Q25: 0.000
- Median: 0.000
- Q75: 0.872

### brenda

**Summary:**
- Total parameters: 3198
- Active (shrinkWeight > 0.3): 41 (1.3%)
- Near prior (|σ_post - σ_prior| < 0.1): 3198 (100.0%)
- Mean deviation from prior: 0.27 σ
- Variance ratio (σ_post / σ_prior): 1.00

**Change magnitude (|log(final/prior)|):**
- Q25: 0.000
- Median: 0.000
- Q75: 0.000

### custom

**Summary:**
- Total parameters: 207
- Active (shrinkWeight > 0.3): 0 (0.0%)
- Near prior (|σ_post - σ_prior| < 0.1): 207 (100.0%)
- Mean deviation from prior: 0.00 σ
- Variance ratio (σ_post / σ_prior): 1.00

**Change magnitude (|log(final/prior)|):**
- Q25: 0.000
- Median: 0.000
- Q75: 0.000

### unlabelled

**Summary:**
- Total parameters: 254
- Active (shrinkWeight > 0.3): 53 (20.9%)
- Near prior (|σ_post - σ_prior| < 0.1): 201 (79.1%)
- Mean deviation from prior: 1.01 σ
- Variance ratio (σ_post / σ_prior): 1.10

**Change magnitude (|log(final/prior)|):**
- Q25: 0.000
- Median: 0.000
- Q75: 0.000

---

## Key Diagnostic Traces

### RMSE Evolution (every 5 generations)
```
Gen | RMSE   | Improvement
----|--------|------------
  0 | 8.6001 | -
  5 | 2.5657 | 12.2%
 10 | 1.2807 | 2.9%
 15 | 0.9933 | 0.0%
 20 | 0.8972 | 0.0%
 25 | 0.8715 | 0.0%
 30 | 0.8715 | 0.0%
 31 | 0.8715 | 0.0%
```

### Acceptance Rates (every 5 generations)
```
Gen | Overall | New Proposals | Boundary Violations
----|---------|---------------|--------------------
  5 |  30.0% |  13.0%       |   2.7%
 10 |  29.9% |   9.8%       |   2.4%
 15 |  30.0% |   5.5%       |   3.0%
 20 |  29.9% |   3.0%       |   3.6%
 25 |  29.9% |   2.2%       |   4.0%
 30 |  29.9% |   1.0%       |   4.4%
```

### Directedness & Diversity (every 5 generations)
```
Gen | Directedness | Diversity | Sparsity %
----|--------------|-----------|----------
  5 | 1.00        | 3.118     |  99.2%
 10 | 1.00        | 5.634     |  99.1%
 15 | 0.94        | 7.553     |  96.5%
 20 | 0.88        | 8.315     |  92.9%
 25 | 0.79        | 8.601     |  89.9%
 30 | 0.72        | 8.972     |  88.4%
```

---

## Potential Issues Detected

### ⚠️ WARNING: Very Low Proposal Acceptance
- Late-stage new proposal acceptance: 1.7% (should be 5-15%)
- **Likely cause:** Proposals too aggressive
- **Recommendation:** Proposal width adaptation may need adjustment

---

## Detailed Diagnostic Data

<details>
<summary>Click to expand full diagnostic traces</summary>

### Complete RMSE Trace
```
0,8.600052
1,6.994339
2,4.125616
3,3.607907
4,2.921939
5,2.565691
6,1.941535
7,1.724222
8,1.378978
9,1.318582
10,1.280676
11,0.993332
12,0.993332
13,0.993332
14,0.993332
15,0.993332
16,0.993332
17,0.897186
18,0.897186
19,0.897186
20,0.897186
21,0.897186
22,0.871475
23,0.871475
24,0.871475
25,0.871475
26,0.871475
27,0.871475
28,0.871475
29,0.871475
30,0.871475
31,0.871475
```

### Acceptance Rate Trace
```
Gen,Overall,NewProposals,BoundaryViolations
1,0.2997,0.1000,0.0000
2,0.3000,0.1375,0.0731
3,0.3000,0.1300,0.0371
4,0.2994,0.1263,0.0312
5,0.2997,0.1300,0.0274
6,0.2995,0.1288,0.0255
7,0.2995,0.1225,0.0244
8,0.2995,0.1113,0.0238
9,0.2994,0.1017,0.0237
10,0.2993,0.0983,0.0241
11,0.2998,0.0933,0.0251
12,0.2992,0.0750,0.0261
13,0.2999,0.0733,0.0273
14,0.2999,0.0483,0.0288
15,0.2998,0.0550,0.0296
16,0.2998,0.0375,0.0315
17,0.2988,0.0425,0.0325
18,0.2984,0.0250,0.0336
19,0.2995,0.0325,0.0347
20,0.2995,0.0300,0.0363
21,0.2995,0.0300,0.0372
22,0.2995,0.0175,0.0377
23,0.2995,0.0250,0.0383
24,0.2995,0.0175,0.0389
25,0.2995,0.0225,0.0397
26,0.2995,0.0150,0.0406
27,0.2995,0.0150,0.0418
28,0.2995,0.0225,0.0429
29,0.2995,0.0200,0.0441
30,0.2995,0.0100,0.0441
31,0.2995,0.0050,0.0449
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
11,0.9978
12,0.9861
13,0.9787
14,0.9505
15,0.9359
16,0.9297
17,0.9071
18,0.8852
19,0.8865
20,0.8820
21,0.8798
22,0.8858
23,0.8370
24,0.8251
25,0.7865
26,0.7772
27,0.7762
28,0.7459
29,0.7305
30,0.7220
31,0.7130
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
21,0.0000,0.0000,0.5458
22,0.0000,0.0000,0.6515
23,0.0000,0.0000,0.6486
24,0.0000,0.0000,0.7185
25,0.0000,0.0000,0.7468
26,0.0000,0.0000,0.7656
27,0.0000,0.0000,0.7818
28,0.0000,0.0000,0.7891
29,0.0000,0.0000,0.8207
30,0.0000,0.0000,0.8396
31,0.0000,0.0000,0.8723
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
29,0.0000,0.0000,0.0000
30,0.0000,0.0000,0.0000
31,0.0000,0.0000,0.0000
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
