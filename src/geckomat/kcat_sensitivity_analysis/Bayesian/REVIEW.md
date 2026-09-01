# Bayesian kcat-tuning: review and path forward

## Context

This branch (`fix/bayesianTuning`) was a 2025 attempt at incrementally fixing the
Bayesian kcat-tuning implementation. In the meantime, `develop4`'s
`bayesianSensitivityTuning.m` and `abc_max.m` were independently rewritten at much
larger scale (roughly 900 changed lines in `bayesianSensitivityTuning.m` alone since
this branch's fork point), and appear to have absorbed this branch's concrete fixes
already: the biomass-correction commit at the tip of this branch is present in
`develop4` verbatim, including comments.

Given that, this branch is not a useful vehicle for further incremental patches. The
plan instead is to develop the method further in Python first, then port the
finished design back into MATLAB, rather than keep iterating on the MATLAB
implementation directly. This document records what a fresh look at the current
`develop4` implementation turned up, as input for that Python redesign.

**Scope note:** the findings below are against `develop4`'s current
`src/geckomat/kcat_sensitivity_analysis/Bayesian/*`, not against this branch's own
(superseded) versions of those files.

## High severity — confirmed bugs

### 1. The documented `kcatSigmaLog` override crashes

In `bayesianSensitivityTuning.m`:

```matlab
if isempty(kcatSigmaLog)
    sigma0log                 = sigma0logDefault * ones(size(ecModel.ec.kcat));
    sigma0log(uniqKcatParams) = sigma0logSource(kcatSourceIdx(uniqKcatParams));
end
```

`sigma0log` is only assigned inside this branch, yet it's referenced throughout the
rest of the function (`sigmaLogTrace = sigma0log`, `devFromPrior = ... ./ sigma0log`,
etc.). The docstring explicitly advertises a custom `kcatSigmaLog` as supported —
its own Examples section shows
`bayesianSensitivityTuning(ecModel, kcatSigmaLog, modelAdapter)` — but calling it
that way skips the block entirely and crashes with "Undefined variable 'sigma0log'".
Even if that crash were fixed, nothing ever does `sigma0log = kcatSigmaLog`, so the
override isn't wired up to have any effect either. The feature is broken twice over.

### 2. Hardcoded personal debug path in the sampling loop

`bayesianSensitivityTuning.m`, inside the `parfor` that evaluates each generation's
proposals:

```matlab
parfor j = 1:N
    %Uncomment for Vera, remove before release
    addpath('/apps/Arch/software/Gurobi/11.0.2-GCCcore-12.3.0/matlab/')
```

The comment says "remove before release"; it wasn't. Runs unconditionally on every
worker, for every sample, every generation.

### 3. Asymmetric empty-guard in `loadBayesianData.m`

```matlab
bayData.fluxData = loadFluxData(fullfile(basePath,'data','bayesianFluxData.tsv'));
if ~isempty(bayData.fluxData)
    bayData.fluxData.biomass = modelAdapter.params.bioRxn;
end
bayData.maxGrate = loadFluxData(fullfile(basePath,'data','bayesianMaxGrowth.tsv'));
bayData.maxGrate.biomass = modelAdapter.params.bioRxn;   % no guard
```

`loadFluxData` returns `[]` when its file is missing or unreadable (caught in its
own try/catch). The `fluxData` case handles that correctly. The `maxGrate` case
doesn't: `emptyArray.biomass = ...` is legal MATLAB and silently turns `[]` into a
1x1 struct containing only `.biomass`. That struct then passes `abc_max.m`'s
`~isempty(bayData.maxGrate)` check, so `rmsecal` runs against it and fails on a
missing field instead of the graceful skip a genuinely-absent
`bayesianMaxGrowth.tsv` should produce.

## Medium severity

### 4. `addCarbonNum.m` hardcodes the biomass reaction by name

```matlab
model.excarbon(strcmp(model.rxnNames,'growth')) = 41;
```

Everywhere else in this pipeline the biomass reaction comes from the model adapter
(`data.biomass` / `modelAdapter.params.bioRxn`). Here it's a literal string match.
If a model's biomass reaction isn't named exactly `'growth'`, this is a silent
no-op — `excarbon` for biomass stays 0 instead of 41, quietly breaking the RMSE
carbon-weighting for any adapter that doesn't use that exact name.

### 5. Silent array-shrinking landmine in `getcarbonnum` (nested in `addCarbonNum.m`)

```matlab
EXmetsIdx = zeros(length(exrxn),1);
for k = 1:length(EXmets(1,:))
    EXmetsIdx(k) = find(EXmets(:,k));
end
```

If `EXmets(:,k)` ever has zero nonzero entries, `find` returns `[]`, and
`EXmetsIdx(k) = []` deletes element k instead of erroring — silently shrinking the
array and shifting every later index. Currently unreachable in practice (inputs
come from `getExchangeRxns`, which guarantees exactly one nonzero entry per
exchange reaction), but fragile: a future change to exchange-reaction
classification would corrupt data silently rather than raise a clear error.

## Low severity / code quality

- **`getrSample.m` is dead code.** Nothing in the current codebase calls it.
  `bayesianSensitivityTuning.m` samples generation 1 with an inline `lognrnd` call
  and later generations via `proposeLowRankMixture`; the 'uniform' method in
  `getrSample` isn't used anywhere.
- **`abc_max.m`'s `rmsecal`** recomputes `ismember(data.conds, data.exchMets)` (and
  its missing-carbon-source validation) on every inner-loop iteration despite not
  depending on the loop variable — wasteful in a real hot path (once per sample,
  hundreds of samples per generation).
- **Magic number `99`** for the failed-model RMSE penalty
  (`rmseList(isnan(rmseList)) = 99;`) is undocumented as to why 99 is "clearly
  worse" than any real RMSE.
- **`addCarbonNum.m` vendors a ~250-line generic chemical-formula parser**
  (`getElementalComposition`, attributed "Siu Hung Joshua Chan May 2017") nested
  three levels deep inside a single-purpose file. Worth its own top-level file, and
  may duplicate a utility that already exists elsewhere in the COBRA/RAVEN
  ecosystem.

## Path forward

Given the scale of rework this implies, and that this area has already been
independently rewritten at least once since this branch forked, the plan is:

1. Prototype the ABC-SMC kcat-tuning method fresh in Python (geckopy), where
   iteration is cheaper and the statistics are easier to unit-test in isolation.
2. Once the method design is validated there, port the finished implementation
   back into MATLAB, rather than continuing to patch the current MATLAB code
   in place.

This branch and this document are the starting reference point for that Python
work — not a queue of MATLAB patches to apply directly.
