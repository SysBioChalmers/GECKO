# Upgrading from GECKO 3 to GECKO 4 (MATLAB)

GECKO 4 was developed alongside **geckopy**, the Python port of GECKO,
and several changes were made so that the two toolboxes build and store
ecModels identically. This guide is for users with existing GECKO 3
scripts and models. Most changes are transparent; the ones that can
affect your code or results are listed first.

## Your existing models still load

GECKO 4 reads GECKO 3 ecModels. You do not need to convert anything by
hand:

- **YAML.** RAVEN's `readYAMLmodel` still reads the old
  `---` / `!!omap` files. `writeYAMLmodel` now emits the new canonical
  format (see below), so **loading an old model and saving it upgrades
  it automatically**.
- **Protein direction.** `loadEcModel` detects the old reverse-direction
  protein reactions and flips them to the new forward convention on load
  (`flipLegacyProtDirection`). Saving then writes them in the forward
  convention.

In other words: load a GECKO 3 model in GECKO 4, save it, and it is now
a GECKO 4 model. (geckopy does the same on its side, so the file is then
also directly usable from Python.)

## Changes that can affect your scripts

### 1. Protein-pool reactions now run forward

This is the most visible change. In GECKO 3 the `usage_prot_<id>`
reactions and `prot_pool_exchange` ran in **reverse**: protein was drawn
as a negative flux, the reactions had `lb = -1000` (`ub = 0`), and the
available protein capacity lived in the **negative lower bound**.

GECKO 4 uses the **forward** direction: positive flux, `ub = 1000`
(`lb = 0`), capacity in the **upper bound**. Adjust scripts that:

- read usage from the pool: `sol.x(prot_pool_exchange)` is now a
  **positive** number;
- set the pool size: use `setParam(model,'ub','prot_pool_exchange',X)`
  instead of `setParam(model,'lb','prot_pool_exchange',-X)`;
- constrain a single enzyme: set the upper bound (`ub = conc`) instead
  of the lower bound (`lb = -conc`);
- read enzyme usage: it is the flux / `ub`, not `-lb`.

Models built with GECKO 3 are flipped automatically on load (above), so
you only need to update **code** that hard-codes the reverse convention.

### 2. Corrected results from bug fixes

Several GECKO 3 bugs were fixed. For these steps GECKO 4 produces
different — and more correct — numbers than GECKO 3, so don't be alarmed
if a re-run doesn't match an old one:

- `sigmaFitter` now returns the model fitted to the **optimal** sigma,
  not the last trial value (sigma = 1.0).
- `getStandardKcat` now applies the **per-subsystem** mean kcat; GECKO 3
  fell back to the global standard kcat unless every subsystem matched.
- `getECfromGEM` now actually assigns EC codes; the GECKO 3 validation
  regex silently discarded every EC string.
- `findECInDB` no longer emits duplicate EC codes (e.g.
  `"1.1.1.1;1.1.1.1"`).
- `findMaxValue`'s wildcard EC branch now matches real codes.
- `fuzzyKcatMatching` no longer escalates wildcards past the
  fully-wildcarded form (which could error out).
- `writeDLKcatInput` no longer indexes out of bounds when `ecRxns` is a
  subset.
- `readDLKcatOutput` matches substrate names case-insensitively.

### 3. BRENDA database files renamed and reformatted

GECKO 3 shipped `max_KCAT.txt` / `max_SA.txt` / `max_MW.txt` (5 columns,
`EC`-prefixed codes, `organism//taxonomy//keggcode` triples). GECKO 4
ships `max_kcat.tsv` / `max_sa.tsv` / `max_mw.tsv`, refreshed from the
BRENDA bulk JSON, with bare EC numbers, plain organism names, a `#`
release header, and a `references` (PMID) column. `loadBRENDAdata` reads
the new files. Remove any stale `max_*.txt` from your project's `data/`
folder.

### 4. OpenKineticsPredictor: submit/fetch instead of write/read

`writeOpenKineticsPredictorInput` and `readOpenKineticsPredictorOutput`
are **removed**. They required manually uploading a CSV at the OKP
website and downloading the result. GECKO 4 replaces them with:

- `submitOpenKineticsPredictor` — builds the input and submits a
  prediction job via the OKP REST API; stores the job id in
  `data/OKP_job.txt`.
- `fetchOpenKineticsPredictor` — checks the job, downloads the result on
  completion, and parses it into a `kcatList`.

You need a free OKP API key (generated on the website). Provide it via
the function argument, the `OKP_API_KEY` environment variable, or
`data/okpApiKey.txt` (do not put it in the model adapter). The per-row
`kcatList.kcatSource` now reports the actual provenance (e.g. `CataPro`,
`BRENDA`, `Sabio-RK`).

### 5. `ec.genes` is sorted alphabetically

`makeEcModel` now sorts `ec.genes` (and the parallel `ec.enzymes`,
`ec.mw`, `ec.sequence`) alphabetically for a stable order across runs.
The values are unchanged, but the **order** of entries in the `ec`
struct differs from GECKO 3.

## Additive changes (won't break anything)

- **KEGG as an additional protein/EC source.** `makeEcModel` (stage 7)
  can fall back to KEGG for genes UniProt cannot match, and
  `getECfromDatabase` can fall back to KEGG for EC numbers. A new
  optional `params.okp` / KEGG configuration controls this.
- **Download helpers split out.** `downloadKEGG.m` and
  `downloadUniProt.m` are now standalone functions; `loadDatabases`
  still auto-downloads when a file is missing, so existing calls are
  unaffected.
- **Bayesian kcat tuning improvements** in `bayesianSensitivityTuning`
  (posterior-sample reporting, optional pruning of insignificant kcats,
  refined acceptance strategy) and updated documentation.
- **Refreshed reference data** (paxDB, CatPred predictions).

## YAML format details

GECKO 4 / RAVEN write the canonical, cobrapy-compatible YAML: a flat
mapping (no outer `---` sequence, no `!!omap` tags) with the
GECKO-specific data under top-level keys. RAVEN still reads the old
format, so nothing breaks. The benefit: ecModels now exchange directly
between GECKO 4 (MATLAB) and geckopy (Python) with no conversion step.

## Moving to Python?

If you want to run the same workflow in Python, see the geckopy
migration guide (`docs/migrating_from_gecko_matlab.md` in the geckopy
repository), which maps GECKO 4 functions and conventions to their
geckopy equivalents.
