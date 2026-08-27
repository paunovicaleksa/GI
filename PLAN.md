# Nanoplastics scRNA-seq Pipeline — Build Plan

Working plan for the Transcriptomics project. Project scope, task list, and
methodological constraints live in `CLAUDE.md`; this file records **how** we
intend to execute them and **why** each choice was made. Update it as decisions
change — it should stay in sync with the notebook.

## Context

We analyze 4 PBMC samples exposed to carboxylated polystyrene nanoparticles
(40nm, 200nm, 40+200nm mixture, control), all from one donor. We follow a
SanBomics YouTube Scanpy tutorial as the "how-to-code-it" layer, while `CLAUDE.md`
is the "what/why" layer. The intent is a learning exercise: mirror the tutorial
closely, and deviate only where `CLAUDE.md` or a practical constraint requires it.

`single_cell_analysis_complete_class.ipynb` (190 cells) is that tutorial — the
SanBomics COVID-lung walkthrough covering doublet removal via scVI/SOLO, QC,
normalization, single-sample clustering, multi-sample integration via scVI,
marker-based manual annotation, composition counts, DE via diffxpy/scVI, GO
enrichment via gseapy, and gene-signature scoring. **It is read-only reference
material — we never edit it.**

### What the data actually contains

Inspecting the 4 `.h5ad` files turned up something `CLAUDE.md` did not anticipate:
they are **already Seurat-processed**, not merely UMI-deduplicated.

| Slot | Contents |
|---|---|
| `layers['counts']` | Raw integer counts — our real starting point |
| `X` | Already normalized + log-transformed (max ≈ 5.5), float64 sparse |
| `obs` | `nCount_RNA`, `nFeature_RNA`, `seurat_clusters`, `RNA_snn_res.0.3` |
| `obs` (annotation) | `predicted.id`, `predicted.celltype`, `prediction.score.*` |
| `obsm` | `X_pca` (50), `X_ref.pca` (30), `X_ref.umap`, `X_umap` |
| `var` | Seurat `FindVariableFeatures` outputs |
| `uns` / `raw` | empty / None |

Cell counts: 8,729 (S1) / 12,676 (S2) / 6,157 (S3) / 6,516 (S4), across ~22k genes.
Mitochondrial genes use the `MT-` uppercase human convention. **No sample or
condition column exists in any file** — we add it ourselves at load time.

The `X_ref.pca` / `X_ref.umap` slots plus the `prediction.score.*` naming
(Seurat's `TransferData()` output) mean these files already carry **a completed
Azimuth run against the PBMC reference**, done by someone else.

The companion `*_CoDi_KLD.csv` files are a second, independent per-cell cell-type
prediction (`CoDi`, plus `_dist` / `_contrastive` variants with confidence
scores), one row per cell. `CLAUDE.md` designates these optional.

## Empirical findings carried forward

Observations produced by the pipeline itself that are **not yet explained**, and that
later sections must test rather than assume. Kept here, near the top, so they are not
buried in the per-section outcome notes where they were first recorded.

### 200nm is the least similar sample to control — the first size-dependent signal

Before correcting anything, scanorama (Section 8) scores how much genuine population
overlap it finds between every pair of samples. Measured, at `knn=20`:

| pair | alignment score |
|---|---|
| mixture ↔ control | 0.925 |
| 40nm ↔ mixture | 0.904 |
| 40nm ↔ 200nm | 0.809 |
| 200nm ↔ mixture | 0.809 |
| 40nm ↔ control | 0.731 |
| **200nm ↔ control** | **0.682 — lowest of all six** |

**The ordering is size-dependent**: control resembles the mixture most (0.925), then 40nm
(0.731), and 200nm least (0.682). This is the first quantitative hint of a particle-size
effect anywhere in the pipeline, and it arrives **before any differential expression
test**, from a method whose only job was batch correction.

**Why this is a hypothesis and not a result:**
- The score reflects **overall population structure** — roughly, how well cell types
  correspond between two samples. It cannot distinguish "different cell-type proportions"
  from "same cell types, shifted expression". Tasks 4 and 5 separate those; this does not.
- **200nm is also the shallowest-sequenced sample** (median 4,390 counts against control's
  7,736). Lower depth degrades neighbour matching on its own, so part of the 0.682 may be
  technical. With one donor that cannot be separated here.
- It is one number per pair, with no null distribution and no replication.

**Where it gets tested properly:** Section 12 (does 200nm differ in cell-type
*proportions*?), Section 13 (does it differ in *expression* within cell types?), and
Section 15, which has to decide what is 40nm-specific, 200nm-specific, shared, or
mixture-only. If those independently reproduce this ordering, the alignment scores become
corroborating evidence. If they do not, this was depth.

**Counter-signal to hold alongside it:** the mixture contains 200nm particles yet
resembles control *most* (0.925). A simple dose-response story does not predict that.

### Two further signals in the same direction (Section 9)

Clustering at `resolution=0.4` produced 16 clusters whose sample composition is otherwise
near-proportional — almost every enrichment value sits between 0.85 and 1.35, which is the
strongest confirmation yet that Section 8 corrected batch without erasing biology. **Three
clusters break that pattern, and all three break it the same way:**

| cluster | cells | 40nm | 200nm | control | mixture | max/min |
|---|---|---|---|---|---|---|
| **11** | **1,264** | 0.53 | **1.59** | **0.52** | 0.96 | **3.1x** |
| **12** | 17 | 0.69 | **1.58** | **0.31** | 0.98 | 5.1x |
| **15** | 5 | 0.00 | **2.15** | 0.00 | 1.11 | — |

Enriched in 200nm, depleted in control. **Cluster 11 matters most**: 1,264 cells is enough
for real statistical power in Sections 12-13, and it is visually an isolated island on the
UMAP. Cluster 15 is five cells and can only ever be an anecdote.

### Cluster 15 — an acute inflammatory program in 200nm-containing samples only

The five-cell cluster that survives at every resolution was expected to be residual
doublets, a contaminant, or dying cells. On the evidence it is none of those.

```
samples:      200nm 4, mixture 1, 40nm 0, control 0
median counts 7,197  (dataset 5,454)
median genes  2,808  (dataset 2,118)
median mito%  2.68   (dataset 3.88)   <- BELOW average
```

Top genes against the rest of the dataset: `IL8` 5.92 vs 1.13, `SERPINB2` 4.65 vs 0.36,
`IL1B` 4.50 vs 0.44, `SOD2` 4.29 vs 0.76, `EREG` 3.92, `THBS1` 3.91, `CXCL2` 3.77,
`NAMPT` 3.54, `CCL3` 3.33, `CXCL1` 3.34, `IL6` 3.14 vs 0.13.

That is a coherent NF-kB-driven cytokine/chemokine program — and `CLAUDE.md` named these
genes in advance ("TNF/IL-6/IL-8 are literature-plausible candidates … check pathway-level
enrichment (NF-kB/cytokine signaling)", concentrated in myeloid cells). The **below-average
mitochondrial content rules out dying cells**; these are transcriptionally hyperactive, not
damaged.

**Why it is still only a hypothesis:**
- **n = 5.** No statistical test is possible.
- **High counts + high genes is also the doublet signature.** The counter-argument — that a
  doublet shows two lineages' markers while this shows one coherent program — needs
  Section 10's markers to verify rather than assert.
- **These are also neutrophil genes.** `IL8`, `CXCL1`, `CXCL2`, `SERPINB2` are high in
  granulocytes. PBMC preps mostly exclude them, but not perfectly.

**One thing it does settle:** the depth confound cannot explain *this* cluster. Its cells
have 7,197 median counts, well above the dataset median — shallow sequencing does not
produce unusually deep cells. That confound remains live for the alignment scores and for
cluster 11, but not here.

**`IL8` is the deferred symbol problem becoming load-bearing.** This dataset uses pre-2017
HGNC names, so MSigDB's `CXCL8` will not match our `IL8`. The gene at the top of the most
interesting result in the pipeline is exactly the one that will silently fail to match in
Sections 14 and 16. The old→new symbol map is no longer a tidy-up task.

### Section 10 identifies cluster 11 — it is the monocyte population

Marker annotation reframes the finding above. Cluster 11 is not an anomalous inflammatory
cluster: it is **the only monocyte population in the dataset** (`CD14`, `LYZ`, `S100A8`,
`S100A9`, `FCGR3A`, `MS4A7`, `CST3`; classical-monocyte panel 1.17 against < 0.20 for every
other panel). Clusters 12 (17 cells) and 15 (5 cells) are also myeloid — the expression
dendrogram, given no cell-type information, groups 11, 12 and 15 on one branch.

So the observation restates as:

> **Monocytes are roughly three times more abundant, as a fraction of cells, in 200nm
> (enrichment 1.59) than in control (0.52).**

That is a **composition** result and belongs to **Task 4 / Section 12**, not to DE. It lands
where `CLAUDE.md` predicted: *"Expect effects concentrated in monocytes/myeloid cells
(phagocytic — actually take up particles), not lymphocytes."*

**Do not conflate this with the DE question.** Cluster 11's high inflammatory-panel score
(3.37 against a cross-cluster mean of ~0.5) is substantially just what monocytes do —
`IL8` and `IL1B` are baseline monocyte genes. Whether 200nm monocytes are *more* inflammatory
than control monocytes is a within-cell-type comparison and is Section 13's job. Section 12
answers "how many", Section 13 answers "how different"; `CLAUDE.md` warns explicitly against
merging those two questions.

**Sample caveat:** monocytes are only **4.45%** of cells here, low against a typical PBMC
10-30%. Worth stating in the report as a property of this dataset.

### Other open observations

- **Doublet burden is flat at ~11–14% across all four samples** (Section 2) rather than
  scaling with cell count as the 10x loading model predicts. Test in Section 12 by
  computing each sample's homotypic fraction (Σpᵢ² over cell-type proportions) once
  labels exist.
- **The mito filter removes slightly more from exposed samples** (7.2–8.2%) than from
  control (6.5%) (Section 3) — small, but consistent in direction across all three
  exposures. Technical damage, or genuine exposure-induced mitochondrial stress? Feeds
  `CLAUDE.md` bonus #1.
- **Sequencing depth is confounded with exposure** (Section 3): control 7,736 > mixture
  6,276 > 40nm 5,436 > 200nm 4,390. Already accounted for in Section 7 (no `regress_out`)
  and Section 8 (`sigma=0.1`), and a stated limitation for every later comparison.

## Decisions and rationale

### Reprocessing scope
Reprocess fully from `layers['counts']`, ignoring Seurat's `X`, QC calls, HVGs,
PCA/UMAP, and `seurat_clusters` entirely. Treat the files as if they were a plain
counts matrix. This is the point of Tasks 1-2 and of the exercise: our own QC
threshold choices, normalization, HVG selection, batch correction, and clustering.

### Batch correction — scanorama
`CLAUDE.md` requires picking and justifying one method. We use **scanorama**
rather than `mnnpy`: `mnnpy` is unmaintained, and a dry-run install only proves
metadata resolves, not that its C/Cython extensions build against numpy 2.2.6.
scanorama is MNN-family (mutual-nearest-neighbors based internally) and actively
maintained. This replaces the tutorial's scVI integration step.

### Annotation (Task 3) — run by us, in Python
The pre-existing `predicted.celltype` / `prediction.score.*` / `X_ref.*` fields
are valid Azimuth output, but **inherited**. We demote them to **verification
only**, alongside `CoDi`, and run annotation ourselves. Two tools, both
`pip install`, both minutes to run:

- **`panhumanpy`** (Pan-Human Azimuth, Satija Lab) — `AzimuthNN` class operating
  directly on AnnData; model weights auto-download to `~/.cache/panhumanpy/`.
  This satisfies `CLAUDE.md`'s explicit "Azimuth" requirement: it *is* Azimuth,
  current version, no R needed. **Caveat to record in the notebook:** it is the
  pan-human model (23 tissues / 380 cell types), not the legacy PBMC-specific
  reference, which the Azimuth site states is no longer maintained.
- **`celltypist`** (`Immune_All_Low` / `Immune_All_High`) — immune-specialised
  second opinion. It fits this pipeline unusually well: its required input is
  log1p-normalised to 10,000 counts/cell, exactly what our Section 5 produces;
  and its `majority_voting` mode accepts **custom over-clustering**, so we feed
  it our own Leiden clusters rather than letting it re-cluster independently.

With our own marker-gene dot-plot annotation (also required by `CLAUDE.md`),
that gives three independent label sources to reconcile per cluster, plus two
external ones as verification.

**Routes considered and rejected**, recorded so we don't re-litigate them:
- The Azimuth Zenodo record (4546839) ships only `ref.Rds` (Seurat S4 — `pyreadr`
  cannot read S4 objects) and `idx.annoy` (vectors only, no labels or gene names).
  At 53 MB for 162k cells, `ref.Rds` is almost certainly an embedding + metadata
  rather than expression, i.e. built to feed Azimuth's own projection step.
- R/Seurat + Azimuth is explicitly permitted by `CLAUDE.md` and remains a valid
  escalation **if the pan-human labels prove too coarse on fine immune subtypes**.
  Its cost is environment debugging (Ubuntu ships R 4.1.2; Seurat compiles ~100
  dependencies; the h5ad↔Seurat conversion via SeuratDisk is the usual stall
  point) rather than anything that teaches us scRNA-seq.

### scVI avoided throughout
Scrublet (`sc.pp.scrublet`, built into scanpy) instead of scVI/SOLO for doublets;
scanorama instead of scVI for integration; pseudobulk DESeq2 via `pydeseq2`
instead of scVI/diffxpy for DE. This keeps one consistent story rather than scVI
in some steps and not others.

## Environment

The project must run on **both Linux and native Windows**. That drove three changes
from the original Linux-only setup; each is recorded below with its reason.

### Python 3.12 on both platforms

This plan originally specified the Linux venv's 3.10.12. That version cannot be
matched on Windows — python.org never shipped a Windows installer past **3.10.11**,
because 3.10 went security-only (source-release) from 3.10.12 onward.

3.12 is the version both platforms can actually agree on:

- **Floor:** `scanpy` and `numba` both declare `requires_python >= 3.10`.
- **Ceiling:** `tensorflow` 2.17 (pulled in by `panhumanpy`) publishes no wheels
  above `cp312`. So 3.13 would silently rule out our Azimuth annotation route.
- Every compiled dependency has a `cp312` wheel for Windows and Linux alike.

### Setup

`.venv/` is local and git-ignored — never commit it.

```bash
# Linux
python3.12 -m venv .venv
source .venv/bin/activate
```

```powershell
# Windows
py -3.12 -m venv .venv
.\.venv\Scripts\Activate.ps1
```

Then, on both:

```bash
python -m pip install --upgrade pip setuptools wheel
python -m pip install -r requirements.txt
python -m ipykernel install --user --name=nanoplastics-scrnaseq \
  --display-name="Python (nanoplastics-scrnaseq)"
```

Install `requirements.txt` in **one** `pip install` pass rather than adding packages
incrementally, so the resolver sees every constraint at once (see numpy, below).

**A C++ compiler is required on both platforms.** `scanorama` depends on `annoy`,
which publishes no wheels for Windows or Linux and always builds from source:
MSVC Build Tools (C++ x86/x64 workload) on Windows, `gcc`/`g++` on Linux. If the
build fails on Windows, retry from the "x64 Native Tools Command Prompt for VS 2022"
so `vcvarsall` is on PATH — do not swap batch-correction method to work around it.

### Requirements files

`requirements.txt` is now a **direct-dependency list**, not a `pip freeze`. A freeze
cannot be shared across platforms: on Windows the graph includes `pywinpty` and
`tensorflow-intel`, neither of which installs on Linux. The direct list resolves
correctly on both.

| File | Purpose |
|---|---|
| `requirements.txt` | ~20 direct dependencies — what both of us install from |
| `requirements-lock-linux-py310.txt` | the original author's Linux/3.10 freeze, kept for reference |
| `requirements-lock-windows-py312.txt` | full freeze from this machine, for exact reproduction |

The three previously-outstanding packages — `scanorama` (Section 8), `panhumanpy`
and `celltypist` (Section 11) — are now **in** `requirements.txt`. `panhumanpy`
downloads its model weights on first use, not via pip. `harmonypy` was unused and
has been dropped.

### Consequence: numpy is 1.26.4, not 2.2.6

`panhumanpy` pins `tensorflow==2.17` and `scikit-learn==1.6.0` exactly, and TF 2.17
requires `numpy>=1.26.0,<2.0.0`. The old Linux freeze pinned `numpy==2.2.6`, so
adding these packages rolls numpy back across the 1.x/2.x boundary. Verified that
nothing else objects — every other constraint is a lower bound or a wide range:

| Package | Constraint | OK at numpy 1.26 / sklearn 1.6 |
|---|---|---|
| scanpy 1.11.5 | `numpy>=1.24.1`, `scikit-learn>=1.1.3` | ✓ |
| anndata 0.11.4 | `numpy>=1.23` | ✓ |
| scipy 1.15.3 | `numpy>=1.23.5,<2.5` | ✓ |
| numba 0.67.0 | `numpy>=1.22,<2.6` | ✓ |
| celltypist | `scikit-learn>=0.24.1` | ✓ |

1.26.4 is the newest release under TF's ceiling, so this is the resolver's best
answer rather than a fallback. This is the re-freeze this plan always called for.

### Data files

The `.h5ad` files are git-ignored — at 167–272 MB each they exceed GitHub's
100 MB file limit, so they must be obtained separately.

## Notebook structure

`nanoplastics_pipeline.ipynb` at the repo root, using the `nanoplastics-scrnaseq`
kernel.

| # | Section | Tutorial cells | Deviation from tutorial |
|---|---|---|---|
| 1 | Load 4 `.h5ad` files from `layers['counts']`; add `obs['Sample']` per file | 0-3 | Tutorial loads raw CSV; we load the h5ad counts layer |
| 2 | Doublet detection per sample, pre-merge | 4-22 (scVI/SOLO) | **Scrublet** (`sc.pp.scrublet`) instead. Its automatic threshold is unusable here — see "Section 2 outcome" below — so we threshold at fixed detection sensitivity |
| 3 | QC, **split into 3a and 3b**. 3a: compute mito% (`MT-`) and ribo% (KEGG_RIBOSOME — the tutorial's Broad URL 301-redirects to `gsea-msigdb.org` but resolves when redirects are followed, which `pd.read_table` does by default), plot all four distributions, filter **nothing**. 3b: apply per-sample MAD thresholds | 23-45 | Thresholds picked **per sample** from its own distribution. The 3a/3b split is forced by that requirement — a per-sample threshold cannot be written down before four distributions have been looked at. **Ribo% is computed and plotted but never filtered on**; the tutorial's `ribo < 2` is inert on its own data and catastrophic on ours — see "Section 3 outcome" |
| 4 | Concatenate the 4 samples after per-sample doublet+QC, then `sc.pp.filter_genes(min_cells=10)` on the merged object | 70-74, 85 | Same concat pattern, no CSV loop needed. **Gene filtering lives here**, a step `PLAN.md` previously assigned to no section: after the join a gene needs 10 cells across all 28,436 rather than 10 within every sample. It must also precede Section 5, since `normalize_total`'s divisor is a sum over the genes present |
| 5 | Normalization: `normalize_total(target_sum=1e4)` + `log1p`, stored explicitly as `layers['normlog']` | 46-52, 88 | Tutorial uses `.raw`; `CLAUDE.md` requires keeping both `counts` and `normlog` layers. `target_sum=1e4` is **not a free choice** — CellTypist (Section 11) documents log1p-normalized to 10,000 counts/cell as its required input, so scanpy's median-based default would quietly break Task 3 |
| 6 | HVG selection: `seurat_v3` flavor, `n_top_genes=3000`, `layer='counts'`, `batch_key='Sample'`, **`subset=False`** | 91 (commented-out alt) | The tutorial's commented-out multi-batch line already matches `CLAUDE.md`'s preferred route — use it over the tutorial's dispersion-based default. **`subset=False` deviates from the tutorial** (its cell 57 discards non-HVGs and relies on `.raw`): Sections 10, 13 and 16 need the full 20,352-gene matrix for marker genes, pseudobulk DE and Hallmark scoring. `sc.tl.pca` reads `var['highly_variable']` automatically, so flagging suffices |
| 7 | Scale + PCA. Scaling runs on a **separate HVG-only copy**; `.X` on the main object stays `normlog` for good. **No `regress_out`** | 57-61 | Two deviations: tutorial regresses out `total_counts`/`pct_mt`/`pct_ribo` and scales in place — see "Section 7 decisions" |
| 8 | Batch correction: **scanorama** (`sigma=0.1`, `knn`/`alpha` left at defaults) → `X_scanorama`; both sides of the tradeoff measured | 92-98 (scVI) | scanorama instead of scVI. `sigma` chosen from a measured sweep rather than the default — see "Section 8 outcome" |
| 9 | Neighbors (`use_rep='X_scanorama'`), UMAP, Leiden — resolution justified by stability | 62-68, 97-99 | Same mechanics, different embedding input |
| 10 | Marker genes: `rank_genes_groups`, manual cluster→cell-type dict + dot plot | 102-112 | Matches the tutorial directly |
| 11 | **Automated annotation, run by us**: `panhumanpy` + `celltypist` (`majority_voting` fed our Leiden clusters). Reconcile both against marker labels per cluster | — (new) | The Azimuth half of Task 3, run by us in Python rather than inherited or via R |
| 11b | **Verification**: cross-tabulate final labels against `predicted.celltype` and `CoDi` | — (new) | External labels are a check on our result, never the result itself |
| 12 | Composition analysis: cell-type proportions across the 4 samples | 119-127 | Generalized from 2 conditions to 4 samples |
| 13 | DE: pseudobulk (counts summed per cell-type × sample) via `pydeseq2`, each exposed sample vs. control | 129-158 | Pseudobulk DESeq2 instead of per-cell diffxpy/scVI, per `CLAUDE.md`'s pseudoreplication concern |
| 14 | Pathway enrichment: GSEA primary (`gseapy`, MSigDB Hallmark/KEGG/Reactome), ORA secondary | 159-165 | Tutorial uses ORA (`gp.enrichr`) only; `CLAUDE.md` prefers GSEA |
| 15 | Size-specific effects: unique-to-40nm / unique-to-200nm / shared / mixture-only | — (new, Task 6) | Not in tutorial — nanoplastics-specific |
| 16 | **Per-cell signature scoring** (`sc.tl.score_genes`, Hallmark `TNFA_SIGNALING_VIA_NFKB`): compare score *distributions* across samples and cell types; UMAP colored by score | 172-178 | Same mechanic as the tutorial's `datp_sig` scoring, different gene set — **promoted to Task 6 core**, see below. Skip the tutorial's `mannwhitneyu` (cell 177) |

Sections 1-9 → Tasks 1-2; 10-11b → Task 3; 12 → Task 4; 13-14 → Task 5;
15-16 → Task 6.

### Why Section 16 is core, not bonus

`CLAUDE.md` lists signature scoring as bonus candidate #6, yet simultaneously
describes it as complementing Task 5/6 "without betting on one gene" — which is
core Task 6 reasoning. It is **not redundant** with Section 14's GSEA, because
the two operate at different units of analysis:

- **GSEA** runs on a ranked gene list from a **pseudobulk** contrast → one
  enrichment statistic per (cell type × contrast). Every cell in a group is
  collapsed to a single number *before* testing.
- **`score_genes`** assigns a score to **every individual cell**, enabling four
  things GSEA cannot do even in principle:
  1. Reveal **within-cell-type heterogeneity** — particle uptake is stochastic,
     so a bimodal score within monocytes would be direct evidence of
     uptake-dependent response. Pseudobulk averages that away entirely.
  2. UMAP colored by score, showing where the response localizes.
  3. Compare all cell types on one common scale — a direct test of `CLAUDE.md`'s
     myeloid-not-lymphocyte expectation.
  4. Compare all 4 samples at once, rather than as three pairwise contrasts.

The gene sets aren't equivalent either: GO Biological Process has NF-κB-adjacent
terms, but they're redundant and noisy, whereas Hallmark
`TNFA_SIGNALING_VIA_NFKB` is a curated ~200-gene set built for this readout.

**Scope limit:** the *significance testing* genuinely is duplicated — scoring
adds no new p-value beyond GSEA. Its value is **descriptive and visual only**.
Do not port the tutorial's cell 177 (`mannwhitneyu` across per-cell scores): with
one donor, cells are not independent replicates, so that test hits exactly the
pseudoreplication problem `CLAUDE.md` flags for DE. Record the caveat in a
markdown cell instead of reporting a p-value.

The van den Brink dissociation-stress signature (the other half of `CLAUDE.md`'s
bonus #6) stays deferred — `CLAUDE.md` ties it to the mito%/complexity QC rework,
which needs a completed Task 3 first.

## Progress

Notebook is `nanoplastics_pipeline.ipynb`, kernel `nanoplastics-scrnaseq`.
Work happens on branch `windows-env-setup` (not yet merged to `main`).

| Section | Task | Status |
|---|---|---|
| Environment | — | done (Python 3.12.10, `.venv`, kernel registered) |
| 1 — load 4 samples | 1 | **done, verified** |
| 2 — doublet detection | 1 | **done, verified** |
| 3a — QC metrics + distributions | 1 | **done, verified** |
| 3b — apply QC thresholds | 1 | **done, verified** |
| 4 — concat + gene filter | 1-2 | **done, verified** |
| 5 — normalization | 1 | **done, verified** |
| 6 — HVG selection | 2 | **done, verified** |
| 7 — scale + PCA | 2 | **done, verified** |
| 8 — batch correction (scanorama) | 2 | **done, verified** |
| 9 — neighbours, UMAP, Leiden | 2 | **done, verified** |
| 10 — marker genes + cluster annotation | 3 | **stage 1 done**; stage 2 (the cluster→cell-type map) next |
| 11-16 | 3-6 | not started |

Cells: 34,078 loaded → **31,839** after doublet removal → **28,436** after QC (10.7%
removed, spread 1.6 pp across samples).
Genes: 24,380 union → **20,388** after the inner join → **20,352** after `min_cells=10`.

Sections 1-7 have been run **top-to-bottom on a restarted kernel** with no errors,
satisfying the end-to-end requirement in the verification checklist (Sections 1-5 at
execution counts 1-24; Sections 6-7 at 25-30; Section 8 at 31-34; Section 9 at 35-40).
A write-only checkpoint is saved at `results/combined_after_qc.h5ad` (~1 GB, git-ignored).
Nothing auto-loads it — a checkpoint that reloads itself is the same bug as a stale
output. Restore by hand with `sc.read_h5ad` if a kernel dies.

### Section 1 outcome

Loaded from `layers['counts']`; Seurat's `X`, `obsm`, `varm`, `uns` and clustering
discarded. `obs` reduced to `predicted.celltype` + `predicted.celltype.score` (kept as
the Section 11b benchmark) plus our `Sample` label. Barcodes prefixed per sample
(`40nm_AAACCC...`) because all four files reuse the `-1` suffix and would collide on
concat. Counts downcast to `float32` — lossless, as `float32` is exact for integers
below 2²⁴ and the largest count is 14,292.

**Open issue for Section 4:** the four files have **different gene sets** (22,613 /
23,206 / 21,715 / 21,961), because each was gene-filtered independently upstream.
`sc.concat` defaults to `join='inner'`, keeping only genes present in all four. That is
the correct choice — `join='outer'` would fill absent genes with **zeros**, asserting
"not expressed" where the truth is "filtered out", and would fabricate large DE hits in
exactly the size comparison Task 6 depends on. **Resolved at Section 4:** 20,388 of a
24,380 union survive (16.4% lost); `mixture` and `control` cost ~525 genes each, `40nm`
209, `200nm` 117.

**Correction to this note:** it originally said to check `TNF`, `IL6` and **`CXCL8`**
survive. That check false-alarms — `CXCL8` is absent from all four files, not because
the join lost it but because this data uses pre-2017 gene symbols and calls it **`IL8`**.
See "Legacy gene symbols" under Known limitations. All 28 cytokine and cell-type marker
genes checked at Section 4 survive both the join and `min_cells=10`.

### Section 2 outcome

Scrublet's automatic threshold **failed** on this data. It assumes a bimodal score
distribution; ours is unimodal, because only *heterotypic* doublets are detectable
(two same-type cells fused look like one brighter cell of that type). So
`threshold_minimum` drifted into the simulated-doublet tail: on `40nm` it chose 0.705
against a maximum observed score of 0.637, calling **zero** cells, and across samples it
ranged 0.22–0.70.

**Rule adopted:** threshold at fixed detection sensitivity. The simulated doublets are
known positives, so the fraction of them above a cutoff *is* that cutoff's sensitivity;
taking the median simulated score pins it near 50% in every sample.

A shoulder-based rule (local minimum of the observed density) was tried and **rejected**:
`doublet_score` is a fraction of k neighbours and therefore discrete, so the "shoulder"
tracks how many cells populate the discrete values — giving the largest sample the
weakest detection (42.7% vs 49.4%). That bias would distort the very Task 4 composition
comparison it is meant to protect.

Removed **2,239 of 34,078 cells (6.6%)**: 618 / 898 / 358 / 365.

**Finding to revisit in Section 12:** sensitivity-corrected doublet burden is a flat
~11–14% across all four samples, rather than scaling with cell count as the 10x loading
model predicts (4.9%–10.1%). Most likely the channels were loaded comparably and differ
in *recovery*. Treat ~11–14% as an upper bound — simulated doublets are cleaner than
real ones, which inflates measured sensitivity. Once cell-type labels exist, test
directly by computing each sample's homotypic fraction (Σpᵢ² over cell-type
proportions) and checking whether any sample is unusually skewed.

**Dependency worth knowing:** `sc.pp.scrublet` needs `scikit-image`, which is *not* a
hard scanpy dependency. Without it scanpy imports fine and the call fails at runtime.
Now pinned in `requirements.txt`.

### Section 3 outcome

**The rule adopted:** MAD-based outlier detection computed within each sample.
`log1p_total_counts` and `log1p_n_genes_by_counts` at 5 MADs both ends;
`pct_counts_mt` at 3 MADs **upper only**, since low mito% is healthy and a two-sided
rule would flag the best cells. Depth and gene counts use the log scale because both
are log-normal — a symmetric rule on raw counts flags far more cells above the median
than below, purely from skew.

Removed **3,403 of 31,839 cells (10.7%)**, spread 10.0–11.6% across samples. The mito
filter dominates (6.5–8.2%); the low-depth end removes only 1.3–3.2%.

**Per-sample thresholds are justified, not ceremonial.** Median depth differs 1.77×
across samples (4,334–7,674), genes detected 1.47×, mito% 1.56×. One shared threshold
would cut deeper into some samples than others and manufacture the composition
difference Task 4 is meant to measure.

**Two readings of the violin plots were wrong and were overturned by cell profiles.**
Worth recording because it generalises: a violin shows *where* cells are, never *what*
they are.
- The `total_counts` lower bound looked like it sat in a populated region. It removes
  1.1–2.7% per sample, and those cells have a genes-per-count ratio of 0.6–0.7 against
  0.35 in kept cells — nearly every UMI hitting a different gene, the signature of
  sampling ambient RNA rather than a transcriptome.
- `pct_counts_mt` at 3 MADs looked like it cut through the healthy bulk. The cells it
  removes carry ~11% ribosomal content against ~32% in kept cells — a distinct
  compromised population, not healthy cells at the edge.

Ribosomal content, which we deliberately do **not** filter on, was the axis that
separated those populations. Two metrics that overlap completely on one axis can
separate cleanly on another.

**Why the tutorial's thresholds were not inherited.** Traced through its own saved
outputs: `pct_counts_mt < 20` and `pct_counts_ribo < 2` together removed **17 cells
(0.3%)** of its 5,526. Over 99.7% of its cells already sat below 2% ribo, a profile
consistent with single-nucleus data. Ours is whole-cell PBMC with ribo% at 23–30%
median, where the same line would delete ~97% of the dataset. This is `CLAUDE.md`'s
named pitfall — *a threshold present in code but numerically inert* — found in the
source we were adapting from, and the reason every filter here prints before/after
counts and asserts it removed something.

**Feeds bonus #1:** the mito filter removes slightly more from exposed samples
(7.2–8.2%) than from control (6.5%) — small, but consistent in direction across all
three exposures. Technical damage, or genuine exposure-induced mitochondrial stress?
Not resolvable in this pass.

### Sections 4-5 outcome

Concatenated with `join='inner', merge='same'` — `merge='same'` is required or
`sc.concat` drops `.var` entirely and the `mt`/`ribo` flags vanish silently, leaving
both percentages a uniform 0.0 that reads as unusually clean data. Verified 13 MT-
genes and 85 ribo genes survived.

`min_cells=10` removed only **36 genes (0.2%)** — nearly inert, because the upstream
processing already gene-filtered each file hard (which is why the gene sets differ at
all). Measured alternatives: 3 removes 0, 30 removes 1,796, 100 removes 5,019. Chose 10
because the binding constraint is **rare-cell-type annotation** in Task 3, not noise
removal: pDCs are ~0.3% of PBMCs (~85 cells here) and `min_cells=100` would strip their
markers and make that population un-annotatable.

Normalization to `target_sum=1e4` then `log1p`. Layer contract now fixed for the rest of
the pipeline: `layers['counts']` raw integers (max 4,212, unchanged from Section 4),
`layers['normlog']` log1p CP10K unscaled (max 8.47), `.X` equal to `normlog` — and, per the
Section 7 decision below, it stays that way for the rest of the pipeline: scaled values never
enter the main object at all.

### Section 6 outcome

3,000 HVGs flagged of 20,352 genes; `subset=False`, so nothing was discarded.

**The mean–variance correction demonstrably worked:** only **1 of 85** ribosomal genes
and **0 of 13** mitochondrial genes were selected. Those are among the most highly
expressed genes in every cell, so a naive variance ranking would have put them at the
top — excluding them is the whole purpose of `seurat_v3`'s loess fit.

**1,711 genes are variable in all four samples** (`highly_variable_nbatches == 4`). That
gives `n_top_genes` a data-driven floor it otherwise lacks: 3,000 reaches ~1,300 genes
past the reproducible core. This is the softest parameter in the pipeline — unlike
`target_sum=1e4` (forced by CellTypist) or the QC thresholds (derived from the data),
it is conventional. The tutorial's stated heuristic is "3k is fine for 10k cells"; we
have 28,436, so 3,000 is if anything conservative. If a stability check is ever wanted,
the test is whether clustering changes materially at 2,000 vs 3,000 — the same standard
Section 9 applies to Leiden resolution.

**9 of 12 cell-type markers were flagged.** The three missing — `CD3E`, `CD8A`,
`FCGR3A` — are expected rather than a problem: `CD3E` is expressed in every T cell, so
it has a high mean and unremarkable standardized variance. **A good marker is not the
same thing as a variable gene.** HVGs drive the embedding; Section 10 reads all 20,352
genes for markers, which is exactly why `subset=False` matters.

**New dependency:** `flavor='seurat_v3'` requires `scikit-misc` (0.5.2) for its loess
fit. Same failure shape as `scikit-image` for Scrublet — not a hard scanpy dependency,
so scanpy imports fine and the call dies at runtime. Now pinned in `requirements.txt`;
it needs `numpy>=1.26.4`, exactly the version TensorFlow's ceiling pins us to.

### Section 7 outcome

Ran clean in the notebook at execution counts 26-30. `X_pca` is (28,436 x 50), built from
the 3,000 HVGs scaled on a throwaway copy with `max_value=10`; `X_pca_sel` is (28,436 x 30)
and is what Section 8 integrates.

**The layer contract held, and the check is the point.** `.X` prints `min 0.00, max 8.47` —
identical to Section 5's recorded value, and the `min 0.00` is the proof: scaled data
would carry negatives (the unclipped floor was -4.32). `sc.pp.scale` is never called on
`adata` at all, so the tutorial's stash-then-scale re-run trap has nothing to overwrite.

**50 PCs capture 22.5% of total variance, and that is a good number, not a poor one.**
After scaling, total variance is exactly 3,000 units (one per gene). Under pure noise,
50 of 3,000 components would capture 50/3000 = **1.7%**; 22.5% is **13x chance**. PC1
alone holds 8.85% = **265 gene-equivalents** of variance compressed into one axis.

**Two percentages, two denominators — do not conflate them in the report.** 30 PCs capture
**92.1% of the variance in `X_pca`**, which is **20.7% of the total**. Neither number is a
quality score: a pure-noise dataset would still yield a high cumulative-of-captured
figure. The 13x-chance comparison above is the one that says PCA found structure.

**`N_PCS = 30`, chosen from the `drop` column rather than the cumulative curve.** `drop`
(the ratio between consecutive components) falls to ~1.01 by PC20 and stays there —
components within ~1% of each other in size, the signature of noise dimensions. But `drop`
tests distinguishability in *size*, which is a poor rare-population test: pDCs at ~0.3%
(~85 cells) have a real coordinated program contributing almost no total variance, and
Task 3 has to resolve them, as Tasks 5-6 need classical vs non-classical monocytes. So 30
rather than 20 (margin for rare populations) and rather than 40 (ten indistinguishable
components handed to scanorama's mutual-nearest-neighbour matching, in a dataset whose
depth gradient is already confounded with exposure).

**Reversible by design.** `X_pca` retains all 50 components and the selection lives in a
separate key, so revising the cut is one integer and two cells — no PCA recomputation, and
the elbow plot the choice was justified from stays the plot the pipeline produced.

**Expectation to carry into Section 9.** The largest `drop` is not at PC1 but at PC3->4
(2.75, against 1.29 at PC2->3). The first three components are most likely the major
lineage splits (myeloid / T / B-NK), with everything after PC3 substructure within them.
If Leiden does not recover roughly that hierarchy, something upstream is wrong.

**A numeric claim in this file was wrong and has been corrected:** an earlier draft said
PC50 had fallen below one gene's worth of variance. It holds 2.4 gene-equivalents. The
0.79 figure was the average across the whole 2,950-component tail, not PC50 itself.
`drop`, not gene-equivalents, is what marks the boundary.

### Section 7 decisions

Three things PLAN.md previously left open, now settled.

**1. No `regress_out` — dropped entirely, not merely "optional".**

`sc.pp.regress_out` fits, per gene, a linear model of expression against the listed
covariates and keeps only the residuals — deleting whatever variation those covariates
can explain. The tutorial (cell 58) regresses `total_counts`, `pct_counts_mt` and
`pct_counts_ribo`. All three are refused here:

- **`total_counts` is confounded with condition.** Median depth runs control 7,736 >
  mixture 6,276 > 40nm 5,436 > 200nm 4,390 (Section 4 output) — a 1.77x gradient that
  tracks exposure. So "the best-fit line through depth" and "the best-fit line through
  exposure" are nearly the same line, and subtracting one subtracts the other. A real
  exposure effect — say `LYZ` at 4.0 in control monocytes and 3.2 in 200nm monocytes —
  would be re-attributed to depth and flattened toward a common value, silently. This is
  `CLAUDE.md`'s over-correction hazard, arriving one section before Section 8 and through
  a different door. **With one donor it is unresolvable**, so the safe direction is to
  leave the variance in and let Section 8 handle batch structure explicitly.
- **`total_counts` is also largely redundant.** Section 5's `normalize_total` already
  divided every cell by its own depth. Regression is a second, blunter pass at a problem
  already addressed.
- **`pct_counts_ribo`** — the identical argument that stopped us *filtering* on ribo in
  Section 3. Ribosomal content tracks translational activity, i.e. activation; regressing
  it removes activation-associated variance, which is the signal Task 6 is looking for.
- **`pct_counts_mt`** — Section 3b already removed mito outliers at 3 MADs, so little is
  left to correct, and the same "is it damage or is it exposure biology?" ambiguity
  flagged there (bonus #1) applies.

Regressing out covariates that correlate with the biology is discouraged generally;
here it is specifically disqualifying.

**2. Scale an HVG-only copy — `.X` on the main object is never scaled.**

`sc.tl.pca` reads `var['highly_variable']` and ignores the other 17,352 genes, so scaling
the full matrix z-scores 85% of it purely to discard the result. Measured cost of doing
it anyway: the 3,000-gene dense block is **682 MB** — `sc.pp.scale` returns **float64**,
not the float32 the earlier 2.3 GB estimate in this file assumed, so the full matrix
would be ~4.6 GB, not 2.3 GB. On a 21 GB machine with ~10 GB free that is survivable but
pointless, and it would land in every subsequent checkpoint.

The larger gain is not memory. Scaling a separate object means `.X` on `adata` stays
log-normalized permanently, and the mislabeled-layer bug class that `CLAUDE.md` names as
a pitfall — Sections 10, 13 and 16 reading z-scores while believing they read expression
— stops being reachable rather than merely being avoided by care. PCA coordinates are
copied back to `adata.obsm['X_pca']`; the scaled copy is not retained.

**3. `max_value=10` — kept, and now justified from the data rather than by convention.**

Scaling the HVG block without clipping gives a z range of **-4.32 to 108.54**. The
asymmetry is the point: z-scores have a hard floor (a cell cannot express less than zero,
so nothing sits more than 4.32 SD below its gene's mean) and no ceiling at all. Clipping
can therefore only ever touch the top tail.

| cutoff | values above | % of all | genes with >=1 |
|---|---|---|---|
| 5 | 473,500 | 0.555% | 2,774/3000 |
| 10 | 116,207 | 0.136% | 2,081/3000 |
| 20 | 24,539 | 0.029% | 1,444/3000 |
| 50 | 897 | 0.001% | 409/3000 |

So `max_value=10` is **not** numerically inert — `CLAUDE.md`'s named pitfall, checked for
rather than assumed. It alters 116,207 values, ~4 per cell, in 90.6% of cells.

What it clips is the right thing. The most extreme genes are `BPI` (108.5), `APOD`
(101.2), `ANXA3` (100.5), `SUCNR1` (89.5) — near-binary genes expressed in 10-25 cells
each. Their enormous z-scores come from a tiny SD in the denominator, not from a large
biological difference. Unclipped, a single such cell carries `(108/10)^2` = **118x** the
PCA leverage of a cell at z=10, letting a dozen cells bend a principal component. Clipping
leaves those cells still maximal for the gene, so rare populations remain findable in
Section 10; it only removes their ability to dominate the embedding.

5 was rejected: it clips 4x as many values across 92% of genes, and z=5-10 is an ordinary
range for a genuine marker in a genuine cell type.

**4. PC count.** Compute 50, plot the variance-ratio elbow, then explicitly keep the
chosen number. The truncation must happen here, not at `sc.pp.neighbors` as the tutorial
does it (cell 62, `n_pcs=30`): `sce.pp.scanorama_integrate` reads its input basis
wholesale and returns `X_scanorama` at matching dimensionality, so any dimensions still
present at Section 8 are dimensions that get integrated.

*Where 50 came from, stated plainly:* it is **scanpy's own default**
(`sc.settings.N_PCS = 50`, applied whenever `n_comps` is omitted, as the tutorial does at
cell 60). It entered this pipeline by convention, not from the data. Writing it
explicitly changes nothing about the computation — it only makes the number visible
rather than implied.

That is acceptable **because `N_PCS_COMPUTED` is a ceiling, not a threshold**, and the two
carry different burdens of proof. A ceiling has one job: be comfortably larger than the
elbow, so the elbow is visible inside the computed range. Too high costs a little compute;
too low is the real failure, because the truncation would happen before the structure ends
and nothing in the output would reveal it. **Verified after the fact:** PC50's variance
ratio is 0.00079 against PC1's 0.08848 (112x smaller) and its `drop` is ~1.01, so
consecutive components are indistinguishable in size from one another. The elbow is far
inside the range. (In gene-equivalents — total variance is exactly 3,000 units after
scaling, one per gene — PC1 holds 265 and PC50 still holds 2.4; the tail *as a whole*
averages 0.79 per component, but no individual PC has crossed below one gene's worth by
PC50. `drop`, not gene-equivalents, is what marks the boundary.) Had `drop` still been >1.1 at PC50 with the
cumulative curve climbing steeply, the correct response would have been to recompute at
100 rather than to truncate.

`N_PCS` — the number actually carried into Section 8 — is a **threshold** and gets no such
latitude: it decides what reaches integration, clustering and annotation, so it is chosen
from the elbow output and its reasoning recorded alongside it. Section 3's inherited
`ribo < 2` is the cautionary case of a threshold accepted as though it were a harmless
default.

This is also why the section computes 50 and slices to `X_pca_sel` rather than calling
`sc.tl.pca(n_comps=N_PCS)` directly. Recomputing PCA to change the cut would leave the
elbow plot and the selection inconsistent — the plot the choice was justified from would
no longer be the plot the pipeline produced.

**Hazard to avoid when writing Section 7.** The tutorial's pattern is
`adata.layers['normlog'] = adata.X.copy()` immediately before `sc.pp.scale`. Do **not**
copy that — Section 5 already set `normlog`, and the line is a mutate-what-you-read trap:
run the cell a second time and `.X` is already scaled, so `normlog` gets silently
overwritten with z-scores. Nothing errors; Sections 13 and 16 would then read scaled data,
which `CLAUDE.md` explicitly forbids for DE and `score_genes`. Same shape as the Section 3b
re-run hazard. Decision 2 above removes the trap at its root — `sc.pp.scale` is never called
on `adata`, so there is nothing for a re-run to overwrite.

### Section 8 outcome

Ran clean at execution counts 31-34. `adata.obsm['X_scanorama']` is (28,436 x 30);
`X_pca_sel` is left untouched so Section 9 can cluster the uncorrected embedding as a
comparison.

**The batch effect was mild before we touched it — the finding that shaped the section.**
Uncorrected mixing is 0.480, which is **66% of the 0.727 ceiling**. Cells already sat
beside the right cell types largely regardless of channel. Correction here is a touch-up,
not a rescue, and that is why the gentle end of the range won.

**Two metrics, both built only from our own data.** `CLAUDE.md` names the
correction-vs-conservation axis, but a UMAP cannot separate "well mixed" from
"over-corrected" — both look mixed.

- **Mixing:** mean fraction of a cell's 15 nearest neighbours drawn from a *different*
  sample. The ceiling is 1 − Σp² = **0.727**, not 1.0, because neighbours include a cell's
  own sample at its natural share. Scoring against 0.727 avoids rating the largest sample
  worst-mixed purely for being largest — the size bias that got a candidate Scrublet
  threshold rejected in Section 2.
- **Containment:** of each cell's 15 nearest **own-sample** neighbours before correction,
  what share remain within its 100 nearest own-sample neighbours after. Restricting to one
  sample is what makes it a fair conservation test **with no cell-type labels**: inside a
  single sample there is no between-sample difference for correction to legitimately fix.

**The inherited `predicted.celltype` labels were deliberately not used**, though they
would have been a convenient yardstick. Tuning integration on them would let them shape
the embedding, then the clusters, then our annotation — leaving Section 11b's "independent
check" comparing our result against labels that helped produce it. Keeping them out here
is what keeps that check meaningful.

**A metric was designed, tested and thrown away.** The first conservation attempt was
exact 15-NN set overlap. Calibration killed it: jitter at 10% of the neighbourhood radius —
far too small to reorganise anything — already halved it, because in dense regions many
cells sit near-equidistant, so *which* 15 are closest is unstable. Containment asks the
weaker, correct question. **Calibrated against known-benign changes** (notebook cell 47):

| transformation | containment |
|---|---|
| identity | 1.000 |
| jitter at 5% of neighbourhood radius | 0.989 |
| jitter at 10% | 0.872 |
| jitter at 25% | 0.468 |
| drop PCs 21-30 (20 dims) | 0.978 |
| random reshuffle (floor) | 0.014 |

Recorded because the near-miss generalises: `retention = 0.218` for scanorama's defaults
looked like proof of over-correction and would have produced a confident, well-argued
paragraph defending a conservative setting **on the basis of an artifact**. A number with
no reference scale is not evidence.

*(The `drop PCs 21-30` row at 0.978 also re-tests Section 7 from a different angle: those
components carry almost no structure, exactly as the `drop → 1.01` reading predicted.)*

**Of the three parameters, only one is real.**

- **`knn` — not a lever.** At 10 / 20 / 40, mixing came out 0.640 / 0.640 / 0.636. A
  four-fold change moves the result by 0.004. Default 20 kept.
- **`alpha` — inert.** Scanorama's own printed alignment scores bottom out at **0.682**,
  so no `alpha` below that changes which pairs merge, and anything above it starts
  refusing genuine overlaps. A "conservative alpha" would be a parameter that looks
  principled and does nothing — `CLAUDE.md`'s named pitfall, the same shape as the
  tutorial's inert `ribo < 2`. Default 0.1 kept.
- **`sigma` — the only real knob**, and it does not behave as expected.

**The full sigma sweep (nine values; the notebook re-runs three).**

| sigma | mixing | of ceiling | containment |
|---|---|---|---|
| uncorrected | 0.480 | 66% | 1.000 |
| 0.05 | 0.619 | 85% | 0.983 |
| **0.1 — chosen** | **0.625** | **86%** | **0.978** |
| 0.25 | 0.627 | 86% | 0.953 |
| 0.5 | 0.620 | 85% | 0.884 |
| 1 | 0.613 | 84% | 0.770 |
| 2 | 0.624 | 86% | 0.682 |
| 5 | 0.638 | 88% | 0.590 |
| 15 — default | 0.640 | 88% | 0.547 |
| 50 | 0.637 | 88% | 0.554 |

**It is barely a tradeoff at all, which is the opposite of the expectation going in.**
Mixing spans only 84-88% across the entire range; containment spans 0.98 to 0.55. Small
`sigma` is not a cautious sacrifice of mixing for safety — it delivers *the same* mixing
with far less disturbance. The default `sigma=15` costs half the local structure for two
percentage points of mixing.

`sigma=0.1` is a **located** value, not a lucky one: above `sigma≈5` the parameter
saturates entirely (5, 15 and 50 all give 88% / ~0.55) and below `sigma≈0.25` containment
plateaus near 0.98. It sits on the flat top of both curves.

*Interpretation, offered as such:* `sigma` sets how broadly each correction is smoothed
across cells, and large values exist for batch effects that hit different cell types
differently. This batch effect is mild, and Section 3 characterised it as largely a
uniform per-sample depth offset. A small, tightly targeted correction is the right shape
of fix; broad smoothing mostly moves cells that did not need moving.

**Final result:** mixing 66% → **86%** of ceiling, containment **0.978** — better than the
5%-jitter reference. `CLAUDE.md`'s requirement (well mixed without erasing structure) is
met as two numbers rather than an impression. Section 9 supplies the visual half.

**Mechanism worth knowing:** the alignment-score matrix printed identically on all three
runs. Matching and correcting are separate stages in scanorama — `knn`/`alpha` decide
which populations correspond, `sigma` decides how far they move. That is why one parameter
can be wholly inert while another controls everything.

**Biological by-product:** see "Empirical findings carried forward" near the top of this
file — the alignment scores rank 200nm as the least similar sample to control.

### Section 9 outcome

Ran clean at execution counts 35-40. `adata.obs['leiden']` holds 16 clusters at
`resolution=0.4`, built on a neighbour graph over `X_scanorama` with `n_neighbors=15`.

**`use_rep='X_scanorama'` is asserted, not assumed.** Omitted, `sc.pp.neighbors` falls back
to `adata.obsm['X_pca']` — 50 *uncorrected* components — discarding Sections 7 and 8 while
running perfectly and plotting beautifully. The cell reads back
`adata.uns['neighbors']['params']['use_rep']` rather than trusting the argument was typed.

**One of the two resolution methods failed, and this is recorded rather than glossed.**

- **Plateau test: no signal.** Cluster count rises almost linearly across
  0.1-1.5 (8, 12, 14, 16, 19, 20, 23, 26, 28, 32) with no flat span. By the standard set out
  in the section itself, a straight diagonal privileges no resolution. It contributed
  nothing, so the choice rests on **one** line of evidence rather than two agreeing ones.
- **Subsampling stability: decisive.** Mean ARI against a re-clustered 90% subsample peaks
  at **0.3 (0.913) and 0.4 (0.912)**, falls sharply to 0.786 at 0.5, and declines to 0.68 by
  1.5.

**`resolution=0.4` chosen.** 0.3 and 0.4 are tied within the noise of three repeats, so the
tie broke on what later tasks need: 16 clusters rather than 14 at no measurable stability
cost; Task 3 must split classical (`CD14`) from non-classical (`FCGR3A`) monocytes, which is
where `CLAUDE.md` expects the effect; and over-splitting is recoverable in Section 10 by
reading markers while under-splitting is invisible.

**Bias worth recording:** stability structurally favours fewer, larger clusters, because
those reproduce more reliably under subsampling. Taken as `argmax` it would have argued for
`resolution=0.1` and the loss of every rare population — including the pDCs that Section 4's
`min_cells=10` was chosen to protect. A sound metric can still be systematically biased
toward one kind of answer, which is why the sweep also tracked the smallest cluster and the
count under 50 cells.

**Cluster sizes** run 6,849 (24.1%) down to 5 (0.02%), with four clusters under 300 cells
(7, 9, 12, 14) — rare populations survived, which was the point.

**Integration confirmed a second time, independently of Section 8's metrics.** Sample
enrichment per cluster sits between 0.85 and 1.35 almost everywhere. Under-corrected batch
would have produced sample-specific clusters and values near 0 and 4. The exceptions are
recorded under "Empirical findings carried forward".

**A diagnostic of ours was mis-specified and has been corrected.** The composition cell
originally flagged clusters where *maximum* enrichment exceeded 2.0x. That tests whether one
sample is over-represented, when the informative quantity is the **spread** between the most
and least represented sample. It flagged only cluster 15 (5 cells) while missing cluster 11
— **1,264 cells, 200nm 1.59 against control 0.52, a 3.1x spread** — which is where any real
statistical power lies. Now reports max/min.

**QC diagnostics.** UMAP coloured by `Sample` shows the four interleaved with no monochrome
island, though note this plot cannot really distinguish "mixed" from "the last-drawn colour
covers everything" — it confirms Section 8's measured 86%, it is not itself the evidence.
UMAP by QC metric shows no region defined by mitochondrial content. There *is* a
depth-position relationship, but that is expected biology rather than artifact: monocytes are
large and transcriptionally busy, naive T cells and platelets are small and RNA-poor. Whether
those regions are cell types or artifacts is decided in Section 10 by marker coherence.

### Section 10 outcome — stage one (markers computed; mapping still to write)

Ran clean at execution counts 41-44. `rank_genes_groups` with `method='wilcoxon'`,
`layer='normlog'`, `use_raw=False`, `pts=True`, into `uns['rank_leiden']`.

**Provisional identities**, from the curated marker dot plot:

| cluster | cells | identity | cluster | cells | identity |
|---|---|---|---|---|---|
| 0 | 4,044 | cytotoxic CD8 T | 8 | 852 | B |
| 1 | 1,259 | CD8 effector memory T (`GZMK`) | 9 | 263 | T subset — confirm |
| 2 | 898 | MAIT / `KLRB1`+ T | 10 | 1,484 | NK (CD16+) |
| 3 | 5,773 | naive CD8 T | **11** | **1,264** | **monocytes** |
| 4 | 6,849 | naive CD4 T | 12 | 17 | myeloid, unresolved |
| 5 | 929 | **regulatory T** (`IL2RA` 57/5, `IKZF2`, `RTKN2`) | 13 | 1,153 | B |
| 6 | 3,332 | CD4 memory T | 14 | 139 | dendritic cells |
| 7 | 175 | NK | 15 | 5 | activated monocytes |

**Absent populations, confirmed by empty marker columns:** platelets (`PPBP`/`PF4`),
neutrophils (`CSF3R`; `FCGR3B` not in the dataset at all), proliferating cells
(`MKI67`/`TOP2A`), plasma cells. Empty columns are informative and normal for PBMCs.

**`pct_nz_group` vs `pct_nz_reference` earned its place immediately.** Cluster 3's top-ranked
genes by score were `RPS12` 100/99, `RPS3A` 100/98, `RPS28` 100/99 — ribosomal genes present
in every cell, which topped the ranking on a small difference in *level*. Ranking by score
alone would have presented them as markers. Only `NELL2` (93/43) discriminates.

**This is not an argument for having removed ribosomal genes earlier**, and the question was
raised and settled: Section 6's `seurat_v3` selection already excluded them from the HVG set
(1 of 85 chosen), so they never influenced PCA, neighbours, UMAP or Leiden. They appear here
only because marker discovery deliberately reads all 20,352 genes — the reason `subset=False`
was chosen. Deleting them from the object would remove activation-associated signal (the same
argument that stopped us filtering on ribo% in Section 3), change every cell's
`normalize_total` denominator, and strip Section 13's gene set. Their prominence in cluster 3
is itself biology: naive T cells are small and RNA-poor, so ribosomal transcripts are a larger
*fraction* of their library. **The fix is presentational only** — optionally hide ribo/mito
genes from the printed marker list, never from the data.

**Three attempts at an automated lineage/doublet verdict all failed, and the feature has been
removed rather than tuned again.** The failures were biological, not parametric:
1. Scoring each panel against its **cross-cluster average** makes abundant lineages invisible.
   T cells span 8 of 16 clusters, so the T average is 0.68 and cluster 4's score of 1.00 is
   only 1.47x — under a 2.0 threshold it printed "unresolved" for the clearest signal in its
   row. **The more common a cell type, the more invisible it became.**
2. Ratios on near-zero panels are meaningless: cluster 11 printed `Neutrophil 11.3x` on an
   absolute value of **0.01**.
3. "Two panels elevated" is not a doublet. Classical and Non-classical monocyte are one
   lineage; and cytotoxic T cells genuinely express the NK panel (`NKG7`, `GNLY`), so a T+NK
   cluster is real biology.

The cell now prints the absolute table and the top three panels per cluster, ranked, with no
verdict. **The dot plot is the instrument for this question**; a threshold rule over panel
averages is not.

**Three independent lines of evidence say the clusters are sound**, none of them a heuristic:
the curated dot plot shows textbook marker patterns; the expression **dendrogram** — given no
cell-type information — independently grouped the two B clusters, the three myeloid clusters,
and the T/NK clusters; and every violin in the stacked-violin plot is **unimodal**, where a
cluster of doublets or a mixed population would be bimodal. **There are no doublet clusters.**

**Stage two still to do:** confirm the CD4/CD8 split for clusters 3 and 4 with explicit
per-cluster means rather than dot sizes; resolve cluster 9; write the cluster→cell-type
mapping with 12 and 15 labelled unresolved rather than given confident names.

## Working process

Per section: read the corresponding tutorial cells first (what it does and why),
then write the adapted version in the notebook, discussing parameter choices as
we go. Every threshold and parameter gets a markdown justification — the rubric
explicitly asks for it. Work incrementally, one section at a time, verifying
output before moving on.

Concretely, for each section: **explain what we are about to do and why → write it →
run it → show the real output against what was expected → stop for review.** One
section at a time, never several batched together. Where a choice has more than one
defensible answer — QC thresholds above all — surface the tradeoff and let the user
decide rather than picking silently.

"Section" always means a row of the table below (1–16); "Task" always means an entry in
`CLAUDE.md`'s task list (1–6). Do not introduce a third numbering scheme.

Verifying a section means executing the notebook, not eyeballing the code. A cell that
runs without error proves nothing about whether it did anything — Section 2's first
attempt filtered zero cells while completing successfully. Run it with:

```bash
.venv/Scripts/python -m jupyter nbconvert --to notebook --execute \
  nanoplastics_pipeline.ipynb --output-dir <scratch> --output check.ipynb \
  --ExecutePreprocessor.timeout=5400
```

Sections 1–2 take ~10 minutes end to end, almost entirely Scrublet.

## Verification checklist

- **Every filtering step**: print before/after cell counts. `CLAUDE.md`'s
  explicit pitfall is a threshold that is present in code but numerically inert —
  confirming the code ran is not the same as confirming it filtered anything.
- **Layer integrity**: `layers['counts']` stays raw integers and
  `layers['normlog']` stays log-normalized/unscaled throughout. Confirm plots
  reference the layer they claim to.
- **After integration (8-9)**: UMAP colored by `Sample` should show samples
  well-mixed — not one cluster per sample — *without* erasing Leiden cluster
  structure. This is the batch-correction vs. biological-conservation tradeoff
  `CLAUDE.md` warns about; over-correction would erase the size-dependent signal
  we're looking for.
- **After annotation (10-11b)**: check per-cluster agreement across our three own
  sources (markers, `panhumanpy`, `celltypist`) and reconcile disagreements
  explicitly, recording the reasoning — none is automatically authoritative. Then
  separately report agreement against `predicted.celltype` and `CoDi`. Low
  agreement is a signal to revisit clustering resolution or marker
  interpretation, **not** a reason to silently adopt the external labels.
- **CellTypist inputs**: confirm it received log1p-normalised, 10,000-target-sum
  data (its documented requirement) and that `over_clustering` actually got our
  Leiden labels rather than defaulting to its own clustering.
- **Before any DE step**: confirm it operates on `layers['normlog']` or raw
  pseudobulk counts, never on scaled/regressed values. Same applies to
  `sc.tl.score_genes` — it must see `normlog`, or its expression-matched
  control-gene binning is meaningless.
- **Section 16**: check whether the per-cell score distribution within each cell
  type is unimodal or bimodal *before* summarizing with a mean — bimodality is a
  finding in itself, and a mean would hide it.
- **End to end**: the notebook should run top-to-bottom without manual patching
  once a section is finalized.

## Known limitations to state explicitly in the report

- **Single donor, no biological replication.** All 4 samples come from one donor,
  so there is no donor-level replication. Pseudobulk DE has n=1 per group. This
  must be flagged, not papered over.
- **First-pass QC is a working baseline.** Per `CLAUDE.md`, once Tasks 1-6 run
  end-to-end on standard QC, revisit with the mito%/complexity joint check before
  finalizing.
- **Pan-human vs. PBMC-specific reference.** `panhumanpy` uses the pan-human
  model, which may be less sharp on fine immune subtypes than a dedicated PBMC
  reference would be.
- **Legacy gene symbols — action required before Sections 14/16.** The data uses
  pre-2017 HGNC nomenclature throughout. A 14-pair spot check of well-known renames
  returned 14/14 legacy, 0/14 current: `IL8` (not CXCL8), `SEPT7`, `SEPT9`,
  `HIST1H1C`, `HIST1H4C`, `GPR56`, `CD97`, `FYB`, `C10orf54`, `MLLT4`, `KIAA0101`,
  `CCDC109B`, `WHSC1`, `FAM19A2`.

  Harmless for QC — ribosomal and mitochondrial gene families were never renamed, and
  KEGG_RIBOSOME still matches 85 of its 88 genes (the 3 misses are tissue-specific
  paralogs and a pseudogene genuinely absent from PBMCs). **But it fails silently in
  enrichment and scoring**: MSigDB ships current symbols, so `gseapy` and
  `sc.tl.score_genes` will simply not match our `IL8` against Hallmark's `CXCL8`, and
  will report enrichment over only the genes that happened to match, with no warning.
  Hallmark `TNFA_SIGNALING_VIA_NFKB` contains `CXCL8`, and `CLAUDE.md` names IL-8 as a
  candidate of interest. **Build an old→new symbol map before Section 14.**
- **Upstream cell filtering already happened, beyond UMI dedup.** Minimum values are
  ~500 counts and ~200 genes in every sample — a hard cut, not a natural taper, visible
  as a flat truncation at the bottom of Section 3a's violins. Almost certainly a Seurat
  `subset(nCount_RNA > 500 & nFeature_RNA > 200)`. `CLAUDE.md` states only that UMI
  dedup happened upstream; this is more than that, and it is why our own lower bound
  removes so little.
- **Sequencing depth is confounded with exposure.** Median counts run control 7,674 >
  mixture 6,090 > 40nm 5,269 > 200nm 4,334 — a 1.77× spread lining up with condition.
  With one donor this cannot be separated into biology (exposure suppressing
  transcriptional output) versus technical difference in capture between channels.
  **Direct input to Section 8:** depth differences confounded with condition are
  precisely the case where over-correction erases the size-dependent signal `CLAUDE.md`
  warns about, so the batch-correction method has to be chosen with this in view rather
  than rediscover it.
