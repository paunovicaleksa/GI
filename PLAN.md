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
| 3 | QC: mito% (`MT-`), ribo% (KEGG_RIBOSOME — tutorial's Broad URL now 301-redirects to `gsea-msigdb.org` but resolves correctly when redirects are followed, which `pd.read_table` does by default; verified returning the gene list), per-sample distribution inspection | 23-45 | Same method, but thresholds picked **per sample** from its own distribution, not one number reused everywhere |
| 4 | Concatenate the 4 samples after per-sample doublet+QC | 70-74 | Same concat pattern, no CSV loop needed |
| 5 | Normalization: `normalize_total` + `log1p`, stored explicitly as `layers['normlog']` | 46-52, 88 | Tutorial uses `.raw`; `CLAUDE.md` requires keeping both `counts` and `normlog` layers |
| 6 | HVG selection: `seurat_v3` flavor on raw counts, `batch_key='Sample'` | 91 (commented-out alt) | The tutorial's commented-out multi-batch line already matches `CLAUDE.md`'s preferred route — use it over the tutorial's dispersion-based default |
| 7 | Scale + optionally `regress_out` (total_counts, pct_mt, pct_ribo), PCA | 57-61 | Same as tutorial |
| 8 | Batch correction: **scanorama** → integrated embedding | 92-98 (scVI) | scanorama instead of scVI |
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
| 3 — QC metrics + filtering | 1 | **next** |
| 4-16 | 1-6 | not started |

Cells: 34,078 loaded → **31,839** after doublet removal.

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
exactly the size comparison Task 6 depends on. Still to do at Section 4: measure how
many genes survive the intersection, and check no biologically load-bearing gene
(`TNF`, `IL6`, `CXCL8`) is lost.

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
