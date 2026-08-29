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

### Task 4 result — monocytes are the only size-dependent composition shift

Section 12, on the outer join at resolution 0.3. Two measures, because proportions are closed
and neither measure is neutral: percentage points favour abundant types, CLR favours rare ones.
**Trust rows that survive both.**

| | 40nm | 200nm | mixture |
|---|---|---|---|
| Monocyte, percentage points | **-1.93** | **+3.46** | -0.16 |
| Monocyte, CLR | **-0.48** | **+0.44** | -0.17 |

Monocytes fall with 40nm particles, rise with 200nm (4.95% -> 8.41%, a 70% relative increase),
and are unchanged with both together. It clears the within-sample sampling floor in the two
single-size samples and not in the mixture.

**Why this row and not the others.** It is the only row whose two single-size exposures move in
**opposite directions**, and closure cannot produce that — a row inflated purely by other rows
shrinking moves the same way whenever composition shifts. The depth confound also runs against
it: 200nm is the shallowest sample and monocytes are RNA-rich, so capture bias would push this
result the other way.

**Weaker observations.** `CD4 T naive` falls in all three exposures (the only row consistent
across every exposure, but small in CLR). `B memory` falls in both single-size samples, not the
mixture. `Basophil` rises in all three with the largest CLR values in the table - on **15
control cells**, so an observation only, and already excluded from Section 13's DE.

**Two things that look like findings and are not.** `CD8 T naive` shows the largest raw changes
(+7.30 / +3.28 pp) but collapses to +0.06 CLR at 200nm: it is the arithmetic shadow of the
monocyte column growing, and it is also the row most exposed to the depth confound, since naive
T cells are RNA-poor. Everything else sits inside the sampling floor.

**The limitation that governs all of it:** one sample per condition. Nothing here separates
exposure from donor or channel. These are **described differences, not attributed effects**,
and that qualification travels with every number above.

**Open question for Task 6:** the mixture contains both particle sizes and shows neither
component's monocyte effect. An additive model does not predict that.

### Ambient RNA contaminates every non-myeloid cell type, and it cannot be corrected here

Found in Section 13, after the pseudobulk fold changes put a myeloid inflammatory program at
the top of *every* cell type's list. It is contamination, and the per-cell evidence is
unambiguous:

| `IL8` | mean copies among expressing cells | share of those carrying only 1-2 copies |
|---|---|---|
| **Monocyte** | **573.8** | **0.01** |
| B / CD4 T / CD8 T / CD8 T naive / MAIT / NK / Treg | 1.9 - 4.7 | **0.74 - 0.78** |

Monocytes carry a median of 478 copies per cell; lymphocytes carry 2, and three quarters of
the lymphocytes carrying any carry one or two. That is a Poisson sprinkle of an abundant
transcript across every droplet - free mRNA from lysed cells, captured beside whatever cell the
droplet held.

**Pseudobulk is structurally blind to it.** 3,000 lymphocytes at 2 copies each sums to 6,000,
indistinguishable in a pseudobulk column from 50 cells at 120 copies each. The step that solved
pseudoreplication destroyed the per-cell evidence separating soup from signal, which is why this
surfaced at Section 13 and not earlier.

**It dominates the rankings worst where the biology is quietest.** Myeloid inflammatory genes
take **15 of the top 30** upregulated genes in CD4 T at 200nm, 16 of 30 in CD8 T naive, 13 in
Treg - against **8 of 30 in the monocytes that produced them**. In a monocyte the same soup is a
rounding error beside 500 real copies and competes with a genuine multi-gene response; in a
lymphocyte there is nothing else, so the soup *is* the top of the list.

**No threshold separates it.** Three criteria were tried and each failed for a different reason:
counts among expressing cells is confounded with sequencing depth (it drops `ACTB` from CD8 T,
which are shallower); expression relative to the lowest-expressing cell type removes uniformly
expressed genes such as `B2M` by construction; and the two combined still miss `IL8` itself,
because a gene abundant enough to dominate the soup clears any threshold low enough to keep real
genes.

**Why it is unfixable on this data.** SoupX and CellBender do not infer the soup - they
**measure** it from the empty droplets in CellRanger's *raw* matrix. Our files are
`filtered_feature_bc_matrix`, CellRanger's filtered output, so the empty droplets were discarded
before the data reached us. The measurement cannot be reconstructed from the cells that
survived filtering. **For any future project of this shape, request the raw matrices at the
start.**

**What this constrains, precisely.** Pathway-level inflammatory enrichment must not be reported
for any lymphocyte population: a myeloid *program* appearing across seven unrelated cell types
at once, in monocyte proportions, including genes lymphocytes have no route to expressing
(`MMP1`, a matrix metalloproteinase; `CSF3`, G-CSF), is contamination - coherence identifies it,
not magnitude. It does **not** license calling any individual gene in any individual cell an
artefact: about 16 lymphocytes carry double-digit `IL8` at normal depth with no myeloid markers,
and T cells are documented to produce IL-8 in some subsets. The monocyte result is unaffected.

Section 13 therefore reports the diagnostic (`mean_nz`, `frac_le2`) beside every fold change
rather than filtering, so each gene can be judged on evidence instead of on a threshold that
does not exist.

### Task 5 result — the inflammatory response is monocyte-specific, and only monocytes

Section 14, preranked GSEA on Section 13's fold changes, Enrichr `KEGG_2021_Human`, 8 cell
types x 3 contrasts.

**Monocytes, clean on the ambient diagnostic:**

| pathway | contrast | NES | FDR | ambient_frac |
|---|---|---|---|---|
| IL-17 signaling | 40nm | **+2.69** | **0.000** | **0.00** |
| IL-17 signaling | mixture | +2.39 | 0.000 | 0.06 |
| IL-17 signaling | 200nm | +2.30 | 0.000 | 0.05 |
| TNF signaling | 40nm | +2.34 | 0.000 | 0.10 |
| TNF signaling | mixture | +2.25 | 0.000 | 0.15 |

`TNF signaling` is the pathway `CLAUDE.md` named in advance as a candidate. `IL-17 signaling`
is essentially uncontaminated in all three exposures.

**Clean pathways by cell type** (FDR < 0.25, `ambient_frac` <= 0.5):

| | 200nm | 40nm | mixture |
|---|---|---|---|
| Monocyte | 42 | 59 | 58 |
| B | 10 | 13 | 8 |
| CD4 T, CD8 T, CD8 T naive, MAIT, NK, Treg | **0** | **0** | **0** |

**899 significant pathways across the six lymphocyte types, none of them clean.** Every one is
ambient-driven. This is `CLAUDE.md`'s prediction - effects concentrated in phagocytic myeloid
cells, not lymphocytes - confirmed by measurement rather than by absence of evidence.

B cells pass the filter but not convincingly: their clean pathways sit at `ambient_frac`
0.38-0.47, against the 0.5 cutoff, where monocytes sit at 0.00-0.15. Direction there is `B cell
receptor signaling` **down** at 40nm and mixture.

### Composition and expression dissociate, in two places

| | monocyte proportion (Sec 12) | monocyte IL-17 response (Sec 14) |
|---|---|---|
| 40nm | **-1.93 pp** | **NES +2.69** (strongest) |
| 200nm | **+3.46 pp** | NES +2.30 |
| mixture | -0.16 (flat) | NES +2.39 |

**The expression response is uniform across all three exposures; the composition response is
not.** 40nm has fewer monocytes that are more activated; 200nm has more monocytes, also
activated; the mixture has unchanged numbers and a full response.

So the answer to Task 6's shared-versus-size-specific question **depends on which readout is
asked**: composition is size-specific (down at 40nm, up at 200nm, flat in the mixture), while
activation is shared across all three. This is exactly why `CLAUDE.md` requires Tasks 4 and 5
not be conflated, and it is the strongest structural finding in the analysis.

**The mixture also breaks the depth confound.** It is the deepest exposed sample (6,373 median
counts against 200nm's 4,440), so if sequencing depth drove the monocyte signal the mixture
should be weakest. It is not - NES +2.39 against 200nm's +2.30. That is the one monocyte
contrast where the confound can be ruled out rather than merely stated.

### Gene-set construction is a parameter, not a lookup

Recorded because it changed a headline. The same ranking, run against MSigDB's **KEGG MEDICUS**
(which replaced the retired `c2.cp.kegg` collection), gave a best monocyte result of q = 0.01.
Run against Enrichr's **`KEGG_2021_Human`** it gives q = 0.000.

MEDICUS fragments the inflammatory program into 4-gene disease-framed modules
(`PATHOGEN_HCMV_US28_TO_GNB_G_PI3K_NFKB_SIGNALING_PATHWAY`); small sets give unstable
enrichment scores. `KEGG_2021_Human`'s `IL-17 signaling pathway` holds 94 genes and produces a
far more stable score from identical input.

MEDICUS is also badly redundant: six "different" B-cell pathways at q < 0.01 proved to be one
87-gene set with **98% median pairwise overlap**. The initial reading of that was inflation -
six names for one finding. The larger cost was the opposite: **fragmenting one real program
into pieces too small to reach significance individually.**

**Cost of the switch:** Enrichr serves undated copies, so the library cannot be pinned to a
release the way an MSigDB GMT can. The access date is printed instead. Also worth recording:
`gseapy.get_library()` fails on this machine with an SSL certificate error (it uses `requests`
and certifi's bundle; the chain resolves only through the OS store), so the library is fetched
with `urllib` and passed to `prerank` as a dict.

### A third artefact: sequencing depth in the pathway rankings

Distinct from ambient RNA, and neither diagnostic finds the other.

Run with all genes, CD4 T's top hit was `TRANSLATION_INITIATION` at NES 2.13, FDR 0.005 - and
it **passed the ambient check cleanly** (`frac_le2` 0.03). Two things were wrong: the set is 78
of 80 ribosomal proteins and zero initiation factors, and ribosomal share tracks sequencing
depth, which is confounded with exposure:

| ribosomal % of pseudobulk | control | mixture | 40nm | 200nm |
|---|---|---|---|---|
| median depth | 7,852 | 6,373 | 5,477 | 4,440 |
| CD4 T | 29.0 | 30.1 | 34.4 | 34.5 |
| CD8 T naive | 30.9 | 32.4 | 38.4 | 38.6 |

The depth ordering, reversed. In a shallower sample low-expression genes drop out entirely and
the highest expressors take a larger share of what remains; size-factor normalisation does not
undo a change in composition.

**Response:** ribosomal and mitochondrial genes are dropped from the ranking handed to GSEA -
not from the data, and not from Section 13's fold changes. Narrower than Section 3's decision,
which declined to filter *cells* on ribo%.

**Not a fix.** The artefact is compositional: remove the ribosomal block and the next tier of
high expressors inherits the behaviour. Watch any pathway whose leading edge is dominated by
highly expressed genes.

**And a false alarm worth recording, because the wrong test looked convincing.** B-cell
oxidative phosphorylation was suspected of being the same artefact - its raw share falls with
depth in all eight cell types. But raw share mixes the depth effect with the biology. Against
each cell type's *own* typical gene, OXPHOS is **+0.371 in B cells (80th percentile) and -0.306
in monocytes (37th)**. A depth artefact moves every cell type the same way; a direction reversal
is biology. The three artefacts in this analysis each needed a different comparison - ambient
needed per-cell distributions, the ribosomal artefact needed across-sample shares, and this one
needed within-cell-type relative fold change.

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

Clustering at `resolution=0.4` produced 18 clusters whose sample composition is otherwise
near-proportional — almost every enrichment value sits between 0.85 and 1.35, which is the
strongest confirmation yet that Section 8 corrected batch without erasing biology. **The
clusters that break that pattern break it in the same direction:**

| cluster | cells | 40nm | 200nm | control | mixture | max/min |
|---|---|---|---|---|---|---|
| **9** (monocytes) | **1,466** | 0.56 | **1.47** | 0.80 | 0.86 | **2.6x** |
| **14** | 188 | 0.29 | **1.47** | 1.40 | 0.59 | 5.0x |
| **16** | 21 | 0.75 | **1.66** | 0.25 | 0.79 | 6.8x |
| **17** | 5 | 0.79 | **1.61** | 0.00 | 1.11 | — |

**Cluster 9 matters most**: 1,466 cells is enough for real statistical power in Sections
12-13, and it is visually an isolated island on the UMAP. Cluster 17 is five cells and can
only ever be an anecdote. Clusters 14 and 16 turned out to be **doublets** (Section 10), so
their skew is an artifact rather than a finding.

> **Numbering note.** All cluster IDs here are from the post-revision run. Before the Section 3
> revision the monocytes were cluster 11 (1,264 cells, 200nm 1.59 vs control 0.52, 3.1x) and
> the five-cell cluster was 15.

### Cluster 17 (was 15) — an acute inflammatory program in exposed samples only

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
the monocyte cluster, but not here.

**`IL8` is the deferred symbol problem becoming load-bearing.** This dataset uses pre-2017
HGNC names, so MSigDB's `CXCL8` will not match our `IL8`. The gene at the top of the most
interesting result in the pipeline is exactly the one that will silently fail to match in
Sections 14 and 16. The old→new symbol map is no longer a tidy-up task.

### Section 10 identifies the monocyte population (cluster 9, was 11)

Marker annotation reframes the finding above. Cluster 11 is not an anomalous inflammatory
cluster: it is **the only monocyte population in the dataset** (`CD14`, `LYZ`, `S100A8`,
`S100A9`, `FCGR3A`, `MS4A7`, `CST3`; classical-monocyte panel 1.17 against < 0.20 for every
other panel). Clusters 12 (17 cells) and 15 (5 cells) are also myeloid — the expression
dendrogram, given no cell-type information, groups 11, 12 and 15 on one branch.

So the observation restates as:

> **Monocytes are 7.38% of the 200nm sample against 4.03% of control and 2.81% of 40nm.**

Final percentages on the Section 11 labels, monocyte cells over that sample's total:

| | control | 40nm | 200nm | mixture |
|---|---|---|---|---|
| before the Section 3 revision | 2.33% | 2.35% | 7.06% | 4.25% |
| **after** | **4.03%** | **2.81%** | **7.38%** | 4.34% |

Two things changed and both matter:

1. **The 200nm-vs-control ratio fell from 3.0x to 1.8x.** The revision restored control
   monocytes that the upper depth bound had removed preferentially (53.7% of them, against
   25.5% in 200nm), so roughly a third of the original effect size was filter-made. The
   direction survives; the magnitude does not.
2. **40nm now sits *below* control.** Before the revision 40nm and control were
   indistinguishable (2.35% vs 2.33%). The ordering is now 200nm >> mixture ~ control > 40nm,
   a size-specific pattern rather than a dose response — and Task 6's actual question.

Still one donor and still descriptive. Section 12 tests it properly.

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

### Concatenation join (Section 4) — outer, changed from inner

The four files carry different gene sets (22,613 / 23,206 / 21,715 / 21,961) because each was
gene-filtered independently upstream at what is evidently `min.cells=3` — among the genes that
differ, the lowest detection anywhere is exactly 3 cells. **A gene "absent" from a sample was
therefore seen in at most 2 of its cells.** It is a near-zero measurement, not missing data.

That bounds the error each join makes, and the bound decides it:

- **outer**'s zero-fill misstates at most 2 cells per sample, always toward zero — it can
  understate a difference, never manufacture one.
- **inner** deletes the gene from all four samples, so it can never be tested at all.

The earlier rationale here argued outer "invents" a fold change. It does not: the invented
quantity is bounded by 2 cells out of several thousand, while inner's loss is total.

**Measured cost of the inner join:** 3,992 genes dropped (16.4% of the union). 2,419 of those
are absent from control but present in an exposed sample, and 288 are detected in 10+ cells.
The most-detected is `PTGES` (66 cells), prostaglandin E synthase — an NF-κB-driven
inflammation gene of exactly the kind Tasks 5-6 exist to find. **Switching to outer nets
+1,475 genes after `min_cells=10`: 20,387 → 21,862.**

`min_cells=10` on the merged object is the right place to make this call — one explicit
threshold applied to all four samples at once, stated in the notebook. The inner join was
making the same call implicitly, with a rule inherited from whoever preprocessed the files.

### Clustering resolution (Section 9) — 0.3, changed from 0.4

The outer join changed the stability sweep, and the change reversed the choice:

| resolution | clusters | stability (inner) | stability (outer) |
|---|---|---|---|
| 0.2 | 11 | 0.892 | 0.834 |
| 0.3 | 15 | **0.798** | 0.834 |
| 0.4 | 18 | 0.863 | **0.758** |

0.4's original justification was that its gap to 0.2 was within the method's own noise, since
0.3 dipped below both neighbours. **Under the outer join that dip is gone** — 0.2 and 0.3 tie
at 0.834 and 0.4 sits clearly below both. Read cold, nobody picks 0.4 from this table.

Keeping it would have meant defending a parameter with knowledge from the run being replaced.
Checked instead of assumed: at 0.3 the **only** population that merges away is the CD8
effector / effector-memory split. MAIT, Treg, Basophil, both B subsets and both monocyte
groups all stay separate. That split is also the annotation's most contested call (Section
11b: 911 of 1,789 cells disputed by both external schemes), and `COARSE` merges it into
`CD8 T` before Section 13 tests anything — so it costs nothing downstream.

**What stability cannot do:** it falls as resolution rises simply because more clusters mean
more boundaries to disagree about — 0.1 scores highest at 0.945 on 9 clusters, far too few for
PBMCs. It rules out unstable regions; it cannot pick a resolution on its own. The other half
of the criterion is whether the clustering resolves the populations Task 3 requires, and that
requirement comes from the project spec, not from a previous run.

### Cluster renumbering — what caught it

The outer join reordered every Leiden cluster while leaving the cluster *count* unchanged, so
the first assert in the map cell (count matches) passed and the hand-written `CELL_TYPES` was
silently wrong — the cluster it called `NK` was holding the monocytes. **The `CD14`-highest-in-
Monocyte sentinel in that cell is what caught it**, and without it every monocyte in Sections
12-16 would have been labelled NK.

Rule this establishes: any change upstream of Section 9 invalidates `CELL_TYPES`, and the map
must be re-derived rather than renumbered. Cross-tabulating the new clustering against a
preserved earlier labelling makes that a confirmation rather than a from-scratch job, which is
why `results/annotated_innerjoin.h5ad` is kept.

### Differential expression (Section 13) - pseudobulk, fold changes only, no p-values

The task list names no method for Task 5; `CLAUDE.md`'s DESeq2 instruction is accumulated
reasoning in its decisions section, not the project description. So the method was argued from
scratch.

**Pseudobulk is kept**, on its own merits rather than on that instruction: per-cell testing
treats 29,157 cells as independent observations of an exposure when they are four channels from
one donor, and returns p-values that answer "did we sequence a lot of cells?"

**Pseudobulk then collapses the design to n = 1 per condition**, and no method recovers a
p-value from that:

- `pydeseq2` per cell type **refuses outright** - *"there are no replicates to estimate the
  dispersion"*. Correct behaviour, and it is what ruled the original plan out.
- A pooled interaction design (`~ Sample + inC + Sample:inC`, 24 residual df) *is* estimable,
  but the replication it borrows is **between cell types**, not between donors, so its p-values
  answer a different question.
- Splitting cells into pseudo-replicates manufactures replication and re-introduces exactly the
  pseudoreplication pseudobulk exists to avoid.
- A fixed literature dispersion (edgeR's no-replicate workaround) states its assumption honestly
  but still assumes the quantity that matters.

**Decision: report effect sizes only.** Median-of-ratios size factors, then
`log2((exposed + k) / (control + k))` with `k` swept rather than asserted. A `padj` column would
be a number requiring permanent disclaimer - and the disclaimer lives in a markdown cell while
the column travels into tables and figures alone. Section 14's GSEA consumes a ranked list and
needs no significance threshold.

**One reporting trap found and fixed:** counting genes with `|log2FC| > 1` per cell type
correlates with cell count at about **-0.86** (fewer cells, noisier pseudobulk, more spurious
large ratios). That table is printed with its cell counts beside it and marked not comparable
across rows.

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

**Task 3 scope — settled, not open.** Asked directly whether "Azimuth PBMC reference" in the
task list requires running that reference ourselves: **it does not.** The deliverable is
`panhumanpy` + `celltypist` run by us, compared against the labels inherited with the `.h5ad`
files. So the pan-human/PBMC-specific gap is an **accepted scope decision** to state in the
report, not an unresolved shortfall — and the R/Seurat escalation above stays unexercised
unless something else forces it.

**There is no project description beyond `CLAUDE.md`'s task list.** Confirmed directly. Its
six lines are the whole spec — there is no rubric, point weighting, or fuller assignment
document to defer to. Where the task list is silent, the call is ours and belongs in this file
with its reasoning, rather than being resolved by guessing at unstated requirements.

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
| 9 — neighbours, UMAP, Leiden | 2 | **done, verified** (`RESOLUTION = 0.3`, 15 clusters) |
| 10 — marker genes + cluster annotation | 3 | **done, verified**; prose rewritten for the new numbering |
| 11 — automated annotation (CellTypist + panhumanpy) | 3 | **done, verified** |
| 11b — verification vs. inherited labels | 3 | **done**; new `EXTERNAL_TO_COMMON` assert not yet exercised |
| 12 — composition (Task 4) | 4 | **done**; effect measure implemented, cell not yet re-run |
| 13 - pseudobulk DE (Task 5) | 5 | **written**, cells not yet run |
| 14 - pathway enrichment (Task 5) | 5 | **written**, cells not yet run |
| 15 - size-specific effects (Task 6) | 6 | not started |

**Sections 4-12 have been walked through cell by cell** and their prose corrected against the
outer-join / resolution-0.3 run. The recurring defect was a number hardcoded into markdown
that the code recomputes correctly — the code was right throughout. Two exceptions worth
naming, because they are different failure classes:

- **Section 11's CellTypist markdown described an action the code no longer takes**, and
  argued for the opposite conclusion to the one now recorded (per-cell reassignment of the
  naive T cluster to CD4). Stale *numbers* are visible by comparing against output; a removed
  *action* leaves no trace in output at all. It was caught by reading adjacent cells in order,
  not by checking figures.
- **`PLAN.md` cited an agreement figure ("excluding NK, 93.2%") that no cell computes.** It was
  removed rather than given code, because the by-label disagreement table already shows what
  it was for, and an "agreement excluding X" headline invites the coarsening it hides:
  agreement rises monotonically as categories are dropped.

Cells: 34,078 loaded → **31,839** after doublet removal → **29,157** after QC (8.4%
removed, spread 1.6 pp across samples). **Unchanged by the join**, as expected.
Genes: 24,380 union → **21,862** after `min_cells=10` on the outer join
(was 20,388 → 20,387 under the inner join).

Preserved for before/after comparison, both git-ignored:
`results/annotated_innerjoin.h5ad` (the full inner-join object, all labels) and
`results/executed_innerjoin.ipynb` (every printed output of that run).

Sections 1-7 have been run **top-to-bottom on a restarted kernel** with no errors
(Sections 1-5 at execution counts 1-24; Sections 6-7 at 25-30; Section 8 at 31-34;
Section 9 at 35-40).

**That claim did not extend to Section 9, and the gap was real.** A cold
`nbconvert --execute` on 2026-08-29 failed there with
`NameError: name 'RESOLUTION' is not defined`. Cell 57 defines `RESOLUTIONS` — plural, the
sweep list — while the cell running the final clustering read `RESOLUTION`, singular, which
**nothing in the notebook ever assigned**. It had been living in the interactive kernel's
memory, left behind by a cell since edited or deleted, so every interactive run succeeded
while the notebook was never actually reproducible. Fixed by adding `RESOLUTION = 0.4` beside
the `sc.tl.leiden` call — the value that cell's own comment already justifies.

This generalises past the one name: **a section can be recorded here as verified on the
strength of output a cold kernel cannot reproduce.** It is also the likely explanation for
interactive runs completing much faster than cold ones — a warm kernel is skipping work, not
running it faster.

**With that fix, Sections 1-11 now pass a cold run end to end** (~23 min, 2026-08-29),
reproducing every number recorded above exactly: 29,157 cells, 20,387 genes, 18 Leiden
clusters at resolution 0.4, 94.0% Azimuth-vs-CoDi, 88.1% / 87.9% against the inherited labels
and `CoDi`, 45.7% CellTypist-vs-panhumanpy on cluster 3, and the same 8 testable coarse groups.
Identical reproduction was not a given — Leiden numbering is not stable across runs, so a
re-ordered clustering would have invalidated `CELL_TYPES` silently; the `CD14`-highest-in-
Monocyte assert in that cell is what would have caught it, and it passed.
A write-only checkpoint is saved at `results/combined_after_qc.h5ad` (~1 GB, git-ignored).
Nothing auto-loads it — a checkpoint that reloads itself is the same bug as a stale
output. Restore by hand with `sc.read_h5ad` if a kernel dies.

A **second** write-only checkpoint is saved at `results/annotated.h5ad` by a new cell at the
end of Section 11b. It exists because the first one is not enough: it predates clustering and
annotation, so `obs['leiden']` and `obs['cell_type']` live nowhere on disk without it, and
Sections 12-16 all group by `cell_type`. A dead kernel therefore used to cost a full re-run of
Sections 5-11 — scanorama's sigma sweep, Leiden's resolution sweep and both annotation models
included. This was discovered the hard way when a jupyter-lab restart discarded the session
that produced the Section 11b results recorded above.

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
`log1p_total_counts` and `log1p_n_genes_by_counts` at 5 MADs both ends (**revised to
lower-only — see the revision below**);
`pct_counts_mt` at 3 MADs **upper only**, since low mito% is healthy and a two-sided
rule would flag the best cells. Depth and gene counts use the log scale because both
are log-normal — a symmetric rule on raw counts flags far more cells above the median
than below, purely from skew.

Removed **2,682 of 31,839 cells (8.4%)**, spread 7.7–9.4% across samples. The mito filter
dominates (6.5–8.2%); the low-depth end removes 1.1–3.2%. (Figures after the revision below;
the original two-sided rule removed 3,403 cells, 10.7%.)

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

### Section 3 revision — the upper depth bound was removed

**Found while checking Section 10's annotation**, by asking a question Section 3 could not ask
at the time: how many cells of *each type* did QC remove? The overall removal rate was even
across samples (10.0–11.6%) and looked healthy. Per cell type it was not.

| cell type (inherited labels) | before QC | after QC | kept |
|---|---|---|---|
| Cytotoxic T | 7,285 | 6,375 | 87.5% |
| CD4 T | 22,153 | 18,819 | 85.0% |
| B | 2,607 | 1,955 | 75.0% |
| **CD14 monocyte** | **2,023** | **1,280** | **63.3%** |

**Cause: the upper tail of the two depth rules.** It flagged 14.0% of monocytes and 0.0% of
cytotoxic T cells. That is not doublet detection — it is a size filter. Monocytes carry a
median 11,476 counts against 3,663 for a cytotoxic T cell, so a bound set by a population that
is 65% T cells removes monocytes for being monocytes. Of the 743 monocytes lost, 639 were
QC-flagged and only 104 came from Scrublet; **309 were caught by the upper depth/genes tail
alone**, with no mito involvement.

**Why it matters beyond lost cells.** The bound bit unevenly across samples — 42.6% of control
monocytes against 12.3% of 200nm ones, because control monocytes sit at a median 26,701 counts.
Task 4 compares monocyte *proportions* between those samples, so a filter that removes them at
different rates per sample manufactures part of the difference it is meant to measure.

**The guard that should have caught it.** Section 3b already checks that the removal rate is
even across samples, and it read 10.0–11.6%. That statistic is computed over all cells, so a
bias confined to the 6% of the dataset that is monocytes averages away. Not an inert check —
a check at the wrong resolution, which is the subtler version of `CLAUDE.md`'s named pitfall.

**Fix applied:** both depth rules go to `'lower'` only; `pct_counts_mt` unchanged at 3 MADs
upper. Doublet removal is left to Section 2, which has per-cell evidence instead of a size
proxy — what `CLAUDE.md` means by calling a max-genes cutoff a weak sole doublet proxy.
Recovers 309 monocytes, 150 B, 247 CD4 T, 2 cytotoxic T.

**Still outstanding after the fix:** monocytes still lose 16.5%, now dominated by the mito
rule. That is bonus #1 (mito%/complexity joint check) and is unchanged by this revision.

**Consequence:** everything from Section 3b onward was re-run. Leiden numbering is not stable
across runs, so all cluster IDs recorded below are from the post-revision run.

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

3,000 HVGs flagged of 20,387 genes; `subset=False`, so nothing was discarded.

**The mean–variance correction demonstrably worked:** only **1 of 85** ribosomal genes
and **0 of 13** mitochondrial genes were selected. Those are among the most highly
expressed genes in every cell, so a naive variance ranking would have put them at the
top — excluding them is the whole purpose of `seurat_v3`'s loess fit.

**1,711 genes are variable in all four samples** (`highly_variable_nbatches == 4`). That
gives `n_top_genes` a data-driven floor it otherwise lacks: 3,000 reaches ~1,300 genes
past the reproducible core. This is the softest parameter in the pipeline — unlike
`target_sum=1e4` (forced by CellTypist) or the QC thresholds (derived from the data),
it is conventional. The tutorial's stated heuristic is "3k is fine for 10k cells"; we
have 29,157, so 3,000 is if anything conservative. If a stability check is ever wanted,
the test is whether clustering changes materially at 2,000 vs 3,000 — the same standard
Section 9 applies to Leiden resolution.

**9 of 12 cell-type markers were flagged.** The three missing — `CD3E`, `CD8A`,
`FCGR3A` — are expected rather than a problem: `CD3E` is expressed in every T cell, so
it has a high mean and unremarkable standardized variance. **A good marker is not the
same thing as a variable gene.** HVGs drive the embedding; Section 10 reads all 20,387
genes for markers, which is exactly why `subset=False` matters.

**New dependency:** `flavor='seurat_v3'` requires `scikit-misc` (0.5.2) for its loess
fit. Same failure shape as `scikit-image` for Scrublet — not a hard scanpy dependency,
so scanpy imports fine and the call dies at runtime. Now pinned in `requirements.txt`;
it needs `numpy>=1.26.4`, exactly the version TensorFlow's ceiling pins us to.

### Section 7 outcome

Ran clean in the notebook at execution counts 26-30. `X_pca` is (29,157 x 50), built from
the 3,000 HVGs scaled on a throwaway copy with `max_value=10`; `X_pca_sel` is (29,157 x 30)
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

Ran clean at execution counts 31-34. `adata.obsm['X_scanorama']` is (29,157 x 30);
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

### Section 9 outcome (post-revision re-run)

`adata.obs['leiden']` holds **18 clusters** at `resolution=0.4`, on a neighbour graph over
`X_scanorama` with `n_neighbors=15`, 29,157 cells.

**`use_rep='X_scanorama'` is asserted, not assumed.** Omitted, `sc.pp.neighbors` falls back to
`adata.obsm['X_pca']` - 50 *uncorrected* components - discarding Sections 7 and 8 while
running perfectly and plotting beautifully.

**One of the two resolution methods failed, and this is recorded rather than glossed.**

- **Plateau test: no signal.** Cluster count rises almost linearly across 0.1-1.5
  (9, 13, 16, 18, 21, 21, 24, 28, 30, 36) with no flat span.
- **Subsampling stability: usable, not decisive.** Mean ARI against a re-clustered 90%
  subsample reads 0.856 / **0.892** / 0.798 / **0.863** / 0.773 across 0.1-0.5. The top is 0.2
  with 0.4 second - but **0.3 sits below both its neighbours**, and stability cannot genuinely
  be non-monotonic in resolution. That dip is the method's own noise across three repeats, of
  order +/-0.05, which makes 0.892 and 0.863 inseparable.

**Before the QC revision this sweep peaked cleanly at 0.3 (0.913) and 0.4 (0.912).** Recovering
~700 previously discarded cells flattened the peak. The same metric on nearly the same data
gave a much more confident-looking answer before, and the extra confidence was not real.

**`resolution=0.4` kept.** The tie broke on what later tasks need: 18 clusters rather than 13
keeps MAIT (890) and the basophils (157) as their own populations rather than absorbed into
neighbours; over-splitting is recoverable in Section 10 by reading markers, under-splitting is
invisible.

**The tiny clusters are not a granularity artifact.** `smallest` is 5 cells at *every*
resolution tested including 0.1, and 3 clusters sit under 50 cells even at 0.1. Coarser
clustering merges the large groups, not the specks.

**Bias worth recording:** stability structurally favours fewer, larger clusters. Taken as
`argmax` it would have argued for `resolution=0.2` and the loss of MAIT and the basophils.

**Cluster sizes** run 7,371 (25.3%) down to 5 (0.02%), six clusters under 200 cells.

**A diagnostic of ours was mis-specified and has been corrected.** The composition cell
originally flagged clusters where *maximum* enrichment exceeded 2.0x. The informative quantity
is the **spread** between the most and least represented sample. Now reports max/min.

### Section 10 outcome (post-revision re-run)

`rank_genes_groups` with `method='wilcoxon'`, `layer='normlog'`, `use_raw=False`, `pts=True`.
`adata.obs['cell_type']` written from `CELL_TYPES`.

| cluster | cells | identity | cluster | cells | identity |
|---|---|---|---|---|---|
| 0 | 3,995 | CD8 T effector | 9 | 1,466 | **Monocyte** |
| 1 | 1,789 | CD8 T effector memory | 10 | 130 | T unresolved |
| 2 | 48 | unresolved (stress program) | 11 | 890 | MAIT (`SLC4A10` 64%) |
| 3 | 7,371 | **T naive, CD4/CD8 unresolved** | 12 | 1,211 | B naive |
| 4 | 5,209 | CD4 T naive | 13 | 157 | **Basophil** (was "DC") |
| 5 | 917 | **Treg** (`FOXP3` 33%) | 14 | 188 | **doublet** |
| 6 | 3,077 | CD4 T memory | 15 | 29 | **doublet** |
| 7 | 1,065 | B memory (`CD27` 42%) | 16 | 21 | **doublet** |
| 8 | 1,589 | NK | 17 | 5 | myeloid, unresolved |

**Three clusters are doublets - 238 cells (0.8%).** Settled by doublet score, not argument:
Scrublet score ~2x the dataset median (0.147 / 0.145 / 0.169 against 0.072) with 36% / 34% /
52% of their cells in the top score decile, *and* median depth 22,927 / 23,853 / 31,762 against
a dataset median of 5,541. A large real cell type is high on depth alone - cluster 9
(monocytes) has median depth 12,861 but only 21% in the top score decile.

**Why they are here at all, and why that is an improvement.** Section 2 thresholded Scrublet at
50% sensitivity by design, so roughly half of all doublets were always going to survive, and
Scrublet cannot detect two cells of the same type in one droplet at all. Before the Section 3
revision these leftovers were removed *incidentally* by the upper depth bound - the same cut
that removed 309 real monocytes. One net was catching two unrelated things. Removing doublets
at cluster level is more targeted and is standard practice: doublets share a mixed profile and
cluster together, whereas a per-cell threshold must cut through overlapping distributions.

**Independent support from the external labels.** Clusters 14, 15 and 16 are exactly where the
inherited labels and CoDi disagree with each other most, and where CoDi's top-label share is
lowest (57%, 48%, 62%, against 90-100% for settled clusters). Contrast cluster 13: CoDi's top
share is also low (50%) but its doublet score is ordinary and its depth is *below* the dataset
median - low external agreement there means a real cell type the references cannot name.

**Genuinely unresolved: 183 cells, 0.6%** - clusters 2 (48, a heat-shock/housekeeping stress
program), 10 (130, no marker with a real gap) and 17 (5, myeloid).

**`pct_nz_group` vs `pct_nz_reference` earned its place immediately.** Cluster 3's top-ranked
genes by score were ribosomal - `RPS3A` 100/98, `RPS12` 100/99 - present in every cell,
topping the ranking on a small difference in *level*. Ranking by score alone would have
presented them as markers.

**Three attempts at an automated lineage/doublet verdict all failed and the feature was removed
rather than tuned.** (1) Scoring each panel against its cross-cluster average makes abundant
lineages invisible - T cells span 8 of 18 clusters. (2) Ratios on near-zero panels are
meaningless - `Neutrophil 11.3x` on an absolute 0.01. (3) "Two panels elevated" is not a
doublet: classical and non-classical monocyte are one lineage, and cytotoxic T cells genuinely
express the NK panel.

**The basophil correction, and how the error was made.** Cluster 13 was called `Dendritic cell`
on `FCER1A` plus a pDC panel score of 0.47. That panel averages `LILRA4`, `CLEC4C` and `IL3RA`,
and the score came almost entirely from `IL3RA` at 70% while `LILRA4` sat at **0%**. `IL3RA` is
CD123, carried by basophils as well as pDCs. Per gene the cluster reads `HDC` 56%, `CPA3` 33%,
`CLC` 30%, `GATA2` 13%, `CD1C` 0% - basophils. **A panel mean can be driven by a single shared
gene**; the per-gene table caught it and CellTypist prompted the check.

### Sub-clustering outcome - two splits tested, neither adopted

**Cluster 3 holds no separable naive CD4 population.** Sub-clustered at 0.3 it gave 2,144 /
1,773 / 3,454 cells, all three reading `CD4` 2-5% and `CD8B` 41-72% - differing in degree of
CD8 detection, not identity.

**That test is weak for this question and the weakness is structural.** Section 6 flagged
`MS4A1`, `CD14`, `LYZ`, `NKG7`, `GNLY`, `FCER1A`, `PPBP`, `LILRA4` as HVGs but **not** `CD3E`,
`CD8A`, `IL7R` or `FCGR3A`. The neighbour graph is built from HVGs, so lineage genes for the
CD4/CD8 distinction barely influence the geometry Leiden partitions.

**Cluster 9 does not split into classical and non-classical monocytes, and that result is
solid.** 14 sub-clusters came back, all `CD14`+/`S100A8`+/`FCGR3A`- with `IL8` at 99-100% -
over-splitting on a gradient. The clustering-independent check agrees: across all 29,157 cells
`FCGR3A` is detected in 5.6% but only **20 cells carry 5+ counts**, most of them NK cells.

**Neither split adopted.** `sub_3` and `sub_9` remain in `obs` as evidence.

### Section 11 + 11b outcome - Task 3 complete

**CellTypist** (`Immune_All_Low`, 98 immune subtypes, per cell, `majority_voting` over our
Leiden clusters). Feature overlap 71.7% - the legacy-symbol tax, visible as `IL8` present in
our data and absent from the model while `CXCL8` is the reverse. It confirmed the marker
annotation on every settled cluster, named MAIT and Treg correctly, gave the doublet clusters
no coherent answer, and **found the basophils**.

**Pan-Human Azimuth** (`panhumanpy`, 3 broad / 14 medium / 45 fine, feature overlap ~75%).
Coarser on immune biology - misses MAIT entirely, returns `False` for 14% of all cells, rising
to 27-67% in the hardest clusters. That refusal is informative and is a property neither
inherited labelling has. **It holds `cells_meta` by reference to `adata.obs`, so its columns
are written into our object in place** - which is why a column diff after `pack_adata()` finds
nothing.

**11b, against the labels that shipped with the files** - final `cell_type` mapped to the
coarsest vocabulary all schemes can express, the mapping written out in the cell:

| | agreement |
|---|---|
| inherited `predicted.celltype` | **88.8%** of 29,157 cells |
| `CoDi` | **88.7%** |

**The disagreements are concentrated, not diffuse**, which is the result that matters. The
by-label table below is the evidence for that, and it is preferred to any single
"agreement excluding X" headline: agreement rises monotonically as categories are dropped, so
such a number reads as more reassuring than it is.

| our label | disputed | why |
|---|---|---|
| `NK` | **1,624 of 1,630** | their schemes call NK "Cytotoxic T"; `CD3E` is in 7% of these cells |
| `CD8 T effector` | 942 of 5,610 | the function axis - they split cytotoxic from CD4, we split by lineage |
| `Basophil` + `Doublet` | 205 | categories they do not have |
| `B memory` + `B naive` | 192 | minority disagreement inside clusters both schemes otherwise agree on |
| everything else | 292 | MAIT 113, CD4 T memory 96, Monocyte 80, Treg 2, CD4 T naive 1 |

`CD8 T naive` contributes essentially no disagreement, because it maps to `other T` on our
side and they call it `CD4+ T cell`, which maps there too - the relabel changed what we call
those cells, not which external category they fall in.

Diffuse disagreement would mean something systematic was wrong. Concentrated disagreement in
categories the references demonstrably lack is the expected shape.

### The two label columns, and why there are two

`obs['cell_type']` holds the 16 fine labels and is what **Task 4 / Section 12** counts.
`obs['cell_type_coarse']` holds 9 groups plus `Exclude` and is what **Task 5 / Section 13**
tests, because pseudobulk sums each group into one column per sample and control, the smallest
sample, binds.

Rule applied: merge only where a group cannot support a test alone, or where no mechanism would
make the split respond differently.

| coarse group | from | control cells | in DE? |
|---|---|---|---|
| **Monocyte** | Monocyte | 228 | yes - where the effect is expected |
| Basophil | cluster 13 | 15 | **no** - composition only |
| T naive | cluster 3 | 1,268 | yes - lineage unresolved, but a real population |
| CD4 T | CD4 naive + memory | 1,749 | yes |
| CD8 T | CD8 effector + eff. memory | 1,182 | yes |
| Treg / MAIT / NK / B | as labelled | 152 / 163 / 336 / 483 | yes |
| Exclude | Doublet, unresolved | - | dropped |

**Basophils are deliberately not merged into Monocyte.** It would make them testable, but a
granulocyte in the phagocyte pseudobulk labelled "myeloid response" blurs precisely the signal
this project exists to find.

**`MIN_PSEUDOBULK = 50` is a rule of thumb, not a derived threshold.** Below roughly 50 cells a
pseudobulk column reflects which cells happened to be captured more than the sample's biology.
Stated in Section 10 so Section 13 inherits an explicit list rather than discovering the problem
as a failed model fit.

## Additional analyses — candidates outside the project description

The task list has six entries and nothing else. Anything below is an **addition**, to be
chosen from rather than all implemented; the working target is 3-5. Kept here so the core
pipeline stays exactly aligned to the six tasks and additions are a deliberate choice rather
than scope drift.

Inherited from `CLAUDE.md`'s bonus list:

1. **Mito%/complexity joint QC rework**, plus per-cell-type mito% as an exposure readout.
   Needs Task 3 done, and supersedes the first-pass QC once done.
2. **Leiden vs. Louvain robustness comparison.**
3. **scVI native DE module** - not applicable, scVI was avoided throughout.
4. **Marker-gene retention / downsampled-reference accuracy check** on annotation quality.
5. **scANVI vs. Azimuth annotation comparison.**
6. **van den Brink dissociation-stress signature scoring** - tests whether high-mito% cells
   also score high on heat-shock/immediate-early genes, which is what actually resolves
   item 1 (technical artefact vs. real biology).

Added during the build:

7. **MSigDB Hallmark GSEA** (Section 14) and **per-cell Hallmark scoring** (Section 16).
   These are one candidate, not two, and the reason is the pairing: Section 14's GSEA
   collapses every cell in a group to a single ranked list before testing, while
   `sc.tl.score_genes` gives every cell its own value. Run on **the same gene set** they
   answer different questions - a bimodal per-cell score inside monocytes would be direct
   evidence of uptake-dependent response, which no pseudobulk method can show. Hallmark is
   the right set for it: 50 curated non-redundant sets, and the best symbol match of the four
   libraries (90.7% against GO:BP's 75.5%). Excluded from core because the task list names
   GO/Reactome/KEGG and not Hallmark, and because Task 6's deliverable is Section 15.
8. **Ambient RNA correction.** Not implementable on this data - SoupX and CellBender need
   the empty droplets in CellRanger's raw matrix, and we have only the filtered output. Listed
   so the gap is a recorded decision rather than an omission, and as the reason to request raw
   matrices at the start of any future project of this shape.
9. **Uncorrected-embedding clustering control.** Cluster on `X_pca_sel` (preserved by Section
   8) and cross-tabulate against the corrected labels. Measured once during the build:
   monocytes recover at 99.7% purity and NK at 99.2% without any batch correction, so the
   headline result does not depend on scanorama. `Treg` and the B subsets do not separate
   without it. Cheap, and it converts an assurance into a measurement.

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
  once a section is finalized — and the evidence has to be a **cold** run
  (`nbconvert --execute` on a fresh kernel), never an interactive pass. An interactive
  pass cannot prove this, because the thing it fails to catch is precisely what a warm
  kernel supplies.
- **Leaked kernel state**: a cell can work only because a long-lived session still holds
  a variable whose defining cell was edited away. `RESOLUTION` was exactly this, and **no
  runtime check in this list would have caught it** — the cell ran, produced correct
  output, and was recorded above as verified. Reading the cell sources for names that are
  never assigned, or assigned only later, catches it before a cold run does.

## Known limitations to state explicitly in the report

- **Single donor, no biological replication.** All 4 samples come from one donor,
  so there is no donor-level replication. Pseudobulk DE has n=1 per group. This
  must be flagged, not papered over.
- **First-pass QC is a working baseline.** Per `CLAUDE.md`, once Tasks 1-6 run
  end-to-end on standard QC, revisit with the mito%/complexity joint check before
  finalizing.
- **Cluster 2 is naive CD8 T cells, at cluster level only — and this reverses an earlier
  call.** The cluster (6,674 cells, 22.9%) was labelled `T naive (CD4/CD8 unresolved)` on the
  grounds that two per-cell models could not agree. That reasoning weighed the models against
  each other and ignored the genes, and the genes are not ambiguous:

  | | `CD4` | `CD8B` |
  |---|---|---|
  | cluster 2 | **2.9%** | **54.8%** |
  | populations we call CD4 (naive, memory, Treg) | 22-24% | 1.8-3.9% |
  | CD8 T effector | 2.6% | 38.0% |

  `CD4` captures poorly in 10x — it reaches only ~24% even where we are confident — so its
  **absence** is weak evidence. `CD8B` captures normally, and its **presence** at 54.8%
  against a 3.9% background is positive evidence. Note it exceeds the CD8 effector cluster's
  38.0%, which is exactly what naive CD8 cells look like. `CCR7` 93.1% and `IL7R` 94.5% place
  it firmly in the naive compartment.

  **The two models failed differently, and that asymmetry is the finding.** CellTypist called
  it 67% `Tcm/Naive cytotoxic T`; Pan-Human Azimuth called it 87.5% `Naive CD4 T`. On every
  cluster where the genes give an independent answer — MAIT, Treg, Basophil, NK — CellTypist is
  right and `panhumanpy` is weaker or wrong: it labels our CD4 **memory** cluster `Naive CD4 T`
  at 70%, and returns `False` on MAIT and Basophil entirely. Its errors run consistently toward
  CD4, the same direction as the inherited labels — a shared prior rather than a second
  independent measurement. So the 42.1% agreement between the two models is not "two good
  methods disagree, truth unknowable"; it is one method reading `CD8B` and one defaulting.

  **What is still not known, and is not attempted.** 48% of these cells express *neither* `CD4`
  nor `CD8A`/`CD8B`. For half the cluster there is no per-cell information at all, which is why
  the two models correlate below chance. The claim is therefore **cluster-level only** —
  aggregation over 6,674 cells recovers a fraction that no individual cell carries, because
  dropout dilutes a cluster-level fraction but cannot invent one. Per-cell CD4/CD8 assignment
  is not recoverable from this data.

  **Consequence to state explicitly:** this flips the dataset's CD4:CD8 ratio relative to both
  inherited labellings, which call the cluster 100% `CD4+ T cell`. With cluster 2 as CD8, CD8
  is ~42% against CD4 ~34.5% — atypical for PBMCs, where ~2:1 CD4-dominant is the norm. **The
  ratio must be reported as a consequence of this call, not as an independent finding**, and
  the disagreement with the inherited labels stated alongside it.

  `CD8 T naive` is its own coarse group rather than merged into `CD8 T`: same lineage, very
  different activation state, and both clear the pseudobulk size check alone (control n = 1,139
  and 1,152).
- **There are no dendritic cells in this dataset, and the cluster first called DC is
  basophils.** `CD1C` reaches 0-5% in every cluster and `LILRA4` never exceeds 6%, so neither
  cDC2 nor pDC is present. Cluster 13 (157 cells) reads `HDC` 56%, `CPA3` 33%, `CLC` 30%,
  `GATA2` 13% and is labelled `Basophil`. Also absent: platelets (`PPBP` in 257 cells, only 41
  with 5+ counts), neutrophils, plasma cells, proliferating cells.
- **Non-classical (CD16+) monocytes are absent, so Task 3's classical/non-classical split
  cannot be made.** `FCGR3A` is detected in 5.6% of cells but only 20 carry 5+ counts, most of
  them NK cells. Sub-clustering the monocytes returned 14 groups, all `CD14`+/`FCGR3A`-.
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
