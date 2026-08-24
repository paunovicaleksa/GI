# Project: Nanoplastics scRNA-seq Immune Response

## Overview
Single-cell RNA-seq analysis of human PBMCs exposed to carboxylated polystyrene
nanoparticles (PSNPs), investigating how particle size shapes the immune response.

- Sample 1: 40 nm PSNPs
- Sample 2: 200 nm PSNPs
- Sample 3: 40 nm + 200 nm mixture
- Sample 4: control (no exposure)
- All 4 samples from one donor. Pre-processed `.h5ad` files (AnnData format).
- `*_CoDi_KLD.csv` files are a reference for optional comparison, not primary input.

Stack: Python, Scanpy/AnnData. R/Seurat only if doing Azimuth annotation that way.

Implementation reference: following a Sanbomics (YouTube) Scanpy tutorial series for
the actual coding pass. Check what it already covers before re-deriving something from
scratch — e.g. gene signature scoring is expected to already be in there. Treat this repo's
plan as the task/rationale layer and the tutorial as the how-to-code-it layer; reconcile
naming/parameter choices between the two rather than assuming either is automatically right.

## Task list (from project spec)
1. QC & preprocessing — filter cells, normalize, select variable genes. Justify thresholds.
2. Integration & clustering — batch-correction method of choice; UMAP + clusters.
3. Cell type annotation — marker genes + Azimuth PBMC reference.
4. Composition analysis — cell type proportions across the 4 samples.
5. Differential expression — each cell type, each exposed sample vs. control; pathway
   enrichment (GO/Reactome/KEGG).
6. Size-specific effects — unique to 40nm, unique to 200nm, shared, or mixture-only.

## Key methodological decisions and constraints

**QC (Task 1) — first pass**
- Do standard QC: pick thresholds by inspecting each sample's own distribution
  (count depth, genes detected, mito%), not by copying generic literature numbers
  unchecked. That's the only non-negotiable part of the first pass.
- Apply mito% filtering normally in this pass — don't do the joint mito%/complexity
  distress analysis or hold back on filtering high-mito cells yet. That level of care
  is deferred (see Bonus candidates) until the full pipeline runs end-to-end once.
- UMI dedup already happened upstream of these `.h5ad` files — not part of this task.
- Include doublet detection (Scrublet/scDblFinder) as its own step, not just a hard
  max-genes cutoff (weak as a sole doublet proxy).
- **Revisit note:** once Tasks 1-6 run once on standard QC, come back and redo QC
  with the mito%/complexity joint check before finalizing — treat the first full run
  as a working baseline, not the final version.

**Normalization / HVG (Task 1)**
- Standard route: `normalize_total` + `log1p`, then HVG selection via `seurat_v3`
  flavor (regularized NB regression — more stable under zero-inflation than the
  older dispersion/binning `seurat` flavor).
- Pearson residuals are an alternative that folds normalization + HVG weighting into
  one step (favors sparse, cell-type-restricted genes) — good for embedding/clustering,
  NOT for DE testing.
- **DE testing must run on log-normalized data (or raw/pseudobulk counts), never on
  Pearson residuals or scaled/z-scored data.** Keep a `layers["counts"]` (raw) and
  `layers["normlog"]` (log-normalized, unscaled) around — don't overwrite `.X` in place
  with scaled values without preserving both.

**Batch correction (Task 2)**
- Required to pick and justify one method (Harmony, Scanorama, scVI, BBKNN, etc.).
  Tradeoff axis: batch correction vs. biological conservation — don't over-correct
  and erase the actual size-dependent signal.
- If scVI is chosen: it also has a native Bayesian DE module (posterior-based,
  Bayes factor or "change mode") — usable as a bonus complementary DE method with
  near-zero extra cost once the model is trained for integration anyway.

**Clustering (Task 2)**
- Leiden preferred over Louvain (guarantees connected communities). Resolution
  controls granularity — no single correct value, justify by stability.

**Annotation (Task 3)**
- Marker-gene dot plot + Azimuth (built on Seurat's anchor/MNN label transfer).
- scANVI (semi-supervised scVI extension) is a valid alternative/complement to
  Azimuth for a bonus cross-method comparison.

**Composition (Task 4)**
- Distinct question from DE: proportion shift vs. per-cell expression shift.
  Don't conflate — bulk-style reasoning would confound these.

**DE (Task 5)**
- Use pseudobulk (sum counts per cell-type x sample) + DESeq2/edgeR, not per-cell
  testing directly — avoids pseudoreplication (cells within a sample aren't
  independent replicates).
- Only 1 donor — no donor-level replication. Flag this limitation explicitly,
  don't paper over it.
- Expect effects concentrated in monocytes/myeloid cells (phagocytic — actually take
  up particles), not lymphocytes. TNF/IL-6/IL-8 are literature-plausible candidates
  but NOT guaranteed — check pathway-level enrichment (NF-kB/cytokine signaling)
  rather than betting on one gene.

**Pathway enrichment (Task 5/6)**
- ORA (hypergeometric/Fisher's exact, on thresholded DE gene list) vs. GSEA (whole
  ranked gene list, catches distributed sub-threshold effects). Prefer GSEA when a
  single gene might not individually clear significance but a pathway module does.
- `gseapy` in Python; pull MSigDB Hallmark/KEGG/Reactome gene sets rather than
  hand-curating.

**Bonus analysis candidates — status: deferred until a full pipeline pass exists**
These are candidates, not committed work. Don't build any of these into the first
pass through Tasks 1-6. Revisit once the whole pipeline runs end-to-end on standard QC.
1. Mito%/complexity joint QC rework + per-cell-type mito% as an exposure readout
   (needs Task 3 done first, and supersedes the first-pass QC once done).
2. Leiden vs. Louvain robustness comparison.
3. scVI native DE module (if scVI used for Task 2).
4. Marker-gene retention / downsampled-reference accuracy check on annotation quality.
5. scANVI vs. Azimuth annotation comparison.
6. Gene signature scoring (`sc.tl.score_genes`, Tirosh et al. 2016 method — per-cell
   score against expression-matched control genes). Likely already covered in the
   Sanbomics tutorial — check there first. Two separate applications, not one:
   - Score MSigDB Hallmark `TNFA_SIGNALING_VIA_NFKB` per cell, per cell type/sample —
     complements Task 5/6's DE and pathway enrichment without betting on one gene.
   - Score the van den Brink et al. 2017 dissociation-stress signature (heat-shock +
     immediate-early genes) and test whether high-mito% cells also score high on it —
     this is the actual test that resolves item 1 (technical artifact vs. real biology),
     not just a related idea.

## Known pitfalls (from working through a practice notebook on similar data)
- Watch for QC thresholds that are technically present in code but numerically
  inert on the actual dataset (silently filter ~nothing) — always sanity check
  before/after cell counts, not just that the filter code ran.
- Watch for mismatched layers: plotting/testing code must reference the layer it
  claims to (raw counts vs. normlog vs. scaled) — mislabeled axes/layers are an easy,
  easy-to-miss bug.
- UMAP `min_dist` and Leiden `resolution` defaults differ for spatial vs. single-cell
  workflows — don't copy spatial-tutorial defaults blindly (not directly relevant here
  since this project isn't spatial, but a reminder to always state actual parameters
  used, not assumed defaults).

## Style
- Justify every threshold/parameter choice in comments or markdown cells, not just
  in the final report — this project's rubric explicitly asks for justification.
