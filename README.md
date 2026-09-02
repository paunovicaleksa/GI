# Nanoplastics scRNA-seq — Immune Response

Single-cell RNA-seq analysis of human PBMCs from **one donor**, exposed to
carboxylated polystyrene nanoparticles (PSNPs) of two sizes, asking how
**particle size** shapes the immune response.

| Sample | Exposure |
|---|---|
| `40nm` | 40 nm PSNPs |
| `200nm` | 200 nm PSNPs |
| `mixture` | 40 nm + 200 nm together |
| `control` | no exposure |

Full task list, methodological constraints, and rationale for every threshold
and parameter choice are in [`CLAUDE.md`](CLAUDE.md) and worked out in detail in
[`PLAN.md`](PLAN.md). The analysis itself is [`nanoplastics_pipeline.ipynb`](nanoplastics_pipeline.ipynb).

## Headline findings

- **Response is confined to monocytes.** An apparent inflammatory signature in
  every cell type turned out to be ambient RNA (soup), not signal — real once
  isolated to the phagocytic, particle-internalizing population.
- **Pathways, monocytes, all three exposures:** NF-κB / TNF-α / inflammatory-
  response signaling, up. `Phagosome` — the pathway most directly tied to
  particle uptake — runs consistently **down**.
- **Genes:** `MMP1`, `IL24`, `CSF3`, `CSF2` are the largest movers (log2FC
  near +6). `IL8`/`CXCL8` is the top hit by magnitude but is also the gene
  driving the ambient-RNA artifact in other cell types, so it's read
  monocyte-only.
- **By expression, the response is shared** across 40nm/200nm/mixture — same
  ~77 pathways, same direction, in all three. Corroborated independently by
  per-cell signature scoring (no pseudobulk step).
- **By composition, the response is size-specific:** monocyte share falls
  1.9pp at 40nm, rises 3.4pp at 200nm, flat in the mixture. Composition and
  expression genuinely disagree — reading only one axis gives half the
  answer.
- **Mixture:** additive in monocyte expression (single-size contrasts predict
  it, coefficients sum to 0.93, R²=0.72) but not in composition — no
  synergistic/antagonistic gene expression effect from combining sizes.
- **Robust:** the monocyte composition/direction result is unchanged across
  OS, Leiden partition, and with batch correction removed entirely.
- **Limitation:** single donor, one sample per condition — no biological
  replication, so every result above is descriptive, not a significance test.

## Presentation

Slides are in [`nanoplastics_presentation-1.pdf`](nanoplastics_presentation-1.pdf); a recording walking through them is on YouTube: https://youtu.be/Oc8qLShCUYM

## Running it

### 1. Get the data

The four `.h5ad` files (167–272 MB each) are **not in this repo** — they
exceed GitHub's 100 MB file limit and are git-ignored. Obtain them separately
and place them at the repo root, alongside the `*_CoDi_KLD.csv` files (already
included) — the notebook's Section 1 expects these exact filenames:

```
filtered_feature_bc_matrix.h5ad          # 40nm
filtered_feature_bc_matrix_Sample2.h5ad  # 200nm
filtered_feature_bc_matrix_Sample3.h5ad  # mixture
filtered_feature_bc_matrix_Sample4.h5ad  # control
```

### 2. Set up the environment

Python 3.12, on Windows or Linux (both are tested and give the same results
through Section 7; see Section 17 of the notebook for what differs after
that and why it doesn't affect the conclusions).

```bash
python -m venv .venv

# Linux
source .venv/bin/activate
# Windows (PowerShell)
.\.venv\Scripts\Activate.ps1

python -m pip install --upgrade pip setuptools wheel
python -m pip install -r requirements.txt
python -m ipykernel install --user --name=nanoplastics-scrnaseq \
  --display-name="Python (nanoplastics-scrnaseq)"
```

Install `requirements.txt` in one pass rather than adding packages
incrementally, so pip's resolver sees every version constraint at once.
A C++ compiler is required on both platforms (`scanorama`'s `annoy`
dependency builds from source).

### 3. Run the notebook

Launch JupyterLab, open `nanoplastics_pipeline.ipynb`, select the
**Python (nanoplastics-scrnaseq)** kernel, and run top to bottom:

```bash
jupyter lab
```

A full cold run takes roughly 20–25 minutes (Pan-Human Azimuth's first-use
model download is the slowest single step). Two checkpoint files are written
to `results/` (git-ignored) partway through so a kernel restart doesn't force
re-running the clustering and annotation sweeps from scratch.

## Repository contents

| File | What it is |
|---|---|
| `nanoplastics_pipeline.ipynb` | The project notebook — all 6 tasks plus 3 additional robustness/depth analyses |
| `nanoplastics_presentation-1.pdf` | Presentation slides |
| `CLAUDE.md` | Project spec: tasks, required methods, non-negotiable constraints |
| `PLAN.md` | Working log of decisions, rationale, and section-by-section outcomes |
| `single_cell_analysis_complete_class.ipynb` | Reference tutorial notebook (SanBomics) — read-only, not part of the analysis |
| `requirements.txt` | Direct dependencies, installable on Windows or Linux |
| `requirements-lock-*.txt` | Full pinned `pip freeze` per platform, for exact reproduction |
| `*_CoDi_KLD.csv` | Pre-existing per-cell cell-type labels shipped with the data, used only as an annotation cross-check |
