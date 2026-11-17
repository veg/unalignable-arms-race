# Evolutionary Arms Races Across Unalignable Mammalian Host–Virus Interfaces

Sergei L. Kosakovsky Pond, Darren P. Martin, Steven Weaver, and Anton Nekrutenko

Companion study demonstrating unified statistical homology alignment, reticulate Ohmic Wheatstone gene-conversion routing, pseudogene admission, and regularized biophysical selection (`WALRUS` and `S-PRIME` in [HyPhy 3](https://github.com/veg/hyphy-3)) across **25 pre-registered Tier-1 mammalian loci** with experimental PDB co-structures and **7,045 canonical Tier-2 TOGA2 mammalian orthogroups**.

## Repository Structure

- [`manuscript/`](manuscript/) — *Science*-format manuscript ([`main.tex`](manuscript/main.tex), [`main.pdf`](manuscript/main.pdf)), Supplementary Materials ([`supplement.tex`](manuscript/supplement.tex), [`supplement.pdf`](manuscript/supplement.pdf)), merged [`science_combined.pdf`](manuscript/science_combined.pdf), BibTeX bibliography ([`references.bib`](manuscript/references.bib)), vector figures ([`manuscript/figures/`](manuscript/figures/)), generated LaTeX tables ([`table1_summary.tex`](manuscript/table1_summary.tex), [`tables_supplement_generated.tex`](manuscript/tables_supplement_generated.tex)), and preserved bioRxiv v1 draft (`main_biorxiv_v1.*`). Tracked directly against both [veg/unalignable-arms-race-paper](https://github.com/veg/unalignable-arms-race-paper) (`origin`) and Overleaf (`overleaf`), and mirrored in this repository.
- [`scripts/`](scripts/) — End-to-end curation, execution, synthesis, figure/table generation, and Overleaf synchronization scripts:
  - [`synthesize_full_study.py`](scripts/synthesize_full_study.py) — Aggregates WALRUS, S-PRIME, PDB contact, legacy baseline, and Tier-2 TOGA2 statistics into JSON summaries.
  - [`generate_science_figures.py`](scripts/generate_science_figures.py) — Generates all 6 Type-42 vector PDF and 300-dpi PNG figures (`fig1`–`fig4`, `fig_s1`–`fig_s2`).
  - [`generate_science_tables.py`](scripts/generate_science_tables.py) — Generates Table 1 and Supplementary Tables S1–S4 from the verified JSON syntheses.
  - [`sync_overleaf.sh`](scripts/sync_overleaf.sh) — Keeps `manuscript/` in sync between Overleaf and this GitHub repository.
- [`data/`](data/) — Curated unaligned mammalian CDS FASTA inputs (`*_unaligned.fa`, `*_manifest.json`), pre-registered PDB heavy-atom structural distance ground truth (`ground_truth/tier1_25gene_ground_truth.json`), TOGA2 canonical summaries (`toga2_completed_summaries.jsonl.gz`), and assembly metadata (`assemblies_and_species.tsv`).
- [`results/`](results/) — Executed `WALRUS` outputs (`.walrus.json`, interactive `.walrus.html` atlases, `.walrus.msa.fa`, `.walrus.sto`, `.walrus.gfa`), `S-PRIME` JSON outputs (`prime/`), legacy MAFFT+MEME baselines (`legacy_foil/`), and master JSON syntheses (`tier1_25gene_master_synthesis.json`, `tier2_toga2_vip_summary.json`).
- [`figures/`](figures/) — High-resolution PNG previews and archival plot outputs.

## Reproducing Figures, Tables, and Manuscript PDFs

From the repository root:

```bash
python3 scripts/synthesize_full_study.py
python3 scripts/generate_science_figures.py
python3 scripts/generate_science_tables.py

cd manuscript
latexmk -pdf -interaction=nonstopmode main.tex
latexmk -pdf -interaction=nonstopmode supplement.tex
```

## Syncing Between Overleaf and GitHub

`manuscript/` has its own dual-remote repository (`origin` = [veg/unalignable-arms-race-paper](https://github.com/veg/unalignable-arms-race-paper), `overleaf` = Overleaf) while all files inside `manuscript/` are also mirrored directly in `veg/unalignable-arms-race`:

```bash
# Pull updates from Overleaf + veg/unalignable-arms-race-paper and mirror into root index
./scripts/sync_overleaf.sh

# Commit & push manuscript changes to Overleaf, veg/unalignable-arms-race-paper, and veg/unalignable-arms-race
./scripts/sync_overleaf.sh --push
```

> **Note for contributors:** Do not run `git add .` or `git add manuscript` from the repository root while `manuscript/.git` exists; use `./scripts/sync_overleaf.sh` (or stage explicit paths inside `manuscript/`) so Git tracks files rather than creating a nested submodule gitlink.
