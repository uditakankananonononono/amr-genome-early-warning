# amr-genome-early-warning

Do genome-based antibiotic resistance predictions hold up over time? A date-split benchmark, error anatomy and structural early-warning test on 12,841 public *K. pneumoniae* and *E. coli* isolates (NCBI Pathogen Detection).

**Paper:** [paper/paper.pdf](paper/paper.pdf) (20 pages)

## Key results
- Balanced accuracy drops by **0.134** (95% CI 0.109-0.158) on isolates collected 2020+ compared with a random split; this holds in all 20 ML cells.
- Of 1,281 missed resistant isolates, only **3.2%** carry a gene the model never saw; **67.8%** carry known determinants (65.9% after removing intrinsic genes, post-hoc).
- Negative: novelty score does not flag misses (AUROC 0.52).
- Negative (exploratory): newer beta-lactamase alleles are not enriched for active-site-proximal substitutions (OR 0.80, p 0.54).
- Tool: `tool/amr_ew.py`, a predictor that reports honest date-split accuracy plus per-determinant early-warning flags and structural annotations.

## Layout
- `GATES_LOCKED.md` - pre-registered gates (standalone commit 9f770aa). `DEVIATIONS.md` - disclosures.
- `data/` - derived isolate table, AMRFinderPlus tables, PDB structures, `SOURCES.md` (URLs + sha256).
- `code/` - dataset build, analysis (G1-G3), structural layer (G4), figures, appendix, sensitivity.
- `results/` - gate reports and raw metrics. `tool/` - CLI + validation. `paper/` - source + PDF.
- `MANIFEST.sha256` - sha256 of every committed file.

## Reproduce
```
python3 code/build_dataset.py     # needs raw PDG files in ../data/raw (see data/SOURCES.md)
python3 code/analysis.py && python3 code/structural_layer.py && python3 code/g2_sensitivity.py
python3 code/figures.py && python3 code/appendix.py && python3 tool/validate_tool.py
```
AMRProt.fa (AMRFinderPlus DB 2026-08-07.1, sha256 41d5ebf4f807c9590f27de7dd23c9f037afb4d8efd613790c54847fb1fc2896b) is needed in ../data/raw for the structural layer.
