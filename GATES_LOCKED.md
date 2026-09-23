# GATES_LOCKED - amr-genome-early-warning

Locked 2026-09-24 (IST) BEFORE any data processing or results commit. These gates are not edited after this commit. Any deviation is logged in DEVIATIONS.md with reason.

## Question
When AMR is predicted from genomes, do published accuracy numbers (random splits) hold for isolates collected AFTER the training period? Which failures come from resistance determinants the model never saw, and can a structural/allele-novelty signal flag those isolates early?

## Prior art (checked 2026-09-23)
- Temporal validation of genome AMR prediction exists for single species (e.g. NasirNesirli/kleb-amr-project, pre-2023 train / 2023-24 test, K. pneumoniae).
- AMRscope (bioRxiv 10.1101/2025.09.12.672331): PLM risk scores for novel missense variants; splits by gene/organism/homology, not by date.
- amr.watch (PLOS GPH 10.1371/journal.pgph.0005256): genomic AMR trend monitoring, no phenotype prediction benchmark.
- ResFinder 4.0 / AMRFinderPlus: rule-based genotype->phenotype, validated on fixed collections.
Gap targeted: cross-species date-split benchmark on NCBI Pathogen Detection AST isolates that decomposes date-split errors into (a) unseen determinant, (b) known determinant/wrong call, (c) no determinant; plus an early-warning score for (a).

## Data (real public only)
- NCBI Pathogen Detection AMR metadata: Klebsiella PDG000000012.2532, E. coli/Shigella PDG000000004.6316 (sha256 in data/SOURCES.md).
- AMRFinderPlus DB 2026-08-07.1 ReferenceGeneCatalog + allele_counts_by_year/bla.tsv.
- Structures from RCSB PDB / AlphaFold DB for beta-lactamase families in the structural layer.

## Pre-registered analysis
- Drugs: ceftazidime, meropenem, ciprofloxacin, gentamicin, trimethoprim-sulfamethoxazole (most-tested classes). I counted as R only in a sensitivity analysis; primary = R vs S, I excluded.
- Features: presence/absence of AMRFinderPlus genotypes (AMR_genotypes column).
- Models: (M0) rule baseline = any AMRFinderPlus determinant of the drug's class -> R; (M1) L2 logistic regression; (M2) random forest (500 trees). Seed 42.
- Splits: RANDOM = stratified 80/20; DATE = train collection year <= 2019, test >= 2020. Isolates without a year excluded.
- Primary metric: balanced accuracy (BA). Secondary: MCC, sensitivity for R.

## Gates
- G0 Data: >= 3,000 isolates with year + AST across both species; per-drug test sets >= 100 with >= 20 R. Fail -> drop that drug (logged), not the gate.
- G1 Temporal drop: date-split BA vs random-split BA per drug x model, 1,000 bootstrap CI. Claim "temporal degradation" only if mean drop >= 0.02 AND CI excludes 0. Otherwise report as a NEGATIVE: "no measurable temporal drop".
- G2 Error decomposition: every date-split false negative assigned to (a) determinant present but absent from training vocabulary, (b) known determinant, (c) no determinant. Reported as-is.
- G3 Early-warning score: novelty score = count of determinants in the isolate unseen before the cutoff, plus AMRFinderPlus allele first-year > cutoff. Pass only if AUROC for flagging M2 false negatives >= 0.65 with CI lower bound > 0.5. Fail -> reported as negative, not re-tuned.
- G4 Structural layer: for bla alleles first seen after cutoff (bla.tsv / catalog), map substitutions vs family reference onto PDB/AlphaFold structure; test whether active-site-proximal (<= 8 A from catalytic Ser/Zn site) substitutions are enriched among alleles in resistant-misclassified isolates (Fisher exact, alpha 0.05). Descriptive if n < 10.
- G5 Tool: CLI that takes an AMRFinderPlus genotype list + year and returns predicted phenotype plus early-warning flag; tested on held-out date-split isolates.

## Anti-fishing
No new drugs, cutoffs, models or thresholds after results are seen. Negatives stay in the paper.
