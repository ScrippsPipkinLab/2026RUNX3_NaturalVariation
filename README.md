# AlbaoRunx3Manuscript

Code and analysis assets for the manuscript *Variation in RUNX protein transcriptional activity determines T cell exhaustion and the pattern of memory CD8 T cell formation* (Albao *et al.*). The repository collects the single-cell RNA-seq, ATAC-seq, CUT&RUN, genome track and flow cytometry analyses used throughout the study, together with figure exports and helper scripts.

## Repository layout

- `single_cell/` — Scanpy-based notebooks and outputs for preprocessing, clustering, differential expression, pathway enrichment, signature scoring, hashtag demultiplexing and RNA velocity analyses (RPE, MRE and WTE experiments).
- `atac/` — ATAC-seq peak filtering, differential accessibility, peak set enrichment analysis (PSEA), and motif analyses implemented in Jupyter notebooks and exported HTML reports.
- `cutnrun/` — RUNX3 and RUNX1 CUT&RUN analyses: peak exploration and cleanup, motif enrichment, peak annotation, functional enrichment, replicate/PCA analyses, Type 1/Type 2 peak clustering, integration with ATAC-seq, and comparison with published RUNX3 ChIP-seq.
- `flow/` — R notebooks for statistical analysis of flow cytometry data, organized by main figure (`Fig01`–`Fig06`), plus shared helpers in `flow/scripts/`.
- `tracks/` — R scripts for preprocessing and plotting genome tracks (CUT&RUN/ATAC-seq loci, chromatin-associated RNA-seq, and retroviral pseudo-chromosomes).
- `csv/` — Compressed intermediate tables for clustering, differential expression, GSEA outputs, and figure-ready summaries.
- `figures/` — Generated figure assets grouped by panel.
- `reanalysis/` — Supplemental GSEA notebook for comparative datasets.
- `scripts/` — Utility scripts for SEA plotting, volcano plots, and notebook conversion.

## Analysis overview

### Single-cell RNA sequencing

- P14 CD8 T cells were hash tagged with BioLegend TotalSeq A antibodies (A0301–A0309) before pooling and loading 100,000 multiplexed cells across two 10x Genomics Single Cell 3' v3.1 GEMs.
- Hash tag oligo (HTO) libraries were amplified from cDNA with the BioLegend HTO additive primer, while gene expression (GEX) libraries followed the standard 10x protocol. Sequencing depth targeted ~35,000 reads per cell for GEX and ~300 reads per cell for HTO. The RPE/MRE experiment was sequenced on an Illumina NextSeq 2000 and the WTE experiment on an Illumina NextSeq 500.
- CellRanger v7.1.0 with the 2020-A mm10 reference generated count matrices. Quality control in Scanpy v1.9.5 included SoupX ambient RNA correction, Scrublet and scDblFinder doublet removal, and filtering on mitochondrial content (>5%), UMI count (<3,000), and detected genes (<1,250). Cell cycle effects were regressed using S and G2/M gene lists, and Leiden clustering at resolution 1.0 (14–16 clusters) guided downstream analyses.
- Differential expression relied on the Scanpy Wilcoxon rank-sum test (v1.9.6), with GSEA run on Wilcoxon statistics via clusterProfiler v4.14.0 and DOSE v4.0.0 (R v4.4.2). Transcriptome correlations used Spearman coefficients on mean-normalized counts of the top 2,000 variable genes. RNA velocity inputs were produced with velocyto v0.17 and analyzed with scVelo v0.3.1.

### ATAC-seq

- Nuclei from sorted P14 CD8 T cells were permeabilized for in situ Tn5 transposition (Nextera DNA Library Preparation kit), titrated to 1.25 μL enzyme per 5×10⁴ nuclei in a 50 μL reaction. Libraries were PCR-amplified to optimal cycles, quality-checked on an Agilent TapeStation (targeting a 2:1 170–280 bp to 280–390 bp fragment ratio), and sequenced on an Illumina NextSeq 2000 (paired-end 61×2, ~3×10⁷ reads per sample).
- Processing and peak calling used the nf-core/atacseq v2.1.2 pipeline (Nextflow 24.04.2) against mm10 in broad-peak mode. Public datasets (GSE111149, GSE88987, GSE144383, GSE131871, GSE213041) were merged with study data to define consensus peaks. A peak was considered observed in a condition if called in ≥75% of that condition's replicates; peaks observed in at least one condition from at least three datasets (excluding chrY) yielded 47,192 shared peaks from an initial 285,385.
- Differential accessibility employed DESeq2 v1.46.0 (each dataset under its own design), PSEA used clusterProfiler/DOSE on the Wald statistic (10,000 permutations, minimum set size 5, BH-corrected), and motif analysis relied on HOMER v5.1. Peak overlaps with ChIP-seq datasets (RUNX3, GSE50131; TBET, GSE72408) were computed with bedtools v2.31.1 (`intersect -f 0.5 -F 0.5`), with accessibility profiles generated via deepTools v3.5.6.

### CUT&RUN

- Native CUT&RUN for RUNX3 and RUNX1 was performed on 7.5×10⁴–1.0×10⁵ nuclei per sample with the EpiCypher CUTANA CUT&RUN v6 kit, using buffers supplemented with 0.0625% NP-40 and protease inhibitors. Two systems were profiled: day-five P14 effectors transduced with sh*Runx3* or sh*Cd19*, and day-eight unperturbed *Tcf7*⁺, α, β and γ effector subsets.
- Libraries were prepared with the NEBNext Ultra II DNA Library Prep Kit and sequenced on an Illumina NextSeq 2000 (paired-end 2×61, ~1×10⁷ reads per sample).
- Reads were processed with a modified nf-core/cutandrun v3.2.2 pipeline ([ScrippsPipkinLab/cutandrun_fork](https://github.com/ScrippsPipkinLab/cutandrun_fork)) that adds replicate merging and cross-condition consensus peak generation. Reads were aligned to mm10, signal was CPM-normalized, and peaks were called with SEACR v1.3 in stringent mode.
- RUNX3 and RUNX1 peaks were combined and k-means clustered with deepTools v3.5.6 (`plotHeatmap`, k = 2…7). Sample-restricted clusters lacking RUNT motif enrichment (HOMER v5.1) were removed as spurious, retaining 4,172 unified peaks. These were partitioned into Type 1 peaks (RUNX1 binding increases upon RUNX3 suppression) and Type 2 peaks (RUNX1 binding unchanged).
- Downstream analyses include HOMER known-motif enrichment against GC-matched background, peak annotation with ChIPseeker, functional enrichment with rGREAT v2.8.0 against custom gene sets, PCA of normalized signal, deepTools signal profiles at signature-associated loci, integration with ATAC-seq differential accessibility, and overlap with published RUNX3 ChIP-seq peaks (GSE50131).

### Flow cytometry

- Compensation and gating were performed in FlowJo 10.10. Workspaces were imported into R v4.4.2 with a development version of [fcexpr](https://github.com/Close-your-eyes/fcexpr) (frozen fork at *dsalbao/fcexpr*), with population statistics cached by `flow/scripts/flowProcessing.R`.
- Notebooks in `flow/FigXX/` reproduce the statistics and plots for each figure, including RUNX3 dosage phenotyping and cell numbers (Fig 1), fate mapping and rechallenge of RUNX3-perturbed and natural α/β/γ effectors (Fig 2, Extended Fig 4), RUNX/CBFβ shRNAmir phenotyping and CXCR6 staining (Fig 3, Extended Fig 5), RUNX3 domain-mutant phenotyping, accumulation and fate mapping (Fig 4, Extended Fig 7), and RUNX3/TOX epistasis, LCMV Armstrong vs Clone 13 phenotyping, and restimulation (Fig 5–6, Extended Figs 8–9).
- Statistical tests are declared in each figure legend (two-tailed Welch or paired T-tests); where multiple comparisons are made, P values are corrected with the Benjamini–Hochberg procedure.

### Gene and peak resources

- CD8 T cell-specific gene lists follow prior definitions, supplemented by re-analyses of public bulk RNA-seq datasets (Tcf7-reporter, ex-KLRG1, Tox genotypes, terminal-TEM). Reads were mapped with Salmon v1.10.2 (mm39, Ensembl 110) using dataset-specific k-mer parameters, with differential expression performed via DESeq2 v1.46.0.
- ATAC-seq peak lists were generated ad hoc from the consensus procedure above to support PSEA.

### Computational reproducibility

Public pipelines:

- nf-core/atacseq v2.1.2 (commit `1a1dbe52ffbd82256c941a032b0e22abbd925b8a`)
- nf-core/nascent v2.2.0 (commit `02bbefb701598dc96d29da23d90d60a9ffadd18d`)
- Modified nf-core/cutandrun v3.2.2: [ScrippsPipkinLab/cutandrun_fork](https://github.com/ScrippsPipkinLab/cutandrun_fork) (commit `634daa0972310f6e899ba47bb44efb5fcd0bc500`)

Containerized tools:

- Public images:
  - CellRanger v7.1.0 (`docker.io/litd/docker-cellranger:v7.1.0`) for alignment and count matrix generation.
  - bedtools v2.31.1 (`docker.io/staphb/bedtools:2.31.1`) for peak overlaps.
  - deepTools v3.5.6 (`quay.io/biocontainers/deeptools:3.5.6--pyhdfd78af_0`) for signal profiles and heatmaps.
  - Salmon v1.10.2 (`docker.io/combinelab/salmon:1.10.2`) for transcript quantification.
- Custom DockerHub images:
  - HOMER v5.1 (`docker.io/pipkinlab/homer:5.1`) for motif analysis.
  - Scanpy v1.9.5 (`docker.io/pipkinlab/scanpy:1.9.5`) for core scRNA-seq QC and clustering.
  - Scanpy v1.9.6 with scVelo (`docker.io/pipkinlab/scanpy:1.9.6`) for differential expression and RNA velocity.
  - R 4.4.2 with DESeq2 v1.46.0 and plotting stack (`docker.io/pipkinlab/r-deseq2-plotting:4.4.2`) for differential accessibility and visualization.

## Data availability

scRNA-seq, ATAC-seq, and CUT&RUN data are available in the Gene Expression Omnibus under SuperSeries accession [GSE348910](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE348910). Gene and peak sets are available at [ScrippsPipkinLab/CommonGeneSets](https://github.com/ScrippsPipkinLab/CommonGeneSets).
