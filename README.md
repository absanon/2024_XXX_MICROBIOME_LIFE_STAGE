# Convergent enrichment of Gammaproteobacteria along *Aedes aegypti* development across different breeding sites 

## Abstract
**Background:** *Aedes aegypti* mosquitoes are the main vector of pathogens like dengue virus and chikungunya virus. The immature life stages of mosquitoes share the same habitat with a variety of microorganisms in aquatic environments. To better understand the microbial diversity in field-derived *Ae. aegypti*, we analysed simultaneously collected larvae, pupae, and freshly emerged adults from Burkina Faso together with their breeding water via 16S rRNA gene sequencing.

**Results:** We observed a decrease in bacterial diversity and load across the mosquito life stages. At the phylum level, a strong increase in relative abundance of Proteobacteria was found along the mosquito stages. The same 40 amplicon sequence variants were consistently found as most abundant in the adults, regardless of the sample collection site, and all belonged to the Gammaproteobacteria. Our data suggest that these bacteria were not randomly derived by chance from the environment in the mosquito but rather deposited by a female mosquito during oviposition, a transmission route recently coined as “diagonal transmission”. Indeed, our results indicated that there is a selection of Gammaproteobacteria from the breeding water and that these bacterial members are further maintained from larvae to adults. 

**Conclusion:** This study provided new data on the microbiome composition of field-collected *Ae. aegypti*, contributing to an enhanced understanding of the origin and colonization route of the mosquito microbiome, potentially via a diagonal transmission route.

---

## Repository overview

### Main folders
- `data/`  
  Input tables used in downstream analyses (e.g. metadata and 16S qPCR results).
- `figures/`  
  Output figures (PDFs) written by the analysis scripts.
- `public_data/`  
  Helper scripts + metadata to reproduce the public datasets re-analyses (Bennett/Hery/Hernandez/Rodpai), plus a download script.
- `renv/`  
  Reproducible R environment (package versions are managed via **renv**).

### Key scripts (project root)
- `dada2.R`  
  DADA2 processing workflow (filter/trim → denoise → merge → chimera removal → taxonomy assignment → phyloseq object creation).  
  *Note:* this depends on having the raw FASTQ files and reference databases available locally.
- `alpha_diversity.R`  
  Rarefaction-based alpha diversity (Observed/Shannon/Simpson) + plots.
- `beta_diversity.R`  
  Beta diversity analysis (avgdist), PCoA/NMDS, and dbRDA workflow (uses `custom_rldbrda.R`).
- `relative_abundance.R`  
  Relative abundance barplots, alluvial plots, Venn/Upset comparisons, and absolute abundance using 16S qPCR.
- `bf_map.R`  
  Map figure for Burkina Faso / Ouagadougou (requires GIS inputs; paths may need editing).
- `custom_rldbrda.R`  
  Helper functions for dbRDA feature selection / effect sizes.

### Public dataset re-analyses
- `public_data/public_data.R`  
  Combines/compares results across the re-analysed public datasets (alpha diversity etc.).
- `public_data/relative_abundance_*.R`  
  One script per study (e.g. `relative_abundance_bennett.R`, `relative_abundance_hery.R`, …) to produce relative abundance + alluvial plots and save intermediate `.RData` workspaces.
- `public_data/download_public_data.sh`  
  Downloads metadata from SRA and fetches FASTQs.

---

## How to reproduce the analysis

### 1) Get the project and restore the R environment (renv)
This repository uses **renv** to pin package versions.

1. Clone the repository and open the R project file: `2024_XXX_MICROBIOME_LIFE_STAGE.Rproj`
2. Restore packages:

```r
renv::restore()
```

If the renv autoloader is enabled, it should activate automatically when opening the project.

### 2) Reproduce the main manuscript figures (Burkina Faso field dataset)
Most downstream scripts load the precomputed workspace:

- `Sanon_16S_DADA2_data.RData`

```r
source("alpha_diversity.R")
source("beta_diversity.R")
source("relative_abundance.R")
```

Outputs are primarily written to:
- `figures/` (PDF figures)

Inputs expected in:
- `data/metadata.tsv`
- `data/16S_qPCR_results.csv`

### 3) (Optional) Re-run DADA2 from FASTQs
To re-run the full DADA2 pipeline you must provide:
- Raw FASTQ files (location configured in `dada2.R`)
- Reference taxonomy files (e.g. SILVA training set / species assignment; paths currently hard-coded)

After adjusting paths, run:

```r
source("dada2.R")
```

This should recreate the objects that are later saved/loaded (e.g. the phyloseq object) depending on your local setup.

### 4) Reproduce analyses for public datasets
The public dataset scripts are written to run on an HPC-style filesystem. You will typically need to edit the paths to match your machine/environment.

Workflow:
1. Download metadata + FASTQs (requires `esearch/efetch/xtract` from Entrez Direct and `fasterq-dump` from SRA toolkit):

```sh
bash public_data/download_public_data.sh
```

2. Run the per-study re-analysis scripts (after adapting paths), e.g.:

```r
source("public_data/relative_abundance_bennett.R")
source("public_data/relative_abundance_hery.R")
source("public_data/relative_abundance_hernandez.R")
source("public_data/relative_abundance_rodpai.R")
```

3. Run the cross-study summary:

```r
source("public_data/public_data.R")
```