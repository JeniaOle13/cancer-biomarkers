# Cancer-Immunotherapy-Microbiome-Analysis  
**Code repository for the manuscript:** *Ecological and functional stratification of the stool microbiome predicts response to immune checkpoint inhibitors across cancer types*

![](https://github.com/JeniaOle13/cancer-biomarkers/blob/main/figure/sample_map.jpg)
*Global distribution of collected samples (n=624) across countries. Circle size represents sample count per region, while color indicates the proportion of patients R - responsive (or NR - non-responsive) to immunotherapy*

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

🔗 **Access the full study report:** [https://jeniaole13.github.io/cancer-biomarkers/](https://jeniaole13.github.io/cancer-biomarkers/)

## Overview  
This project presents a complete bioinformatics pipeline for identifying gut microbial biomarkers associated with response to immunotherapy in cancer patients. Based on 976 stool metagenomes from 14 public studies, we constructed a catalog of 3,816 non-redundant metagenome-assembled genomes (MAGs), which were used for taxonomic and functional profiling. This study aims to discover gut microbial signatures predictive of immunotherapy success by identifying specific operational genomic units (OGUs) linked to patient response (R) or non-response (NR). We further investigate how microbial prevalence influences outcomes, characterize the functional profiles (CAZy, KEGG) of these marker taxa.

## Results
Analysis of 3,816 metagenome-assembled genomes identified 350 microbial markers, with 149 enriched in responders and 201 in non-responders. Fusicatenibacter saccharivorans emerged as a consistent pan-cancer biomarker for positive response, while non-responders showed enrichment of Proteobacteria and taxa originating from oral/food sources. Functionally, responder-associated microbes possessed expanded CAZy repertoires for complex carbohydrate and mucin degradation.

## Repository Structure
```
data/                       # project data
├── songbird/               # Songbird analysis results for each dataset
markers-report_files/       # Files for Quarto report
README.md                   # Repository description file
markers-report.html         # Quarto HTML file
markers-report.qmd          # Quarto qmd file
```
## Data Availability
In this study, we used open access data from the NCBI/EBI Sequence Read Archives, identified by the following BioProjects accession numbers: [PRJNA397906](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA397906), [PRJEB22893](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB22893), [PRJNA399742](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA399742/), [PRJNA770295](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA770295/), [PRJEB43119](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB43119/), [PRJNA762360](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA762360/), [PRJNA1011235](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1011235/), [PRJNA615114](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA615114/), [PRJNA866654](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA866654/), [PRJNA494824](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA494824/), [PRJEB49516](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB49516/). OGUs catalog assembly pipeline was described in https://github.com/JeniaOle13/Cancer_MAGs. All initial MAGs sequences were deposited in NCBI under accession [PRJNA1196825](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1196825/).

## Studies links:
1) [Frankel, Arthur E., et al. "Metagenomic shotgun sequencing and unbiased metabolomic profiling identify specific human gut microbiota and metabolites associated with immune checkpoint therapy efficacy in melanoma patients." Neoplasia 19.10 (2017): 848-855.](https://www.sciencedirect.com/science/article/pii/S1476558617302385)
2) [Peng, Zhi, et al. "The gut microbiome is associated with clinical response to anti–PD-1/PD-L1 immunotherapy in gastrointestinal cancer." Cancer Immunology Research 8.10 (2020): 1251-1261.](https://aacrjournals.org/cancerimmunolres/article/8/10/1251/466881/The-Gut-Microbiome-Is-Associated-with-Clinical)
3) [Gopalakrishnan, Vancheswaran, et al. "Gut microbiome modulates response to anti–PD-1 immunotherapy in melanoma patients." Science 359.6371 (2018): 97-103.](https://www.science.org/doi/10.1126/science.aan4236)
4) [Matson, Vyara, et al. "The commensal microbiome is associated with anti–PD-1 efficacy in metastatic melanoma patients." Science 359.6371 (2018): 104-108.](https://www.science.org/doi/10.1126/science.aao3290)
5) [Spencer, Christine N., et al. "Dietary fiber and probiotics influence the gut microbiome and melanoma immunotherapy response." Science 374.6575 (2021): 1632-1640.](https://www.science.org/doi/10.1126/science.aaz7015?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed)
6) [Lee, Karla A., et al. "Cross-cohort gut microbiome associations with immune checkpoint inhibitor response in advanced melanoma." Nature Medicine 28.3 (2022): 535-544.](https://www.nature.com/articles/s41591-022-01695-5)
7) [McCulloch, John A., et al. "Intestinal microbiota signatures of clinical response and immune-related adverse events in melanoma patients treated with anti-PD-1." Nature Medicine 28.3 (2022): 545-556.](https://www.nature.com/articles/s41591-022-01698-2)
8) [Tsakmaklis, Anastasia, et al. "TIGIT+ NK cells in combination with specific gut microbiota features predict response to checkpoint inhibitor therapy in melanoma patients." BMC Cancer 23.1 (2023): 1160.](https://bmccancer.biomedcentral.com/articles/10.1186/s12885-023-11551-5)
9) [Liu, Ben, et al. "Exploring gut microbiome in predicting the efficacy of immunotherapy in non-small cell lung cancer." Cancers 14.21 (2022): 5401.](https://www.mdpi.com/2072-6694/14/21/5401)
10) [Heshiki, Yoshitaro, et al. "Predictable modulation of cancer treatment outcomes by the gut microbiota." Microbiome 8 (2020): 1-14.](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00811-2)
11) [Gunjur, Ashray, et al. "A gut microbial signature for combination immune checkpoint blockade across cancer types." Nature Medicine 30.3 (2024): 797-809.](https://www.nature.com/articles/s41591-024-02823-z)

## Authors:
[Vera Kanaeva](https://scholar.google.ru/citations?hl=ru&user=Ie7RMLAAAAAJ) (Lopukhin FRCC PCM, MIPT).

[Evgenii Olekhnovich](https://scholar.google.ru/citations?user=RA9ItlsAAAAJ&hl=ru) (Lopukhin FRCC PCM).

## Publication
This repository contains the complete analytical code for the research paper:

- **Orletskaia et al.** *Ecological and functional stratification of the stool microbiome predicts response to immune checkpoint inhibitors across cancer types.* bioRxiv (2025).  
  DOI: [https://doi.org/10.1101/2025.05.07.652660v2](https://www.biorxiv.org/content/10.1101/2025.05.07.652660v2)

## Related Publications

- **Olekhnovich et al.** *Consistent stool metagenomic biomarkers associated with the response to melanoma immunotherapy.* mSystems 8.2 (2023).
  DOI: [https://doi.org/10.1128/msystems.01023-22](https://doi.org/10.1128/msystems.01023-22)

- **Zakharevich et al.** *Systemic metabolic depletion of gut microbiome undermines responsiveness to melanoma immunotherapy.* Life Science Alliance 7.5 (2024).  
  DOI: [https://doi.org/10.26508/lsa.202302480](https://doi.org/10.26508/lsa.202302480)

## Findings
Financial support for this study was provided by the Russian Science Foundation under the grant #22-75-10029 (https://rscf.ru/project/22-75-10029/).

## License
This project is licensed under the MIT License.
