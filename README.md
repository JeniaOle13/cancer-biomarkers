# Cancer-Immunotherapy-Microbiome-Analysis  
**Code repository for the manuscript:** *"Non-Intestinal Microbial Signatures in Stool as Predictors of Cancer Immunotherapy Outcome"*

![](https://github.com/JeniaOle13/cancer-biomarkers/blob/main/data/sample_map.jpg)
*Global distribution of collected samples (n=624) across countries. Circle size represents sample count per region, while color indicates the proportion of patients R - responsive (or NR - non-responsive) to immunotherapy*

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

🔗 **Access the full study report:** [https://jeniaole13.github.io/cancer-biomarkers/](https://jeniaole13.github.io/cancer-biomarkers/)

## Overview  
This repository contains the computational workflows and scripts used to analyze 624 stool metagenomes from 11 cohorts to identify microbial biomarkers of cancer immunotherapy response. Key findings include:  
- **Blautia species** as keystone taxa linked to immunotherapy success.  
- **Oral/food-derived microbes** (e.g., *Veillonella*, *Klebsiella*) as predictors of non-response.

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
In this study, we used open access data from the NCBI/EBI Sequence Read Archives, identified by the following BioProjects accession numbers: PRJNA397906, PRJEB22893, PRJNA399742, PRJNA770295, PRJEB43119, PRJNA762360, PRJNA1011235, PRJNA615114, PRJNA866654, PRJNA494824, PRJEB49516. OGUs catalog assembly pipeline was described in https://github.com/JeniaOle13/Cancer_MAGs. All initial MAGs sequences were deposited in NCBI under accession PRJNA1196825.

## Authors:
[Vera Kanaeva](https://scholar.google.ru/citations?hl=ru&user=Ie7RMLAAAAAJ) (Lopukhin FRCC PCM, MIPT).

[Evgenii Olekhnovich](https://scholar.google.ru/citations?user=RA9ItlsAAAAJ&hl=ru) (Lopukhin FRCC PCM).

## Findings
Financial support for this study was provided by the Russian Science Foundation under the grant #22-75-10029 (https://rscf.ru/project/22-75-10029/).

## License
This project is licensed under the MIT License.
