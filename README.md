# NHHR and cardiovascular mortality in US adults with diabetes

**Authors**: [Ruil]
**Journal**: BMC Endocrine Disorders (submitted)
**DOI**: [Zenodo 将自动生成]

## Description
This repository contains the R code used to reproduce the analysis for the paper:
"Association between the non-high-density lipoprotein cholesterol to high-density lipoprotein cholesterol ratio and cardiovascular mortality among US adults with diabetes: a prospective cohort study of NHANES 1999–2018".

## Data source
Data are publicly available from:
- NHANES: https://www.cdc.gov/nchs/nhanes/
- NHANES Linked Mortality Files: https://www.cdc.gov/nchs/data-linkage/mortality.htm

## Repository contents
- `final_analysis.R` - Complete R script for data processing, statistical analysis, and figure generation
- `Figure_1_Flowchart.tiff` - Study population selection flow diagram
- `Figure_2_RCS.tiff` - Restricted cubic spline analysis
- `Figure_3_KM.tiff` - Kaplan-Meier survival curves

## Requirements
- R >= 4.1.0
- Required packages: nhanesdata, dplyr, survey, survival, ggplot2, rms, cmprsk, survminer, foreign, readr

## How to run
1. Clone this repository or download the ZIP file
2. Open `final_analysis.R` in RStudio
3. Run the script sequentially

The script will automatically download NHANES data, perform all analyses, and generate the figures reported in the manuscript.

## License
MIT License – you are free to reuse, modify, and distribute this code with attribution.
