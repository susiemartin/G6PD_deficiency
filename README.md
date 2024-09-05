# Undiagnosed G6PD deficiency is common in Black individuals in the UK, and contributes to health inequalities in type 2 diabetes diagnosis

Electronic health record (EHR) code lists, R scripts and PLINK commands required for analysis of this project, as well as model diagnostic plots referenced in the **Supplementary Appendix**.

## EHR code lists for clinical phenotype defintions
EHR code lists for all conditions used in the analysis are given in read v2, read v3, ICD-9, ICD-10 and SNOMED format, as well as self-reported definitions from the UK Biobank 20002 and 20003 fields. Details of the UK Biobank fields, assessment centre visits and values used to define each phenotype are given in **UKB_phenotype_definitions.csv**.

Previously published EHR code lists used here with references:
- 20002_T2D.txt, 20002_diabetes.txt, 20003_diabetes_meds.txt, read_2_diabetes_qof.txt, read_2_diagnosis.txt, read_2_drugs.txt, read_2_exclusion.txt, read_3_diabetes_qof.txt, read_3_diagnosis.txt, read_3_exclusion.txt, SNOMED_diabetes.csv - from Young et al^1.

- MULTIPLY_Snomed_T2D.csv

## R scripts and PLINK commands

PLINK commands to extract whole exome sequencing (WES) data for UK Biobank (on DNAnexus) and Genes & Health (on TRE):


Below R scripts are split in data preparation (**prep_\*.R**) and analysis (**analysis_*.R**).

### UK Biobank
Below R scripts to be used with UK Biobank data on DNAnexus.

Data preparation scripts are to be run in the following order:
- **prep_UKB_genetic_data_step1.R**: Script to prepare genetic data (after applying above PLINK commands)
- **prep_UKB_phenotype_data_step2.R**: Script to prepare phenotype data - uses EHR codes and defines all clinical phenotypes used in analyses
- **prep_UKB_data_cleaning_step3.R**: combines genetic and phenotype data and applies necessary quality control to this data.

Scripts to conduct statistical analyses:
- **analysis_UKB_prevalence.R**: Calculates prevalence of G6PD variants and G6PD deficiency diagnoses
- **analysis_UKB_HbA1c.R**: Comparison of HbA1c and random glucose across genotypes in individuals without diabetes or pregnancy
- **analysis_UKB_age_diagnosis.R**: Comparison of age of type 2 diabetes diagnosis across genotpyes in individuals diagnosed since 2011.

### Genes & Health
Below R scripts to be used with Genes & Health data on TRE. Data preparation scripts will be available soon.

Scripts to conduct statistical analyses:
- **analysis_GnH_prevalence.R**: Calculates prevalence of G6PD variants and G6PD deficiency diagnoses
- **analysis_GnH_HbA1c.R**: Comparison of HbA1c and random glucose across genotypes in individuals without diabetes or pregnancy
- **analysis_GnH_age_diagnosis.R**: Comparison of age of type 2 diabetes diagnosis across genotpyes in individuals diagnosed since 2011.

## UK Biobank version date and fields extracted
The UK Biobank data used for this analysis (stored on DNAnexus) was last updated on March 1st 2024 (total N=502,185). The UK Biobank fields extracted from the cohort dataset (using the 'Table Exporter' tool) into *.csv files that are used in the R scripts are given in **UKB_table_exporter_fields.xlsx**.

## Model diagnostic plots
Model diagnostic plots for the linear regression model used in the analysis of age of type 2 diagnosis and described in the **Supplementary Appendix**. Names of plots reflect the cohort and ethnicity grouping that the model was applied to.

## References
1. Young KG, McGovern AP, Barroso I, et al. The impact of population-level HbA(1c) screening on reducing diabetes diagnostic delay in middle-aged adults: a UK Biobank analysis. Diabetologia 2023;66:300-9.
2. Minassian C, Williams R, Meeraus WH, Smeeth L, Campbell OMR, Thomas SL. Methods to generate and validate a Pregnancy Register in the UK Clinical Practice Research Datalink primary care database. Pharmacoepidemiol Drug Saf 2019;28:923-33.
