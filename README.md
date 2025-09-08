# Undiagnosed G6PD deficiency in Black and Asian individuals is prevalent and contributes to health inequalities in type 2 diabetes diagnosis and complications

This repository contains the electronic health record (EHR) code lists, R scripts and PLINK commands required for analysis for Martin et al. (2025)[^1], as well as model diagnostic plots for the linear regression models detailed in the **Supplemental Material**.

[^1]: Martin, S, Samuel, M, Stow, D, Ridsdale, AM, Chen, J, Young, KG, Green, HD, Genes & Health Research Team, Hattersley, AT, L’Esperance, V, McKinley, TJ, Finer, S, Barroso, I. Undiagnosed G6PD deficiency in Black and Asian individuals is prevalent and contributes to health inequalities in type 2 diabetes diagnosis and complications. Diabetes Care, (in-press), 2025.

## EHR code lists for clinical phenotype defintions
EHR code lists for all clinical phenotypes used in the analysis are given in ICD-9, ICD-10, Read V2, Read CTV3, SNOMED and OPCS formats, as well as self-reported definitions from the UK Biobank 20002 and 20003 fields. For the UK Biobank, EHR data was curated from assessment centre interview (self-reported), primary care (Read V2 and Read CTV3 coded) and secondary care (ICD-9, ICD-10 and OPCS coded) sources[^2]. Details of the UK Biobank self-reported fields, assessment centre visits and values used to define each phenotype are given in **UKB_phenotype_definitions.xlsx**. For Genes & Health, EHR data was curated from primary care (SNOMED coded) and secondary care (ICD-10 coded) sources, in accordance with the type of data available.

[^2]: Green, HG. UKBB Health Care Records, <https://github.com/hdg204/UKBB> (2024).

EHR code lists used here that were previously published elsewhere are listed below with references:
- All files within **EHR_codes/Diabetes**, **20002_cvd.txt**, **20002_hypertension.txt**, **20003_blood_pressure_meds.txt**, **exeter_icd10_*.txt** - from Young et al. (2023; 2024)[^3][^4] **SNOMED_diabetes.csv** contains the SNOMED codes for the ‘DM_COD’ section of the v44 Quality and Outcomes Framework (QOF)[^5].
- **MULTIPLY_*.csv** - from Eto et al. (2023)[^6]
- **read_codes_pregnancy_Minassian_et_al.txt** - from Minassian et al. (2019)[^7]
- **LSHTM_CVD_ICD10_codes.csv** - from Forbes et al. (2018)[^8]
- **Ritchie_20003_*.csv** - from Ritchie et al. (2024)[^9]

[^3]: Young, KG, McGovern, AP, Barroso, I, et al. The impact of population-level HbA(1c) screening on reducing diabetes diagnostic delay in middle-aged adults: a UK Biobank analysis. Diabetologia (2023); 66:300-9.
[^4]: Young, KG. UK Biobank codelists, <https://github.com/drkgyoung/UK_Biobank_codelists> (2024).
[^5]: Quality and Outcomes Framework (QOF) v44, <https://digital.nhs.uk/data-and-information/data-collections-and-data-sets/data-collections/quality-and-outcomes-framework-qof/quality-and-outcome-framework-qof-business-rules/quality-and-outcomes-framework-qof-business-rules-v44-2019-2020-october-2020-release> (2023).
[^6]: Eto, F, Samuel, M, Finer, S. MULTIPLY initiative, <https://github.com/Fabiola-Eto/MULTIPLY-Initiative> (2023).
[^7]: Minassian, C, Williams, R, Meeraus, WH, Smeeth, L, Campbell, OMR, Thomas, SL. Methods to generate and validate a Pregnancy Register in the UK Clinical Practice Research Datalink primary care database. Pharmacoepidemiol Drug Saf (2019); 28:923-33.
[^8]: Forbes H, Langan S. Clinical codelist - CVD ICD-10 codes. London School of Hygiene & Tropical Medicine, London, United Kingdom, 2018.
[^9]: Ritchie SC, Taylor HJ, Liang Y, Manikpurage HD, Pennells L, Foguet C, Abraham G, Gibson JT, Jiang X, Liu Y, Xu Y, Kim LG, Mahajan A, McCarthy MI, Kaptoge S, Lambert SA, Wood A, Sim X, Collins FS, Denny JC, Danesh J, Butterworth AS, Di Angelantonio E, Inouye M. Integrated clinical risk prediction of type 2 diabetes with a multifactorial polygenic risk score. medRxiv (2024).

All other EHR code lists were curated by Susan Martin with Michael Barrington, Miriam Samuel and Inês Barroso.

## R and shell scripts

The shell script **extract_UKB_WES_genotypes.sh** is to be executed using Bash on JupyterLab on DNAnexus. It contains PLINK commands to extract the whole exome sequencing (WES) data for the two genetic variants of interest, and recode the genotypes according to the additive model. Equivalent commands were used to extract the WES data for Genes & Health (on TRE) and are available on request.

The below R scripts are split into phenotype data preparation and analysis. Names of the scripts include the cohort that these scripts are to be applied to - either UK Biobank (on DNAnexus) or Genes & Health (on TRE). Phenotype data preparation scripts for Genes & Health are comparable and are available on request.

R scripts are to be run in the following order (after running the above shell script):
- **phenotype_prep_*.R**: Script to prepare phenotype data - uses EHR codes and defines all clinical phenotypes used in analyses
- **G6PD_analysis_*.R**: Script to combine genetic and phenotype data, apply necessary quality control to data and conduct all statistical analyses.

## UK Biobank version date and fields extracted
The UK Biobank data used for this analysis (stored on DNAnexus) was last updated on March 1st 2024 (total N=502,185). The UK Biobank fields extracted from the cohort dataset into *_participant.csv files (using the 'Table Exporter' tool), and which were then used in the R scripts, are given in **UKB_table_exporter_fields.xlsx**.

## Model diagnostic plots
Model diagnostic plots for the linear regression models used in the analysis of age of type 2 diagnosis and presence of complications, and described in the **Supplemental Material**. Names of plot files refer to the cohort, ethnicity grouping and sex that the model was applied to.

## Matching diagnostic plots
Diagnostic and love plots for the matching of participants by QDiabetes-2018 score A and risk factors used in the QDiabetes-2018 score analysis and matching, and described in the **Supplemental Material**. Names of plot files refer to the ethnicity grouping and sex that the matching was applied to.

Model and matching diagnostic plots not currently accessable will be made available within the coming days.
