# Undiagnosed G6PD deficiency in Black and Asian individuals is prevalent and contributes to health inequalities in type 2 diabetes diagnosis

This repository contains the electronic health record (EHR) code lists, R scripts and PLINK commands required for analysis of this project, as well as model diagnostic plots for the linear regression models detailed in the **Supplementary Appendix**.

## EHR code lists for clinical phenotype defintions
EHR code lists for all clinical phenotypes used in the analysis are given in read v2, read v3, ICD-9, ICD-10 and SNOMED format, as well as self-reported definitions from the UK Biobank 20002 and 20003 fields. For the UK Biobank, EHR data was curated from assessment centre interview (self-reported), primary care (read v2 and v3 coded) and secondary care (ICD-9 and ICD-10 coded) sources<sup>1</sup>. Details of the UK Biobank self-reported fields, assessment centre visits and values used to define each phenotype are given in **UKB_phenotype_definitions.xlsx**. For Genes & Health, EHR data was curated from primary care (SNOMED coded) and secondary care (ICD-10 coded) sources, in accordance with the type of data available.

EHR code lists used here that were previously published elsewhere are listed below with references:
- **20002_T2D.txt**, **20002_diabetes.txt**, **20003_diabetes_meds.txt**, **read_2_diabetes_diagnosis.txt**, **read_2_diabetes_drugs.txt**, **read_2_diabetes_other_types.txt**, **read_2_diabetes_qof.txt**, **read_3_diabetes_diagnosis.txt**, **read_3_diabetes_other_types.txt**, **read_3_diabetes_qof.txt**, **SNOMED_diabetes.csv** - from Young et al.<sup>2,3</sup> SNOMED_diabetes.csv contains the SNOMED codes for the ‘DM_COD’ section of the v44 Quality and Outcomes Framework (QOF)<sup>4</sup>.
- **MULTIPLY_Snomed_T2D.csv** - from Eto et al.<sup>5</sup>
- **read_codes_pregnancy_Minassian_et_al.txt** - from Minassian et al.<sup>6</sup>

All other EHR code lists were curated by Susan Martin with Michael Barrington, Miriam Samuel and Inês Barroso.

## R and shell scripts

The shell script **extract_UKB_WES_genotypes.sh** is to be executed using Bash on JupyterLab on DNAnexus. It contains PLINK commands to extract the whole exome sequencing (WES) data for the two genetic variants of interest, and recode the genotypes according to the additive model. Equivalent commands were used to extract the WES data for Genes & Health (on TRE).

The below R scripts are split into phenotype data preparation and analysis. Names of the scripts include the cohort that these scripts are to be applied to - either UK Biobank (on DNAnexus) or Genes & Health (on TRE). Phenotype data preparation scripts for Genes & Health are comparable and will be available soon.

R scripts are to be run in the following order (after running the above shell script):
- **phenotype_prep_*.R**: Script to prepare phenotype data - uses EHR codes and defines all clinical phenotypes used in analyses
- **G6PD_analysis_*.R**: Script to combine genetic and phenotype data, apply necessary quality control to data and conduct all statistical analyses.

## UK Biobank version date and fields extracted
The UK Biobank data used for this analysis (stored on DNAnexus) was last updated on March 1st 2024 (total N=502,185). The UK Biobank fields extracted from the cohort dataset into *.csv files (using the 'Table Exporter' tool), and which were then used in the R scripts, are given in **UKB_table_exporter_fields.xlsx**.

## Model diagnostic plots
Model diagnostic plots for the linear regression model used in the analysis of age of type 2 diagnosis and described in the **Supplementary Appendix**. Names of plot files refer to the cohort, ethnicity grouping and sex that the model was applied to.

## References
1. Green, H.G. UKBB Health Care Records, <https://github.com/hdg204/UKBB> (2024).
2. Young, K.G., McGovern, A.P., Barroso, I., et al. The impact of population-level HbA(1c) screening on reducing diabetes diagnostic delay in middle-aged adults: a UK Biobank analysis. Diabetologia (2023); 66:300-9.
3. Young, K.G. UK Biobank codelists, <https://github.com/drkgyoung/UK_Biobank_codelists> (2024).
4. Quality and Outcomes Framework (QOF) v44, <https://digital.nhs.uk/data-and-information/data-collections-and-data-sets/data-collections/quality-and-outcomes-framework-qof/quality-and-outcome-framework-qof-business-rules/quality-and-outcomes-framework-qof-business-rules-v44-2019-2020-october-2020-release> (2023).
5. Eto, F., Samuel, M., Finer, S. MULTIPLY initiative, <https://github.com/Fabiola-Eto/MULTIPLY-Initiative> (2023).
6. Minassian, C., Williams, R., Meeraus, W.H., Smeeth, L., Campbell, O.M.R., Thomas, S.L. Methods to generate and validate a Pregnancy Register in the UK Clinical Practice Research Datalink primary care database. Pharmacoepidemiol Drug Saf (2019); 28:923-33.

