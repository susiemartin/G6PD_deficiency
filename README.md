# Undiagnosed G6PD deficiency is common in Black individuals in the UK, and contributes to health inequalities in type 2 diabetes diagnosis

Electronic health record (EHR) code lists and R scripts required for analysis of this project, as well as model diagnostic plots referenced in the Supplementary Appendix.

## EHR code lists for clinical phenotype defintions
EHR code lists for all conditions used in the analysis are given in read v2, read v3, ICD-9, ICD-10 and SNOMED format, as well as self-reported definitions from the UK Biobank 20002 and 20003 fields. Details of the UK Biobank fields, assessment centre visits and values used to define each phenotype is given in UKB_phenotype_definitions.csv.

## Scripts

### UK Biobank

PLINK commands to extract whole exome sequencing (WES) data:

Below R scripts to be used with UK Biobank data on DNAnexus. Scripts are split in data preparation (prep_*.R) and analysis (analysis_*.R).

Data preparation scripts are to be run in the following order:
    prep_genetic_data_step1.R: Script to prepare genetic data (after applying above PLINK commands)
    prep_phenotype_data_step2.R: Script to prepare phenotype data - uses EHR codes and defines all clinical phenotypes used in analyses
    prep_data_cleaning_step3.R: combines genetic and phenotype data and applies necessary quality control to this data.

Scripts to conduct statistical analyses:
    analysis_prevalence.R: Calculates prevalence of G6PD variants and G6PD deficiency diagnoses
    analysis_HbA1c.R: Comparison of HbA1c and random glucose across genotypes in individuals without diabetes or pregnancy
    analysis_age_diagnosis.R: Comparison of age of diagnosis of type 2 diabetes across genotpyes in individuals diagnosed since 2011.
    
## UK Biobank version date and fields extracted
The UK Biobank data used for this analysis (stored on DNAnexus) was last updated on March 1st 2024 (total N=502,185). The UK Biobank fields extracted from the cohort dataset (using the 'Table Exporter' tool) into *.csv files that are used in the R scripts are given in UKB_table_exporter_fields.xlsx.
