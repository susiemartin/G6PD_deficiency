# Undiagnosed G6PD deficiency is common in Black individuals in the UK, and contributes to health inequalities in type 2 diabetes diagnosis

Electronic health record (EHR) code lists and R scripts required for analysis of this project, as well as model diagnostic plots referenced in the Supplementary Appendix.

# EHR code lists for clinical phenotype defintions
The UK Biobank data used for this analysis (stored on DNAnexus) was last updated on March 1st 2024 (total N=502,185). Tables showing which UK Biobank fields were extracted using the DNAnexus 'Table Exporter' tool into each of the .csv files read into the R scripts are given in UKB_table_exporter_fields.xlsx.

EHR code lists for all conditions used in the analysis are given in read v2, read v3, ICD-9, ICD-10 and SNOMED format, as well as self-reported definitions from the UK Biobank fields.

# R scripts
# UK Biobank
R scripts to be used with UK Biobank data on DNAnexus. Scripts are split in data prep (prep_*.R) and sub-analysis (analysis_*.R).

Data prep scripts to be run in the following order:

    step1_genetic_data_prep.R - script to prep genetic data (after applying PLINK commands to extract data for R)
    step2_phenotype_data_prep.R - script to prep phenotype data, uses EHR codes and defines all conditions used in analyses
    step3_data_cleaning.R - combines genetic and phenotype data and applies QC to this data.

Scripts for research question analyses:

    Research_Question_1.R - applies analysis for research question 1 (prevalence of G6PD variants and G6PD deficiency / haemolytic anaemia diagnoses)
    Research_Question_2.R - applies main analysis for research question 2 (comparison of HbA1c and glucose in individuals without diabetes)
    Research_Question_3.R - applies analysis for research question 3a (comparison of age of diagnosis of type 2 diabetes).
    
