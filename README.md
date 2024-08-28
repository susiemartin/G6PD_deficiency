# G6PD and age of diagnosis analysis - NEJM letter

Electronic health record (EHR) code lists and R scripts required for analysis of G6PD NEJM project - identification of G6PD carriers, G6PD deficiency / haemolytic anaemia diagnoses, and age of type 2 diabetes diagnosis.

The UK Biobank data stored on DNAnexus used for this analysis was in a project created on March 1st 2024 (total N=X). Tables showing which UK Biobank fields were extracted using the DNAnexus 'Table Exporter' tool into each of the .csv files read into the R scripts are given and denoted as 'UKB_table_exporter_NameOfCsvFile.xlsx'.

EHR code lists for all conditions to be used in analysis are given in read, ICD, OPCS and SNOMED format, as well as self-reported definitions from the UK Biobank 20002/20003 variable fields.

R scripts to be used with UK Biobank data on DNAnexus. Scripts are split in data prep (step*.R) and reseach question analysis (Research_Question*.R).

Data prep scripts to be run in the following order:

    step1_genetic_data_prep.R - script to prep genetic data (after applying PLINK commands to extract data for R)
    step2_phenotype_data_prep.R - script to prep phenotype data, uses EHR codes and defines all conditions used in analyses
    step3_data_cleaning.R - combines genetic and phenotype data and applies QC to this data.

Scripts for research question analyses:

    Research_Question_1.R - applies analysis for research question 1 (prevalence of G6PD variants and G6PD deficiency / haemolytic anaemia diagnoses)
    Research_Question_2.R - applies main analysis for research question 2 (comparison of HbA1c and glucose in individuals without diabetes)
    Research_Question_2_sensitivity.R - applies sensitivity analysis for research question 2 (comparison of HbA1c and glucose in individuals without conditions that affect HbA1c levels)
    Research_Question_2_T2D.R - applies additional analysis for research question 2 (comparison of HbA1c and glucose in individuals with type 2 diabetes)
    Research_Question_3.R - applies analysis for research question 3a (comparison of age of diagnosis of type 2 diabetes).
    
