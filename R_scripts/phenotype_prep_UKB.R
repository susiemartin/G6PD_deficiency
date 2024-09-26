# Script 2 - extract phenotype data from the UK Biobank and generate variables for case definitions
# 3 January 2024

# Load required packages
library(dplyr)
library(stringr)
library(lubridate)

# Read in Harry's script to pull-out primary care data
library(devtools) 
source_url("https://raw.githubusercontent.com/hdg204/UKBB/main/UKBB_Health_Records_New_Project.R")

system("dx cd /G6PD_diabetes_comps")

##### Prepare the phenotype data required for analysis #####

# Read in phenotype data
system("dx download /Phenotype_data/baseline_participant.csv")
alldata <- read.csv("baseline_participant.csv", header = TRUE, stringsAsFactors = FALSE)

# Determine self-reported ethnicity from multiple instances
alldata$ethnicity <- alldata$Ethnic.background...Instance.0
alldata[is.na(alldata$ethnicity),]$ethnicity <- alldata[is.na(alldata$ethnicity),]$Ethnic.background...Instance.1
alldata[is.na(alldata$ethnicity),]$ethnicity <- alldata[is.na(alldata$ethnicity),]$Ethnic.background...Instance.2

# Generate census-based self-reported ethnic group variables to group together smaller ethnicities
alldata$ethnic_group <- NA
alldata[alldata$ethnicity %in% c("White", "British", "Irish", "Any other white background"),]$ethnic_group <- "White"
alldata[alldata$ethnicity %in% c("African", "Black or Black British", "Caribbean", "Any other Black background"),]$ethnic_group <- "Black"
alldata[alldata$ethnicity %in% c("Bangladeshi", "Indian", "Pakistani", "Chinese", "Any other Asian background"),]$ethnic_group <- "Asian"
alldata[alldata$ethnicity %in% c("Mixed", "White and Asian", "White and Black African", "White and Black Caribbean", "Any other mixed background"),]$ethnic_group <- "Mixed"
alldata[!is.na(alldata$ethnicity) & alldata$ethnicity == "Other ethnic group",]$ethnic_group <- "Other"

# Generate sub-ethnic group variable for sensitivity analysis
alldata$sub_ethnic_group <- NA
alldata[alldata$ethnicity %in% c("Bangladeshi", "Indian", "Pakistani"),]$sub_ethnic_group <- "South Asian"
alldata[!is.na(alldata$ethnicity) & alldata$ethnicity == "Chinese",]$sub_ethnic_group <- "East Asian"
alldata[!is.na(alldata$ethnicity) & alldata$ethnicity == "Any other Asian background",]$sub_ethnic_group <- "Other Asian"
alldata[!is.na(alldata$ethnicity) & alldata$ethnicity == "African",]$sub_ethnic_group <- "African"
alldata[!is.na(alldata$ethnicity) & alldata$ethnicity == "Caribbean",]$sub_ethnic_group <- "Caribbean"
alldata[alldata$ethnicity %in% c("Black or Black British", "Any other Black background"),]$sub_ethnic_group <- "Other black"

# Set sex as a factor variable with male as reference level
alldata$Sex <- relevel(as.factor(alldata$Sex), ref = "Male")

### Rename and generate new variables ###

# Rename main variables
colnames(alldata)[colnames(alldata) == "Glycated.haemoglobin..HbA1c....Instance.0"] <- "HbA1c"
colnames(alldata)[colnames(alldata) == "Glucose...Instance.0"] <- "Glucose"
colnames(alldata)[colnames(alldata) == "Fasting.time...Instance.0"] <- "Fasting"
colnames(alldata)[colnames(alldata) == "Body.mass.index..BMI....Instance.0"] <- "BMI"

# Generate single UK Biobank assessment centre variable
alldata$centre <- alldata$UK.Biobank.assessment.centre...Instance.0
alldata[is.na(alldata$centre),]$centre <- alldata[is.na(alldata$centre),]$UK.Biobank.assessment.centre...Instance.1
alldata[is.na(alldata$centre),]$centre <- alldata[is.na(alldata$centre),]$UK.Biobank.assessment.centre...Instance.2
alldata[is.na(alldata$centre),]$centre <- alldata[is.na(alldata$centre),]$UK.Biobank.assessment.centre...Instance.3

# Generate date of birth variable (assuming mid-point of month)
alldata$dob <- as.Date(paste0(as.character(alldata$Year.of.birth), "-", as.character(alldata$Month.of.birth), "-15"), "%Y-%b-%d")

# Generate variable for date of HbA1c reading
alldata$date_centre_0 <- as.Date(alldata$Date.of.attending.assessment.centre...Instance.0)

# Generate single age variable - for age at HbA1c reading
alldata$age <- alldata$Age.when.attended.assessment.centre...Instance.0
alldata[is.na(alldata$age),]$age <- (alldata[is.na(alldata$age),]$date_centre_0 - alldata[is.na(alldata$age),]$dob)/365.25

# Write function to identify any negative calculated ages or times and set to missing
neg_vars <- function(num_variable) {
  for (i in unique(num_variable)) {
    alldata$condition_age <- as.numeric(unlist(alldata[colnames(alldata) == i]))
    if (nrow(alldata[!is.na(alldata$condition_age) & alldata$condition_age < 0,]) > 0) {
      alldata[!is.na(alldata$condition_age) & alldata$condition_age < 0,]$condition_age <- NA
    }
    alldata <- alldata[names(alldata) != i]
    colnames(alldata)[colnames(alldata) == "condition_age"] <- i
  }
  alldata
}

alldata <- neg_vars(num_variable = "age")

### G6PD deficiency diagnoses ###

# Read in blood disorder codes for G6PD deficiency
system("dx download /EHR_codes/G6PD_deficiency/*")

# Identify all IDs with diagnosis of G6PD deficiency in GP/HES data
g6pdcodes <- NULL
for (code in c("read_2", "read_3")) {
  g6pddata <- read.csv(paste0(code, "_G6PD_deficiency.txt"), header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  g6pdcodes <- unique(c(g6pdcodes, g6pddata[,1]))
}
g6pdids <- read_GP(g6pdcodes)

icd10codes <- read.csv("ICD10_G6PD_deficiency.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
icd10codes <- unique(icd10codes$ICD10)
icd9codes <- read.csv("ICD9_G6PD_deficiency.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
icd9codes <- unique(icd9codes$Code)

icd10data <- read_ICD10(icd10codes)
icd9data <- read_ICD9(icd9codes)
g6pdids <- unique(c(g6pdids, icd10data$eid, icd9data$eid))

# Generate marker variable of G6PD deficiency diagnosis
alldata$g6pd <- 0
alldata[alldata$Participant.ID %in% g6pdids,]$g6pd <- 1

### Diabetes (any type) diagnoses (prevalent at baseline visit) - loose definition (for exclusion) ###

# Read in my files of pre-prepared diabetes codes from Katie
system("dx download /EHR_codes/Diabetes/*")

# Identify all IDs with any diabetes diagnosis or prescription of diabetes medication in GP/HES data (and dates of first occurrence)
# Note: No dates linked with ICD-9 codes, so these cannot be used
gpcodes <- NULL
for (code in c("read_2", "read_3")) {
  qofcodes <- read.table(paste0(code, "_diabetes_qof.txt"), header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  readcodes <- read.table(paste0(code, "_diabetes_diagnosis.txt"), header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  gpcodes <- unique(c(gpcodes, qofcodes[,1], readcodes[,1]))
}
icd10codes <- read.table("hes_icd10_diabetes.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
drugcodes <- read.table("read_2_diabetes_drugs.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
icddrugcodes <- read.table("ICD10_diabetes_drugs.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
diabdata <- first_occurence(ICD10 = c(icd10codes[,1], icddrugcodes[,1]), GP = c(gpcodes, drugcodes[,1]), OPCS = "", cancer = "")
colnames(diabdata)[colnames(diabdata) == "date"] <- "diabetes_date_EHR"
alldata <- merge(alldata, unique(diabdata[names(diabdata) != "source"]), by.x = "Participant.ID", by.y = "eid", all.x = TRUE)

# Write function to set UKB placeholder dates to missing or date of birth as appropriate
ukb_dates <- function(date_variable) {
  for (i in unique(date_variable)) {
    alldata$condition_date <- as.Date(unlist(alldata[colnames(alldata) == i]))
    # Find any dates of 1901-01-01 and set to missing - as before date of birth
    if (nrow(alldata[!is.na(alldata$condition_date) & alldata$condition_date == "1901-01-01",]) > 0) {
      alldata[!is.na(alldata$condition_date) & alldata$condition_date == "1901-01-01",]$condition_date <- NA
    }
    # Find any dates of 1902-02-02 or 1903-03-03 and set to date of birth - as within one month of date of birth
    if (nrow(alldata[!is.na(alldata$condition_date) & alldata$condition_date %in% c("1902-02-02", "1903-03-03"),]) > 0) {
      alldata[!is.na(alldata$condition_date) & alldata$condition_date %in% c("1902-02-02", "1903-03-03"),]$condition_date <- alldata[!is.na(alldata$condition_date) & alldata$condition_date %in% c("1902-02-02", "1903-03-03"),]$dob
    }
    alldata <- alldata[names(alldata) != i]
    colnames(alldata)[colnames(alldata) == "condition_date"] <- i
  }
  alldata
}

alldata <- ukb_dates(date_variable = "diabetes_date_EHR")

# Identify those who self-reported having diabetes at baseline assessment centre visit
alldata$diabetes_self <- 0
alldata[(!is.na(alldata$Diabetes.diagnosed.by.doctor...Instance.0) & alldata$Diabetes.diagnosed.by.doctor...Instance.0 == "Yes"),]$diabetes_self <- 1

# Read in lists of self-reported diabetes/diabetes-related condition or diabetes medication codes
selfdiabcodes <- read.table("20002_diabetes.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
selfdrugcodes <- read.table("20003_diabetes_meds.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")

# Read in self-reported age of non-cancer illness diagnosis
system("dx download /Phenotype_data/non_cancer_illness_participant.csv")
extradata <- read.csv("non_cancer_illness_participant.csv", header = TRUE, stringsAsFactors = FALSE, fill = TRUE, na.strings = "")
alldata <- merge(alldata, extradata, by = "Participant.ID")

system("dx download /Phenotype_data/non_cancer_illness_cont_participant.csv")
extradata <- read.csv("non_cancer_illness_cont_participant.csv", header = TRUE, stringsAsFactors = FALSE, fill = TRUE, na.strings = "")
alldata <- merge(alldata, extradata, by = "Participant.ID")

system("dx download /Phenotype_data/treatment_medication_participant.csv")
extradata <- read.csv("treatment_medication_participant.csv", header = TRUE, stringsAsFactors = FALSE, fill = TRUE, na.strings = "")
alldata <- merge(alldata, extradata, by = "Participant.ID")

# Set to missing any negative ages of non-cancer illness diagnosis (many columns here so do in one batch)
alldata[grepl("Interpolated.Age.of.participant.when.non.cancer.illness.first.diagnosed", colnames(alldata))] <- lapply(alldata[grepl("Interpolated.Age.of.participant.when.non.cancer.illness.first.diagnosed", colnames(alldata))], function(x) replace(x, x < 0, NA))

# Identify those who self-reported having diabetes/diabetes-related conditions or taking diabetes medication
# & extract self-reported age that diabetes/diabetes-related condition was diagnosed from corresponding array
eyedis0 <- strsplit(alldata$Eye.problems.disorders...Instance.0, "[|]")
index <- which(sapply(eyedis0, FUN = function(x) "Diabetes related eye disease" %in% x))
alldata$diabetes_self[index] <- 1

medall0 <- strsplit(alldata$Medication.for.cholesterol..blood.pressure..diabetes..or.take.exogenous.hormones...Instance.0, "[|]")
index <- which(sapply(medall0, FUN = function(x) "Insulin" %in% x))
alldata$diabetes_self[index] <- 1

medsub0 <- strsplit(alldata$Medication.for.cholesterol..blood.pressure.or.diabetes...Instance.0, "[|]")
index <- which(sapply(medsub0, FUN = function(x) "Insulin" %in% x))
alldata$diabetes_self[index] <- 1

trtmed0 <- alldata[grepl("Treatment.medication.code...Instance.0", colnames(alldata))]
index <- as.numeric(rownames(trtmed0[apply(trtmed0, 1, function(x) any(selfdrugcodes$Description %in% x)),]))
alldata$diabetes_self[index] <- 1

noncan0 <- alldata[grepl("Non.cancer.illness.code..self.reported...Instance.0", colnames(alldata))]
index <- as.numeric(rownames(noncan0[apply(noncan0, 1, function(x) any(selfdiabcodes$Description %in% x)),]))
alldata$diabetes_self[index] <- 1

# Identify the array/order that any diabetes/diabetes-related condition were reported - to extract corresponding age
noncan <- alldata[grepl("Non.cancer.illness.code..self.reported", colnames(alldata))]
noncan <- noncan[, order(colnames(noncan))]
age <- alldata[grepl("Interpolated.Age.of.participant.when.non.cancer.illness.first.diagnosed", colnames(alldata))]
age <- age[, order(colnames(age))]

# Write function to pull-out corresponding ages for self-reported non-cancer illnesses in array
noncancer_age <- function(disease_codes, condition_name) {
  alldata$condition_age_self <- NA
  allindex <- NULL
  for (i in unique(disease_codes)) {
    index <- which(noncan == i, arr.ind = TRUE)
    age_sub <- age[index]
    index <- data.frame(index)
    index$self_age <- age_sub
    allindex <- rbind(allindex, index)
  }
  allindex <- aggregate(self_age ~ row, data = allindex, min)
  alldata$condition_age_self[allindex$row] <- allindex$self_age
  colnames(alldata)[colnames(alldata) == "condition_age_self"] <- paste0(condition_name, "_age_self")
  alldata
}

alldata <- noncancer_age(disease_codes = selfdiabcodes$Description, condition_name = "diabetes")

# Identify earliest age that diabetes was self-reported
my_min <- function(x) ifelse(!all(is.na(x)), min(x, na.rm = TRUE), NA)
alldata[grepl("Age.diabetes.diagnosed", colnames(alldata))] <- as.numeric(unlist(alldata[grepl("Age.diabetes.diagnosed", colnames(alldata))]))
alldata$diabetes_age_self <- apply(alldata[grepl("Age.diabetes.diagnosed|diabetes_age_self", colnames(alldata))], 1, FUN = my_min)

# Estimate date of self-diagnosis from age of self-diagnosis
alldata$diabetes_date_self <- as.Date(alldata$dob %m+% years(floor(alldata$diabetes_age_self)))

# Generate earliest date of diabetes diagnosis from recorded EHR and self-reported
alldata$diabetes_date <- apply(alldata[c("diabetes_date_EHR", "diabetes_date_self")], 1, FUN = my_min)

# Assign IDs binary diabetes variable marker if diabetes diagnosis prior to baseline assessment centre
alldata$diabetes_loose <- 0
alldata[alldata$diabetes_self == 1 |
          (!is.na(alldata$diabetes_date) & !is.na(alldata$date_centre_0) & (alldata$diabetes_date <= alldata$date_centre_0)),]$diabetes_loose <- 1

### Pregnancy term ###

# Read in outcome of delivery/miscarriage/termination codes - to identify if individuals were pregnant at time of sample
# Note: No dates linked with ICD-9 codes, so these cannot be used
system("dx download /EHR_codes/Pregnancy_outcome/*")
icd10codes <- read.csv("ICD10_pregnancy.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
icd10codes12 <- unique(icd10codes[icd10codes$Term_group == 0,]$Code)
icd10codes24 <- unique(icd10codes[icd10codes$Term_group == 1,]$Code)
icd10codes40 <- unique(icd10codes[icd10codes$Term_group == 2,]$Code)

icd10data12 <- read_ICD10(icd10codes12)
icd10data24 <- read_ICD10(icd10codes24)
icd10data40 <- read_ICD10(icd10codes40)

# Read in outcome of pregnancy read codes from Minassian et al. (2019) and group into term length categories
readcodes <- read.csv("Minassian_read_codes_pregnancy.csv", header = TRUE, stringsAsFactors = FALSE)
readcodes$read.code <- substr(readcodes$read.code, 1, nchar(readcodes$read.code)-2)
readcodes$duration <- NA
readcodes[readcodes$ectopic == 1 | readcodes$molar == 1 | readcodes$blighted_ovum == 1,]$duration <- 0
readcodes[is.na(readcodes$duration) & (readcodes$miscarriage == 1 | readcodes$top == 1 | readcodes$top_probable == 1),]$duration <- 1
readcodes[is.na(readcodes$duration) & (readcodes$delivery == 1 | readcodes$multiple == 1 | readcodes$stillbirth == 1),]$duration <- 2

readdata12 <- read_GP(unique(readcodes[!is.na(readcodes$duration) & readcodes$duration == 0,]$read.code))
readdata24 <- read_GP(unique(readcodes[!is.na(readcodes$duration) & readcodes$duration == 1,]$read.code))
readdata40 <- read_GP(unique(readcodes[!is.na(readcodes$duration) & readcodes$duration == 2,]$read.code))

# Bind HES and GP pregnancy records with end dates together
readdata12$epiend <- readdata12$event_dt
readdata24$epiend <- readdata24$event_dt
readdata40$epiend <- readdata40$event_dt

pregdata12 <- rbind(icd10data12[c("eid", "epiend")], readdata12[c("eid", "epiend")])
pregdata24 <- rbind(icd10data24[c("eid", "epiend")], readdata24[c("eid", "epiend")])
pregdata40 <- rbind(icd10data40[c("eid", "epiend")], readdata40[c("eid", "epiend")])

# Estimate start date of pregnancy based on type of pregnancy outcome
pregdata12$preg_start_date <- as.Date(pregdata12$epiend) %m-% weeks(12)
pregdata24$preg_start_date <- as.Date(pregdata24$epiend) %m-% weeks(24)
pregdata40$preg_start_date <- as.Date(pregdata40$epiend) %m-% weeks(40)
pregdata <- rbind(pregdata12, pregdata24, pregdata40)

# Estimate date of 6-weeks postpartum - when HbA1c becomes reliable again
pregdata$postpartum_date <- as.Date(pregdata$epiend) %m+% weeks(6)
alldata <- merge(alldata, pregdata[c("eid", "preg_start_date", "postpartum_date")], by.x = "Participant.ID", by.y = "eid", all.x = TRUE)
alldata <- ukb_dates(date_variable = c("preg_start_date", "postpartum_date"))

# Generate binary variable for pregnant during baseline assessment centre visit
alldata$pregnant <- 0
alldata[(!is.na(alldata$Pregnant...Instance.0) & alldata$Pregnant...Instance.0 == "Yes") |
          (!is.na(alldata$date_centre_0) & !is.na(alldata$preg_start_date) & !is.na(alldata$postpartum_date) & (alldata$date_centre_0 >= alldata$preg_start_date) & (alldata$date_centre_0 < alldata$postpartum_date)),]$pregnant <- 1

# Remove columns relating to pregnancy dates so can have one row per ID
alldata <- alldata[, -which(names(alldata) %in% c("preg_start_date", "postpartum_date"))]
alldata <- alldata[!duplicated(alldata$Participant.ID), ]

### Type 2 diabetes diagnosis - strict definition ###

# Read in files of ALL diabetes QOF codes from Katie (everyone meeting their strict definition of diabetes (all))
katie_read2 <- read.csv("read_2_diabetes_qof.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
katie_read3 <- read.csv("read_3_diabetes_qof.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
katie_read <- unique(c(katie_read2$read_2, katie_read3$read_3))

## Use Harry's script to identify any diabetes diagnosis (using Katie's QOF codes)
diabgpdata <- first_occurence(ICD10 = "", GP = katie_read, OPCS = "", cancer = "")
colnames(diabgpdata)[colnames(diabgpdata) == "date"] <- "t2d_date_GP"

## Use Harry's script to identify first HES occurrence of T2D diagnosis
# Note: No dates linked with ICD-9 codes, so these cannot be used - but these are "249|250" for reference
system("dx download /EHR_codes/Type_2_diabetes/*")
icd10codes <- read.csv("ICD10_T2D.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
diabhesdata <- first_occurence(ICD10 = icd10codes[,1], GP = "", OPCS = "", cancer = "")
colnames(diabhesdata)[colnames(diabhesdata) == "date"] <- "t2d_date_HES"
diabdata <- merge(diabgpdata[names(diabgpdata) != "source"], diabhesdata[names(diabhesdata) != "source"], by = "eid", all = TRUE)

# Read in self-reported type 2 diabetes
alldata$t2d_self <- 0
index <- as.numeric(rownames(noncan[apply(noncan, 1, function(x) "type 2 diabetes" %in% x),]))
alldata$t2d_self[index] <- 1

# Identify the array/order that type 2 diabetes was reported - to extract corresponding age
alldata <- noncancer_age(disease_codes = "type 2 diabetes", condition_name = "t2d")

# Read in pre-generated table of diabetes-related baseline data
system("dx download /Phenotype_data/diabetes_info_participant.csv")
basedata <- read.csv("diabetes_info_participant.csv", header = TRUE, stringsAsFactors = FALSE)
alldata <- merge(alldata, basedata, by = "Participant.ID", all.x = TRUE)

alldata <- merge(alldata, diabdata, by.x = "Participant.ID", by.y = "eid", all.x = TRUE)
alldata <- ukb_dates(date_variable = c("t2d_date_GP", "t2d_date_HES"))
alldata$diabetes_strict <- 0
alldata[!is.na(alldata$t2d_date_GP) | !is.na(alldata$t2d_date_HES) | alldata$t2d_self == 1,]$diabetes_strict <- 1

# Read in files of non-T1D/T2D exclusion codes from Katie (e.g. genetic, gestational, etc.)
katie_excl2 <- read.csv("read_2_diabetes_other_types.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
katie_excl3 <- read.csv("read_3_diabetes_other_types.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
katie_excl <- unique(c(katie_excl2$Code, katie_excl3$Code))

# Read in ICD-10 non-T1D/T2D exclusion codes
exclcodes <- read.csv("ICD10_diabetes_other_types.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")

## Use Harry's script to identify non-T1D/T2D diabetes diagnoses to exclude (using Katie's exclusion codes)
excludedata <- first_occurence(ICD10 = exclcodes[,1], GP = katie_excl, OPCS = "", cancer = "")

# Generate a column to identify the T1D/T2D diagnoses
alldata$t1d_or_t2d <- 0
alldata[alldata$diabetes_strict == 1 & !(alldata$Participant.ID %in% excludedata$eid),]$t1d_or_t2d <- 1

## Identify insulin use (to later remove those having insulin within one year of diagnosis)

# Extract insulin prescription codes
insulin_codes <- unique(drugcodes[drugcodes$drug_substance == "insulin" | drugcodes$drug_class1 == "insulin" | drugcodes$drug_class2 == "insulin",]$Code)

## Use Harry's script to identify insulin GP prescriptions to exclude later (using Katie's drug codes)
insulin_data <- first_occurence(ICD10 = "", GP = insulin_codes, OPCS = "", cancer = "")
colnames(insulin_data)[colnames(insulin_data) == "date"] <- "insulin_date_GP"
alldata <- merge(alldata, insulin_data, by.x = "Participant.ID", by.y = "eid", all.x = TRUE)
alldata <- ukb_dates(date_variable = c("insulin_date_GP"))

## Use Harry's script to identify insulin hospital records to exclude later (using Katie's drug codes)
icd10codes <- read.csv("ICD10_insulin.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
insulin_data <- first_occurence(ICD10 = icd10codes[,1], GP = "", OPCS = "", cancer = "")
colnames(insulin_data)[colnames(insulin_data) == "date"] <- "insulin_date_HES"
alldata <- merge(alldata, insulin_data, by.x = "Participant.ID", by.y = "eid", all.x = TRUE)
alldata <- ukb_dates(date_variable = c("insulin_date_HES"))

# Read in self-reported treatment/medication data for non-baseline assessment centre visits
system("dx download /Phenotype_data/treatment_med_cont_participant.csv")
extradata <- read.csv("treatment_med_cont_participant.csv", header = TRUE, stringsAsFactors = FALSE, fill = TRUE, na.strings = "")
alldata <- merge(alldata, extradata, by = "Participant.ID", all.x = TRUE)

# Identify those who self-reported as taking insulin at any assessment centre visit and align with date of visit - for very approximate estimate of first use
medcolnames <- c("Medication.for.cholesterol..blood.pressure..diabetes..or.take.exogenous.hormones", "Medication.for.cholesterol..blood.pressure.or.diabetes")

alldata$insulin_date_self_0 <- NA
index <- which(sapply(medall0, FUN = function(x) "Insulin" %in% x))
alldata$insulin_date_self_0[index] <- alldata$date_centre_0[index]
index <- which(sapply(medsub0, FUN = function(x) "Insulin" %in% x))
alldata$insulin_date_self_0[index] <- alldata$date_centre_0[index]
index <- as.numeric(rownames(trtmed0[apply(trtmed0, 1, function(x) "insulin product" %in% x),]))
alldata$insulin_date_self_0[index] <- alldata$date_centre_0[index]
alldata$insulin_date_self_0 <- as.Date(alldata$insulin_date_self_0)

alldata$insulin_date_self_1 <- NA
for (i in paste0(medcolnames, "...Instance.1")) {
  medall <- strsplit(alldata[,i], "[|]")
  index <- which(sapply(medall, FUN = function(x) "Insulin" %in% x))
  alldata$insulin_date_self_1[index] <- alldata$Date.of.attending.assessment.centre...Instance.1[index]
}
trtmed1 <- alldata[grepl("Treatment.medication.code...Instance.1", colnames(alldata))]
index <- as.numeric(rownames(trtmed1[apply(trtmed1, 1, function(x) "insulin product" %in% x),]))
alldata$insulin_date_self_1[index] <- alldata$Date.of.attending.assessment.centre...Instance.1[index]
alldata$insulin_date_self_1 <- as.Date(alldata$insulin_date_self_1)

alldata$insulin_date_self_2 <- NA
for (i in paste0(medcolnames, "...Instance.2")) {
  medall <- strsplit(alldata[,i], "[|]")
  index <- which(sapply(medall, FUN = function(x) "Insulin" %in% x))
  alldata$insulin_date_self_2[index] <- alldata$Date.of.attending.assessment.centre...Instance.2[index]
}
trtmed2 <- alldata[grepl("Treatment.medication.code...Instance.2", colnames(alldata))]
index <- as.numeric(rownames(trtmed2[apply(trtmed2, 1, function(x) "insulin product" %in% x),]))
alldata$insulin_date_self_2[index] <- alldata$Date.of.attending.assessment.centre...Instance.2[index]
alldata$insulin_date_self_2 <- as.Date(alldata$insulin_date_self_2)

alldata$insulin_date_self_3 <- NA
for (i in paste0(medcolnames, "...Instance.3")) {
  medall <- strsplit(alldata[,i], "[|]")
  index <- which(sapply(medall, FUN = function(x) "Insulin" %in% x))
  alldata$insulin_date_self_3[index] <- alldata$Date.of.attending.assessment.centre...Instance.3[index]
}
trtmed3 <- alldata[grepl("Treatment.medication.code...Instance.3", colnames(alldata))]
index <- as.numeric(rownames(trtmed3[apply(trtmed3, 1, function(x) "insulin product" %in% x),]))
alldata$insulin_date_self_3[index] <- alldata$Date.of.attending.assessment.centre...Instance.3[index]
alldata$insulin_date_self_3 <- as.Date(alldata$insulin_date_self_3)

# Calculate earliest estimated date of self-reported insulin use
alldata$insulin_date_self <- apply(alldata[paste0("insulin_date_self_", 0:3)], 1, FUN = my_min)

# Use earliest date of self-reported estimate insulin use and GP prescription/hospital records - and create variable for origin of date (for potential use in model adjustment later on)
alldata$insulin_date <- as.Date(apply(alldata[c("insulin_date_GP", "insulin_date_HES", "insulin_date_self")], 1, FUN = my_min))

alldata$t2d_age_self <- apply(alldata[grepl("Age.diabetes.diagnosed|t2d_age_self", colnames(alldata))], 1, FUN = my_min)

alldata$insulin_first_year <- 0
alldata[(!is.na(alldata$Started.insulin.within.one.year.diagnosis.of.diabetes...Instance.0) & alldata$Started.insulin.within.one.year.diagnosis.of.diabetes...Instance.0 == "Yes") |
          (!is.na(alldata$Started.insulin.within.one.year.diagnosis.of.diabetes...Instance.1) & alldata$Started.insulin.within.one.year.diagnosis.of.diabetes...Instance.1 == "Yes") |
          (!is.na(alldata$Started.insulin.within.one.year.diagnosis.of.diabetes...Instance.2) & alldata$Started.insulin.within.one.year.diagnosis.of.diabetes...Instance.2 == "Yes") |
          (!is.na(alldata$Started.insulin.within.one.year.diagnosis.of.diabetes...Instance.3) & alldata$Started.insulin.within.one.year.diagnosis.of.diabetes...Instance.3 == "Yes"),]$insulin_first_year <- 1

# Calculate first-occurrence age of T2D diagnosis using date of birth and first occurrence of T2D diagnosis
alldata$t2d_age_GP <- as.numeric(alldata$t2d_date_GP - alldata$dob)/365.25
alldata$t2d_age_HES <- as.numeric(alldata$t2d_date_HES - alldata$dob)/365.25
alldata <- neg_vars(num_variable = c("t2d_age_GP", "t2d_age_HES"))

# Set earliest age of T2D diagnosis to minimum of self-reported, GP and HES diagnosis ages
alldata$t2d_age <- alldata$t2d_age_self
alldata$t2d_age_origin <- NA
alldata[!is.na(alldata$t2d_age),]$t2d_age_origin <- "Selfreport"

alldata[(is.na(alldata$t2d_age) & !is.na(alldata$t2d_age_GP)) | (!is.na(alldata$t2d_age_GP) & !is.na(alldata$t2d_age) & alldata$t2d_age_GP < alldata$t2d_age),]$t2d_age_origin <- "GP"
alldata[!is.na(alldata$t2d_age_origin) & alldata$t2d_age_origin == "GP",]$t2d_age <- alldata[!is.na(alldata$t2d_age_origin) & alldata$t2d_age_origin == "GP",]$t2d_age_GP

# Set earliest age of diagnosis to minimum of self-reported and GP diagnosis ages (for age of diagnosis analysis)
alldata$t2d_age_GP_self_origin <- alldata$t2d_age_origin
alldata$t2d_age_GP_self_min <- alldata$t2d_age

alldata[(is.na(alldata$t2d_age) & !is.na(alldata$t2d_age_HES)) | (!is.na(alldata$t2d_age_HES) & !is.na(alldata$t2d_age) & alldata$t2d_age_HES < alldata$t2d_age),]$t2d_age_origin <- "HES"
alldata[!is.na(alldata$t2d_age_origin) & alldata$t2d_age_origin == "HES",]$t2d_age <- alldata[!is.na(alldata$t2d_age_origin) & alldata$t2d_age_origin == "HES",]$t2d_age_HES

# Generate binary variables to indicate whether diagnosis was provided in self-reported, GP or HES data
alldata$t2d_GP <- alldata$t2d_HES <- 0
alldata[!is.na(alldata$t2d_age_GP),]$t2d_GP <- 1
alldata[!is.na(alldata$t2d_age_HES),]$t2d_HES <- 1

# Estimate date of self-reported T2D diagnosis using date of birth and age of diagnosis
alldata$t2d_date_self <- as.Date(alldata$dob %m+% years(floor(alldata$t2d_age_self)))

# Create variable for earliest (recorded or estimated) date of T2D diagnosis
alldata$t2d_date <- alldata$t2d_date_self
alldata[!is.na(alldata$t2d_age_origin) & alldata$t2d_age_origin == "GP",]$t2d_date <- alldata[!is.na(alldata$t2d_age_origin) & alldata$t2d_age_origin == "GP",]$t2d_date_GP
alldata[!is.na(alldata$t2d_age_origin) & alldata$t2d_age_origin == "HES",]$t2d_date <- alldata[!is.na(alldata$t2d_age_origin) & alldata$t2d_age_origin == "HES",]$t2d_date_HES

alldata$t2d_date_GP_self_min <- alldata$t2d_date_self
alldata[!is.na(alldata$t2d_age_GP_self_origin) & alldata$t2d_age_GP_self_origin == "GP",]$t2d_date_GP_self_min <- alldata[!is.na(alldata$t2d_age_GP_self_origin) & alldata$t2d_age_GP_self_origin == "GP",]$t2d_date_GP

# Generate binary variable to identify whether diagnosed with diabetes BEFORE UK Biobank assessment centre
alldata$t2d_before_centre <- 0
alldata[!is.na(alldata$t2d_date) & !is.na(alldata$date_centre_0) & (alldata$t2d_date <= alldata$date_centre_0),]$t2d_before_centre <- 1

# Make binary variable to define if had insulin within one year of diagnosis
alldata$insulin_first_year2 <- 0
alldata$years_to_insulin <- as.numeric(alldata$insulin_date - alldata$t2d_date)/365.25
alldata <- neg_vars(num_variable = c("years_to_insulin"))

alldata[!is.na(alldata$years_to_insulin) & alldata$years_to_insulin < 1,]$insulin_first_year2 <- 1
alldata[alldata$insulin_first_year2 == 1,]$insulin_first_year <- 1

# Create a T2D indicator variable (with and without <35y age restriction)
alldata$t2d_anyage <- 0
alldata[alldata$diabetes_strict == 1 & alldata$t1d_or_t2d == 1 & !is.na(alldata$insulin_first_year) & alldata$insulin_first_year == 0,]$t2d_anyage <- 1

alldata$t2d_over35 <- 0
alldata[alldata$t2d_anyage == 1 & !is.na(alldata$t2d_age) & alldata$t2d_age >= 36,]$t2d_over35 <- 1

# Set age and date of T2D diagnosis missing if T2D binary variable is zero
alldata[alldata$t2d_anyage == 0,]$t2d_age <- NA
alldata[alldata$t2d_anyage == 0,]$t2d_date <- NA

## Save the full complete dataframe with all columns in case need it again later
write.table(alldata, file = "G6PD_phenotype_data_FULL.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
system("dx upload G6PD_phenotype_data_FULL.txt")

### Keep only phenotype codes will need in upcoming analyses ###
alldata <- alldata[c("Participant.ID", "Sex", "ethnicity", "ethnic_group", "sub_ethnic_group", "HbA1c", "Glucose", "BMI", "centre", "age",
                     "g6pd", "diabetes_loose", "pregnant",
                     "t2d_HES", "t2d_GP", "t2d_self", "t2d_date", "t2d_over35", "t2d_anyage", "t2d_age", "t2d_age_origin", "t2d_age_GP_self_min", "t2d_age_GP_self_origin", "t2d_date_GP_self_min")]

### Save final phenotype dataframe for merging with genotypes ###
write.table(alldata, file = "G6PD_phenotype_data.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
system("dx upload G6PD_phenotype_data.txt")

