# Script to run G6PD diagnostic age project analysis in the UK Biobank
# 20 August 2024

# Load and install required packages
install.packages("metafor")
install.packages("patchwork")
library(stringr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(metafor)
library(patchwork)
library(gridExtra)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

system("dx cd /G6PD_diabetes_comps")

## Read in genetic data

# Read in PLINK G6PD variant genotyping results file
system("dx download /PLINK_genotypes/G6PD_diabetes_comps/genotype_files/var_plink_chrX_genotype.txt")
gendata <- read.table(paste0("var_plink_chrX_genotype.txt"), header = TRUE, stringsAsFactors = FALSE)
gendata <- gendata[, -which(names(gendata) %in% c("FID", "PAT", "MAT", "SEX", "PHENOTYPE"))]
colnames(gendata)[colnames(gendata) == "IID"] <- "ID"

## Read in phenotype data to merge together

system("dx download G6PD_phenotype_data.txt")
phenodata <- read.csv("G6PD_phenotype_data.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")

alldata <- merge(gendata, phenodata, by.x = "ID", by.y = "Participant.ID")

colnames(alldata)[colnames(alldata) == "X23.154534419.G.A_A"] <- "23:154534419:G:A"
colnames(alldata)[colnames(alldata) == "X23.154536002.C.T_T"] <- "23:154536002:C:T"

# Remove any IDs without self-reported ethnic group
alldata <- alldata[!is.na(alldata$ethnic_group),]

# Exclude all men with genotypes of 1 (heterozygous) for either X-chr variant - as either Klinefelter syndrome or error
ksids <- unique(alldata[alldata$Sex == "Male" &
                          ((!is.na(alldata$`23:154536002:C:T`) & alldata$`23:154536002:C:T` == 1) | (!is.na(alldata$`23:154534419:G:A`) & alldata$`23:154534419:G:A` == 1)),]$ID)
length(ksids)
alldata <- alldata[!(alldata$ID %in% ksids),]

# Set to missing any uncorrected HbA1c measures > 184 mmol/mol (maximum detectable value)
alldata[!is.na(alldata$HbA1c) & alldata$HbA1c > 184,]$HbA1c <- NA

# Apply UK Biobank HbA1c correction (from Young et al. (2022))
alldata$uncorrected_HbA1c <- alldata$HbA1c
alldata$HbA1c <- (0.9696*alldata$HbA1c) + 3.3595

## Generate summary characteristics of UK Biobank cohort

length(unique(alldata$ID))

length(unique(alldata[alldata$Sex == "Male",]$ID))
length(unique(alldata[alldata$Sex == "Female",]$ID))

length(unique(alldata[alldata$ethnic_group == "White",]$ID))
length(unique(alldata[alldata$ethnic_group == "Black",]$ID))
length(unique(alldata[alldata$ethnic_group == "Asian",]$ID))
length(unique(alldata[alldata$ethnic_group == "Mixed",]$ID))
length(unique(alldata[alldata$ethnic_group == "Other",]$ID))

mean(alldata$age, na.rm = TRUE)
sd(alldata$age, na.rm = TRUE)

mean(alldata$BMI, na.rm = TRUE)
sd(alldata$BMI, na.rm = TRUE)

mean(alldata$HbA1c, na.rm = TRUE)
sd(alldata$HbA1c, na.rm = TRUE)

mean(alldata$Glucose, na.rm = TRUE)
sd(alldata$Glucose, na.rm = TRUE)

length(unique(alldata[alldata$g6pd == 1,]$ID))

## Determine the MAF for the IDs remaining in the cohort (e.g. those with pheno/ethnic data) - for self-reported ethnicity

ethnicities <- c("White", "Black", "Asian", "Mixed", "Other")
mafdata <- data.frame(matrix(NA, nrow = 3*2*length(ethnicities), ncol = 4))
colnames(mafdata) <- c("Variant", "Sex", "Ethnicity", "MAF")
mafdata$Variant <- rep(c("23:154536002:C:T", "23:154534419:G:A"), each = 3*length(ethnicities))
mafdata$Sex <- rep(c("Both", "Male", "Female"), each = length(ethnicities))
mafdata$Ethnicity <- ethnicities
for (i in 1:nrow(mafdata)) {
  subdata <- alldata[!is.na(alldata[mafdata$Variant[i]]),]
  if (mafdata$Sex[i] != "Both") {
    subdata <- subdata[subdata$Sex == mafdata$Sex[i],]
  }
  if (mafdata$Ethnicity[i] %in% c("White", "Black", "Asian", "Mixed", "Other")) {
    subdata <- subdata[subdata$ethnic_group == mafdata$Ethnicity[i],]
  } else {
    subdata <- subdata[!is.na(subdata$sub_ethnic_group) & subdata$sub_ethnic_group == mafdata$Ethnicity[i],]
  }
  mafdata$MAF[i] <- sum(subdata[mafdata$Variant[i]])/(2*nrow(subdata))
}
mafdata <- mafdata[order(mafdata$Variant, decreasing = TRUE),]
mafdata$Prevalence <- 100*mafdata$MAF

write.table(mafdata, file = "MAF_selfreport_pheno_cohort_only.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
system("dx upload MAF_selfreport_pheno_cohort_only.txt")

## Determine the prevalence of G6PD deficiency diagnoses in carriers

myresults <- data.frame(matrix(NA, nrow = 2*length(ethnicities)*5, ncol = 6))
colnames(myresults) <- c("Variant", "Ethnicity", "Sex", "Genotype", "N", "G6PD_deficiency")
myresults$Variant <- rep(c("23:154536002:C:T", "23:154534419:G:A"), each = length(ethnicities)*5)
myresults$Ethnicity <- rep(ethnicities, each = 5)
myresults$Sex <- c(rep("Male", 2), rep("Female", 3))
myresults$Genotype <- c(0,2,0:2)
for (i in 1:nrow(myresults)) {
  subdata <- alldata[alldata$Sex == myresults$Sex[i],]
  subdata <- subdata[subdata$ethnic_group == myresults$Ethnicity[i],]
  subdata <- subdata[!is.na(subdata[myresults$Variant[i]]) & subdata[myresults$Variant[i]] == myresults$Genotype[i],]
  myresults$N[i] <- length(unique(subdata$ID))
  myresults$G6PD_deficiency[i] <- length(unique(subdata[subdata$g6pd == 1,]$ID))
}

write.table(myresults, file = "genotype_g6pd.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
system("dx upload genotype_g6pd.txt")

### HbA1c and random glucose analyses ###

## Determine participants meeting exclusion criteria for HbA1c and glucose analyses

# Find number of starting individuals before applying any exclusion criteria (UK Biobank)
length(unique(alldata$ID))

# Find number with missing variables
length(unique(alldata[is.na(alldata$HbA1c),]$ID))
length(unique(alldata[is.na(alldata$Glucose),]$ID))
length(unique(alldata[is.na(alldata$HbA1c) | is.na(alldata$Glucose),]$ID))

# Find number with HbA1c > 184 mmol/mol (maximum detectable value)
length(unique(alldata[alldata$HbA1c > 184,]$ID))

# Find number whose only measurement is taken after diabetes (any) diagnosis
length(unique(alldata[alldata$diabetes == 1,]$ID))

# Find number whose only measurement was taken during gestation or postpartum period of pregnancy
length(unique(alldata[alldata$pregnant == 1,]$ID))

length(unique(c(alldata[alldata$diabetes == 1,]$ID, alldata[alldata$pregnant == 1,]$ID)))

## Apply exclusion criteria for HbA1c and glucose analyses

# Exclude those with missing HbA1c or random glucose measurements
hba1cdata <- alldata[!is.na(alldata$HbA1c) & !is.na(alldata$Glucose),]

# Exclude those with uncorrected HbA1c > 184 mmol/mol (maximum detectable value)
hba1cdata <- hba1cdata[hba1cdata$HbA1c <= 184,]

# Exclude individuals who visited after diabetes (any) diagnosis
hba1cdata <- hba1cdata[hba1cdata$diabetes == 0,]

# Exclude individuals who visited during gestation or postpartum period of pregnancy
hba1cdata <- hba1cdata[hba1cdata$pregnant == 0,]

length(unique(hba1cdata$ID))

# Recode the male homozygote genotypes from 2 to 1 for ease of the regression models - to keep as continuous variable (makes no difference)
hba1cdata[hba1cdata$Sex == "Male" & !is.na(hba1cdata$`23:154536002:C:T`) & hba1cdata$`23:154536002:C:T` == 2,]$`23:154536002:C:T` <- 1
hba1cdata[hba1cdata$Sex == "Male" & !is.na(hba1cdata$`23:154534419:G:A`) & hba1cdata$`23:154534419:G:A` == 2,]$`23:154534419:G:A` <- 1

hba1cdata$genotype <- NA
hba1cdata[hba1cdata$ethnic_group == "Black",]$genotype <- hba1cdata[hba1cdata$ethnic_group == "Black",]$`23:154536002:C:T`
hba1cdata[hba1cdata$ethnic_group == "Asian",]$genotype <- hba1cdata[hba1cdata$ethnic_group == "Asian",]$`23:154534419:G:A`
hba1cdata$ethnic_group <- factor(hba1cdata$ethnic_group, levels = c("White", "Black", "Asian", "Mixed", "Other"))

## Plot HbA1c and RG violin plots for genotypes for each ethnicity/sex group

# Generate genotype and sample size labels required for plots
hba1cdata$label <- hba1cdata$sub_label <- hba1cdata$label2 <- hba1cdata$sub_label2 <- hba1cdata$N <- hba1cdata$sub_N <- NA
hba1cdata$Sex <- factor(hba1cdata$Sex, levels = c("Male", "Female"))
afrdata <- aisdata <- hba1cdata
afrdata <- afrdata[!is.na(afrdata$`23:154536002:C:T`),]
aisdata <- aisdata[!is.na(aisdata$`23:154534419:G:A`),]

tabdata <- data.frame(matrix(NA, nrow = 2*length(ethnicities)*5, ncol = 5))
colnames(tabdata) <- c("Variant", "Ethnicity", "Sex", "Genotype", "Status")
tabdata$Variant <- rep(c("23:154536002:C:T", "23:154534419:G:A"), each = length(ethnicities)*5)
tabdata$Ethnicity <- rep(ethnicities, each = 5)
tabdata$Sex <- c(rep("Male", 2), rep("Female", 3))
tabdata$Genotype <- c(0, 1, 0:2)
tabdata$Status <- rep(c("No mutations", "Hemizygote", "No mutations", "Heterozygote", "Homozygote"), times = length(ethnicities))

for (i in 1:nrow(tabdata)) {
  ethnicity <- tabdata$Ethnicity[i]
  sex <- tabdata$Sex[i]
  genotype <- tabdata$Genotype[i]
  status <- tabdata$Status[i]
  subdata <- hba1cdata[hba1cdata$ethnic_group == ethnicity & hba1cdata$Sex == sex & !is.na(hba1cdata[tabdata$Variant[i]]) & hba1cdata[tabdata$Variant[i]] == genotype,]
  N <- prettyNum(nrow(subdata), big.mark = ",", scientific = FALSE)
  if (N != 0) {
    if (tabdata$Variant[i] == "23:154536002:C:T") {
      afrdata[afrdata$ethnic_group == ethnicity & afrdata$Sex == sex & afrdata$`23:154536002:C:T` == genotype,]$label <- paste0(status, "\n N = ", N)
      afrdata[afrdata$ethnic_group == ethnicity & afrdata$Sex == sex & afrdata$`23:154536002:C:T` == genotype,]$label2 <- paste0(status, " (N = ", N, ")")
      afrdata[afrdata$ethnic_group == ethnicity & afrdata$Sex == sex & afrdata$`23:154536002:C:T` == genotype,]$N <- nrow(subdata)
    } else if (tabdata$Variant[i] == "23:154534419:G:A") {
      aisdata[aisdata$ethnic_group == ethnicity & aisdata$Sex == sex & aisdata$`23:154534419:G:A` == genotype,]$label <- paste0(status, "\n N = ", N)
      aisdata[aisdata$ethnic_group == ethnicity & aisdata$Sex == sex & aisdata$`23:154534419:G:A` == genotype,]$label2 <- paste0(status, " (N = ", N, ")")
      aisdata[aisdata$ethnic_group == ethnicity & aisdata$Sex == sex & aisdata$`23:154534419:G:A` == genotype,]$N <- nrow(subdata)
    }
  }
}

afrdata <- afrdata[order(afrdata$Sex, afrdata$`23:154536002:C:T`),]
afrdata$label <- factor(afrdata$label, levels = unique(afrdata$label))
afrdata$label2 <- factor(afrdata$label2, levels = unique(afrdata$label2))
afrdata$sub_label <- factor(afrdata$sub_label, levels = unique(afrdata$sub_label))
afrdata$sub_label2 <- factor(afrdata$sub_label2, levels = unique(afrdata$sub_label2))

aisdata <- aisdata[order(aisdata$Sex, aisdata$`23:154536002:C:T`),]
aisdata$label <- factor(aisdata$label, levels = unique(aisdata$label))
aisdata$label2 <- factor(aisdata$label2, levels = unique(aisdata$label2))
aisdata$sub_label <- factor(aisdata$sub_label, levels = unique(aisdata$sub_label))
aisdata$sub_label2 <- factor(aisdata$sub_label2, levels = unique(aisdata$sub_label2))

# Create function to generate data and labels required and produce violin plots
violinplot <- function(variant, ethnicity, measure) {
  if (!(variant %in% c("AFR", "AIS")) | !(ethnicity %in% ethnicities) | !(measure %in% c("HbA1c", "Glucose"))) {
    stop("Missing or incorrect function input provided")
  }
  # Subset the appropriate dataframe
  if (variant == "AFR") {
    subdata <- afrdata
  } else {
    subdata <- aisdata
  }
  subdata <- subdata[subdata$ethnic_group == ethnicity,]
  subdata$measure <- unlist(subdata[colnames(subdata) == measure])
  # Remove groups where N < 5 (as these cannot be plotted)
  subdata <- subdata[subdata$N >= 5 | subdata$N == 0,]
  # Calculate the truncation limits based on MAD
  highlimit <- median(subdata$measure) + (5*mad(subdata$measure))
  lowlimit <- max(c(min(subdata$measure), (median(subdata$measure) - (5*mad(subdata$measure)))))
  nabove <- aggregate(measure ~ label, data = subdata[subdata$measure > highlimit,], length)
  maxvals <- aggregate(measure ~ label, data = subdata, max)
  # Create dataframe containing labels and variables defining location of labels for plots
  maxlabs <- merge(nabove, maxvals, by = "label")
  maxlabs$text <- paste("\U2191 N = ", maxlabs$measure.x, "\n Max = ", round(maxlabs$measure.y, digits = 1))
  maxlabs$measure <- highlimit
  maxlabs <- merge(maxlabs, unique(subdata[c("label", "Sex")]), by = "label")
  # Generate plot for specific measure
  if (measure == "HbA1c") {
    plotcol <- cbPalette[2]
  } else {
    plotcol <- cbPalette[3]
  }
  if (all(subdata$N >= 10)) {
    x <- ggplot(subdata, aes(x = factor(label), y = measure)) +
      geom_violin(fill = plotcol) + geom_boxplot(width = 0.1) +
      facet_grid( ~ Sex, scales = "free") + theme_bw() + theme(legend.position = "none") +
      geom_text(data = maxlabs, aes(label = text), size = 8, hjust = -0.2, vjust = 1) + xlab("Genotype") +
      theme(text = element_text(size = 32))
  } else {
    subdata2 <- subdata
    subdata2[subdata2$N < 10,]$measure <- 1000
    x <- ggplot(subdata[subdata$N >= 10,], aes(x = factor(label), y = measure)) +
      geom_violin(data = subdata2, fill = plotcol) + geom_boxplot(width = 0.1) +
      geom_dotplot(data = subdata[subdata$N < 10,], binaxis = "y", stackdir = "center", dotsize = 0.5, stackratio = 1.5) +
      facet_grid( ~ Sex, scales = "free") + theme_bw() + theme(legend.position = "none") +
      geom_text(data = maxlabs, aes(label = text), size = 8, hjust = -0.2, vjust = 1) + xlab("Genotype") +
      theme(text = element_text(size = 32))
  }
  if (measure == "HbA1c") {
    x <- x + ylab("HbA1c (mmol/mol)") + scale_y_continuous(sec.axis = sec_axis(trans = ~./10.929 + 2.15, name = "HbA1c (%)"), limits = c(lowlimit, highlimit))
  } else {
    x <- x + ylab("Random glucose (mmol/L)") + ylim(lowlimit, highlimit)
  }
  return(x)
}

# African variant (23:154536002:C:T)
png("HbA1c_23:154536002:C:T_black_violin_plot.png", 1400, 850)
violinplot(variant = "AFR", ethnicity = "Black", measure = "HbA1c")
dev.off()

png("RG_23:154536002:C:T_black_violin_plot.png", 1400, 850)
violinplot(variant = "AFR", ethnicity = "Black", measure = "Glucose")
dev.off()

# Asian variant (23:154534419:G:A)
png("HbA1c_23:154534419:G:A_asian_violin_plot.png", 1400, 850)
violinplot(variant = "AIS", ethnicity = "Asian", measure = "HbA1c")
dev.off()

png("RG_23:154534419:G:A_asian_violin_plot.png", 1400, 850)
violinplot(variant = "AIS", ethnicity = "Asian", measure = "Glucose")
dev.off()

system("dx upload *_violin_plot.png")

### Age of T2D diagnoses made since 2011 analyses ###

# Generate binary variable to identify whether T2D was diagnosed after 2011
alldata$t2d_2011 <- 0
alldata[!is.na(alldata$t2d_date) & alldata$t2d_date >= as.Date("2011-1-1"),]$t2d_2011 <- 1

alldata$t2d_2011_GP_self <- 0
alldata[!is.na(alldata$t2d_date_GP_self_min) & alldata$t2d_date_GP_self_min >= as.Date("2011-1-1"),]$t2d_2011_GP_self <- 1

# Generate summary characteristics for cohort
nrow(alldata[alldata$t2d_anyage == 1,])

# Make dataframe of those with T2D only (diagnosed at any time)
alldata <- alldata[alldata$t2d_anyage == 1 & !is.na(alldata$t2d_age),]

alldata <- alldata[alldata$ethnic_group %in% c("Black", "Asian"),]
alldata$ethnic_group <- factor(alldata$ethnic_group, levels = c("Black", "Asian"))
alldata$genotype <- NA
alldata[alldata$ethnic_group == "Black",]$genotype <- alldata[alldata$ethnic_group == "Black",]$`23:154536002:C:T`
alldata[alldata$ethnic_group == "Asian",]$genotype <- alldata[alldata$ethnic_group == "Asian",]$`23:154534419:G:A`
alldata <- alldata[!is.na(alldata$genotype),]

# Restrict the dataframes to those diagnosed by GP or self-reported - and keep the minimum age of these two sources only!
alldata <- alldata[!is.na(alldata$t2d_age_GP_self_min),]

alldata$t2d_age <- alldata$t2d_age_GP_self_min
alldata$t2d_age_origin <- alldata$t2d_age_GP_self_origin

# Generate variable for origin of diagnosis age (centre 1, ..., centre X, GP)
alldata[alldata$t2d_age_origin == "Selfreport",]$t2d_age_origin <- alldata[alldata$t2d_age_origin == "Selfreport",]$centre

# Recode the male homozygote genotypes from 2 to 1 for ease of the regression models - to keep as continuous variable (makes no difference)
alldata[alldata$Sex == "Male" & alldata$genotype == 2,]$genotype <- 1

afterdata <- alldata[alldata$t2d_2011_GP_self == 1,]

length(unique(afterdata$ID))

## Generate box/violin plots of age of diagnosis in Black males

subdata <- afterdata[afterdata$ethnic_group == "Black" & afterdata$Sex == "Male",]
subdata$label <- NA
subdata[subdata$genotype == 0,]$label <- paste0("No mutations \n N = ", length(unique(subdata[subdata$genotype == 0,]$ID)))
subdata[subdata$genotype == 1,]$label <- paste0("Hemizygote \n N = ", length(unique(subdata[subdata$genotype == 1,]$ID)))
subdata <- subdata[order(subdata$genotype),]
subdata$label <- factor(subdata$label, levels = unique(subdata$label))

tabdata <- data.frame(matrix(NA, nrow = 7, ncol = 3))
colnames(tabdata) <- c("Age of type 2 \n diabetes \n diagnosis (years)", "No mutations", "Hemizygote")

tabdata[6:7,1] <- c("Mean (SD)", "Difference of \n means (95% CI)")

tabdata$`No mutations`[1] <- mean(subdata[subdata$genotype == 0,]$t2d_age)
tabdata$`No mutations`[2] <- sd(subdata[subdata$genotype == 0,]$t2d_age)

tabdata$Hemizygote[1] <- mean(subdata[subdata$genotype == 1,]$t2d_age)
tabdata$Hemizygote[2] <- sd(subdata[subdata$genotype == 1,]$t2d_age)

tabdata$`No mutations`[3] <- mean(subdata[subdata$genotype == 1,]$t2d_age) - mean(subdata[subdata$genotype == 0,]$t2d_age)
tabdata$`No mutations`[4] <- t.test(subdata[subdata$genotype == 1,]$t2d_age, subdata[subdata$genotype == 0,]$t2d_age, var.equal = TRUE)$conf.int[1]
tabdata$`No mutations`[5] <- t.test(subdata[subdata$genotype == 1,]$t2d_age, subdata[subdata$genotype == 0,]$t2d_age, var.equal = TRUE)$conf.int[2]

tabdata[,2:3] <- round(tabdata[,2:3], digits = 1)
tabdata$`No mutations`[6] <- paste0(tabdata$`No mutations`[1], " (", tabdata$`No mutations`[2], ")")
tabdata$Hemizygote[6] <- paste0(tabdata$Hemizygote[1], ".0 (", tabdata$Hemizygote[2], ")")
tabdata$`No mutations`[7] <- paste0(tabdata$`No mutations`[3], " (", tabdata$`No mutations`[4], ", ", tabdata$`No mutations`[5], ")")
tabdata <- tabdata[6:7,]

tabdata[3,] <- colnames(tabdata)

subdata2 <- subdata
subdata2[subdata2$genotype == 1,]$t2d_age <- 1000

png("T2D_age_23:154536002:C:T_black_violin_plot.png", 1600, 1400)
p1 <- ggplot(subdata[subdata$genotype == 0,], aes(x = factor(label), y = t2d_age)) +
  geom_violin(data = subdata2, fill = cbPalette[4]) + geom_boxplot(width = 0.1) +
  geom_dotplot(data = subdata[subdata$genotype == 1,], binaxis = "y", stackdir = "center", dotsize = 0.5, stackratio = 1.5) +
  theme_bw() + theme(legend.position = "none") + xlab("Genotype") + ylab("Age of type 2 diabetes diagnosis (years)") +
  ylim(c(range(subdata$t2d_age))) +
  theme(text = element_text(size = 32), plot.margin = unit(c(0.1, 0.1, 0.1, 2.6), "inches"))

pc <- tableGrob(tabdata[3,], rows = NULL, cols = NULL, theme = ttheme_default(base_size = 25, padding = unit(c(1, 10), "mm"), core = list(bg_params = list(fill = "grey"), fg_params = list(fontface = "bold"))))
p2 <- tableGrob(tabdata[1,], rows = NULL, cols = NULL, theme = ttheme_default(base_size = 25, padding = unit(c(1, 15), "mm")))
p3 <- tableGrob(tabdata[2,1:2], rows = NULL, cols = NULL, theme = ttheme_default(base_size = 25, padding = unit(c(1, 10), "mm")))
jn <- gtable_combine(pc, p2, along = 2)
jn <- gtable_combine(jn, p3, along = 2)
jn$widths <- c(0.38, 1, 1)
jn$layout[c(14,16) , "r"] <- 3
grid.arrange(p1, jn, nrow = 2, heights = c(2.5, 0.7))
dev.off()

system("dx upload T2D_age_*_plot.png")

# Apply Mann-Whitney-U test to Black males
wilcox.test(t2d_age ~ genotype, data = subdata)

## Compare age of T2D diagnosis across genotypes for each variant/ethnicity grouping

se <- function(x) sqrt(var(x, na.rm = TRUE)/length(x[!is.na(x)]))

# Find the mean age of T2D diagnosis for each group and compare across genotypes for each variant - for diagnosis since 2011
tabdata <- data.frame(matrix(NA, nrow = 2*5, ncol = 10))
colnames(tabdata) <- c("Ethnicity", "Sex", "Genotype", "Status", "N", "Model.Beta", "Model.SE", "Model.LCI", "Model.UCI", "Model.P")
tabdata$Ethnicity <- rep(c("Black", "Asian"), each = 5)
tabdata$Sex <- c(rep("Male", 2), rep("Female", 3))
tabdata$Genotype <- c(0, 1, 0:2)
tabdata$Status <- c("No mutations", "Hemizygote", "No mutations", "Heterozygote", "Homozygote")
for (i in 1:nrow(tabdata)) {
  subdata <- afterdata[afterdata$ethnic_group == tabdata$Ethnicity[i],]
  subdata <- subdata[subdata$Sex == tabdata$Sex[i],]
  if (tabdata$Genotype[i] == 0) {
    if (length(unique(subdata$genotype)) > 1) {
      z <- lm(t2d_age ~ genotype + t2d_age_origin + ethnicity, data = subdata)
      tabdata$Model.Beta[i] <- summary(z)$coefficients["genotype", "Estimate"]
      tabdata$Model.SE[i] <- summary(z)$coefficients["genotype", "Std. Error"]
      tabdata$Model.LCI[i] <- tabdata$Model.Beta[i] - (1.96*tabdata$Model.SE[i])
      tabdata$Model.UCI[i] <- tabdata$Model.Beta[i] + (1.96*tabdata$Model.SE[i])
      tabdata$Model.P[i] <- summary(z)$coefficients["genotype", "Pr(>|t|)"]
    }
  }
  subdata <- subdata[subdata$genotype == tabdata$Genotype[i],]
  tabdata$N[i] <- length(unique(subdata$ID))
}

write.table(tabdata, file = "sum_model_T2D_age_after2011.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
system("dx upload sum_model_T2D_age_after2011.txt")

# Generate model diagnostic plots for the regression models
mod_m_age_t2d <- lm(t2d_age ~ genotype + t2d_age_origin + ethnicity, data = afterdata[afterdata$ethnic_group == "Black" & afterdata$Sex == "Male",])
png("UKB_black_male_diagnostic_plots.png", 800, 800)
par(mfrow = c(2, 2))
plot(mod_m_age_t2d, id.n = 0)
dev.off()

mod_f_age_t2d <- lm(t2d_age ~ genotype + t2d_age_origin + ethnicity, data = afterdata[afterdata$ethnic_group == "Black" & afterdata$Sex == "Female",])
png("UKB_black_female_diagnostic_plots.png", 800, 800)
par(mfrow = c(2, 2))
plot(mod_f_age_t2d, id.n = 0)
dev.off()

mod_m_age_t2d <- lm(t2d_age ~ genotype + t2d_age_origin + ethnicity, data = afterdata[afterdata$ethnic_group == "Asian" & afterdata$Sex == "Male",])
png("UKB_asian_male_diagnostic_plots.png", 800, 800)
par(mfrow = c(2, 2))
plot(mod_m_age_t2d, id.n = 0)
dev.off()

mod_f_age_t2d <- lm(t2d_age ~ genotype + t2d_age_origin + ethnicity, data = afterdata[afterdata$ethnic_group == "Asian" & afterdata$Sex == "Female",])
png("UKB_asian_female_diagnostic_plots.png", 800, 800)
par(mfrow = c(2, 2))
plot(mod_f_age_t2d, id.n = 0)
dev.off()

system("dx upload UKB_*_diagnostic_plots.png")

## Conduct meta-analysis across the different cohorts/variants/ethnicities groupings

# Read in Genes & Health model results for age of T2D diagnosis after 2011
gnhdata <- data.frame(Study = "G&H", Ethnicity = "South Asian", Sex = c("Male", "Female"),
                      N.0 = c(3621, 3490), N.1 = c(19, 57), N.2 = c(NA, 3),
                      Model.Beta = c(2.948304, 0.484190),
                      Model.SE = c(2.5932410, 1.4448189))
tabdata$Study <- "UKB"
metadata <- tabdata[!is.na(tabdata$Model.Beta) & tabdata$Ethnicity %in% c("Black", "Asian"),]
metadata <- rbind(metadata[c("Study", "Ethnicity", "Sex", "Model.Beta", "Model.SE")], gnhdata[c("Study", "Ethnicity", "Sex", "Model.Beta", "Model.SE")])

# Meta-analyse the model output (difference in age of T2D diagnosis per additional allele (females) / carriers vs non-carriers (males))
myresults <- data.frame(matrix(NA, nrow = 2, ncol = 9))
colnames(myresults) <- c("Sex", "Beta", "SE", "P", "Q", "I2", "hetP", "df", "tau2")
myresults$Sex <- c("Male", "Female")
for (i in 1:nrow(myresults)) {
  subdata <- metadata[metadata$Sex == myresults$Sex[i],]
  meta <- rma(yi = Model.Beta, sei = Model.SE, data = subdata, method = "REML")
  myresults$Beta[i] <- meta$beta[,1]
  myresults$SE[i] <- meta$se
  myresults$P[i] <- meta$pval
  myresults$Q[i] <- meta$QE
  myresults$I2[i] <- meta$I2
  myresults$hetP[i] <- meta$QEp
  myresults$df[i] <- meta$k - 1
  myresults$tau2[i] <- meta$tau2
  if (myresults$Sex[i] == "Male") {
    metam <- meta
  } else {
    metaf <- meta
  }
}

write.table(myresults, file = "meta_T2D_age_after2011.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
system("dx upload meta_T2D_age_after2011.txt")

tabdata <- tabdata[tabdata$Ethnicity %in% c("Black", "Asian"),]
tabdata$N.0 <- tabdata$N.1 <- tabdata$N.2 <- NA
for (ethnicity in c("Black", "Asian")) {
  for (sex in c("Male", "Female")) {
    tabdata[tabdata$Ethnicity == ethnicity & tabdata$Sex == sex,]$N.0 <- tabdata[tabdata$Ethnicity == ethnicity & tabdata$Sex == sex & tabdata$Genotype == 0,]$N
    tabdata[tabdata$Ethnicity == ethnicity & tabdata$Sex == sex,]$N.1 <- tabdata[tabdata$Ethnicity == ethnicity & tabdata$Sex == sex & tabdata$Genotype == 1,]$N
    if (sex == "Female") {
      tabdata[tabdata$Ethnicity == ethnicity & tabdata$Sex == sex,]$N.2 <- tabdata[tabdata$Ethnicity == ethnicity & tabdata$Sex == sex & tabdata$Genotype == 2,]$N
    }
  }
}
tabdata <- tabdata[!is.na(tabdata$Model.Beta),]
tabdata <- rbind(tabdata[c("Study", "Ethnicity", "Sex", "N.0", "N.1", "N.2", "Model.Beta", "Model.SE")], gnhdata[c("Study", "Ethnicity", "Sex", "N.0", "N.1", "N.2", "Model.Beta", "Model.SE")])
tabdata$label <- paste0(tabdata$Study, " ", tabdata$Ethnicity)
tabdata$N.0 <- prettyNum(tabdata$N.0, big.mark = ",", scientific = FALSE)

tabdata[tabdata$N.1 < 5 & tabdata$N.1 != 0,]$N.1 <- "<5"
tabdata[!is.na(tabdata$N.2) & tabdata$N.2 < 5 & tabdata$N.2 != 0,]$N.2 <- "<5"

plotdata <- tabdata[tabdata$Sex == "Male",]
png("forestplot_meta_T2D_age_male_after2011.png", 1100, 400)
par(mar = c(5.1,4.1,0,2.1))
metaplot <- forest(x = plotdata$Model.Beta, sei = plotdata$Model.SE, slab = unique(plotdata$label),
                   header = "Ethnicity", xlab = "Difference in age of T2D diagnosis (years)", ylim = c(-0.5, 6), xlim = c(-34, 26), cex = 1.7, col = cbPalette[c(6:7,4)],
                   ilab = cbind(plotdata$N.0, plotdata$N.1), ilab.xpos = c(-21, -12), digits = 1)
text(mean(metaplot$ilab.xpos[1:2]), 6, "No. of participants", cex = 1.7, font = 2)
text(metaplot$ilab.xpos[1:2], c(5.05, 5), c("No mutations", "Hemizygote"), cex = 1.7, font = 2)
segments(metaplot$ilab.xpos[1]-4.2, 5.5, metaplot$ilab.xpos[2]+4.2, 5.5)
text(-21, rev(metaplot$rows), paste0(plotdata$N.0), col = rep(c(cbPalette[c(6:7,4)]), times = 4), cex = 1.7)
text(-12, rev(metaplot$rows), paste0(plotdata$N.1), col = rep(c(cbPalette[c(6:7,4)]), times = 4), cex = 1.7)
addpoly(metam, row = 0, mlab = "", efac = 2, col = "darkgrey", border = "darkgrey")
text(metaplot$xlim[1], 0, pos = 4, cex = 1.7, bquote(paste(bold(.("Meta-analysis")))))
text(metaplot$xlim[1], -0.5, pos = 4, cex = 1.7, bquote(paste(I^2, " = ", .(fmtx(metam$I2, digits = 1)), "%, ", p[heterogeneity],
                                                              .(fmtp(metam$QEp, digits = 2, pname = " ", add0 = TRUE, sep = TRUE, equal = TRUE)))))
dev.off()

plotdata <- tabdata[tabdata$Sex == "Female",]
png("forestplot_meta_T2D_age_female_after2011.png", 1300, 400)
par(mar = c(5.1,4.1,0,2.1))
metaplot <- forest(x = plotdata$Model.Beta, sei = plotdata$Model.SE, slab = unique(plotdata$label),
                   header = "Ethnicity", xlab = "Difference in age of T2D diagnosis (years)", ylim = c(-0.5, 6), xlim = c(-36, 18), cex = 1.7, col = cbPalette[c(6:7,4)],
                   ilab = cbind(plotdata$N.0, plotdata$N.1, plotdata$N.2), ilab.xpos = c(-26, -19, -12), digits = 1)
text(metaplot$ilab.xpos[2], 6, "No. of participants", cex = 1.7, font = 2)
text(metaplot$ilab.xpos[1:3], c(5.05, 5, 5), c("No mutations", "Heterozygote", "Homozygote"), cex = 1.7, font = 2)
segments(metaplot$ilab.xpos[1]-3.2, 5.5, metaplot$ilab.xpos[3]+3.2, 5.5)
text(-26, rev(metaplot$rows), paste0(plotdata$N.0), col = rep(c(cbPalette[c(6:7,4)]), times = 4), cex = 1.7)
text(-19, rev(metaplot$rows), paste0(plotdata$N.1), col = rep(c(cbPalette[c(6:7,4)]), times = 4), cex = 1.7)
text(-12, rev(metaplot$rows), paste0(plotdata$N.2), col = rep(c(cbPalette[c(6:7,4)]), times = 4), cex = 1.7)
addpoly(metaf, row = 0, mlab = "", efac = 2, col = "darkgrey", border = "darkgrey")
text(metaplot$xlim[1], 0, pos = 4, cex = 1.7, bquote(paste(bold(.("Meta-analysis")))))
text(metaplot$xlim[1], -0.5, pos = 4, cex = 1.7, bquote(paste(I^2, " = ", .(fmtx(metaf$I2, digits = 1)), "%, ", p[heterogeneity],
                                                              .(fmtp(metaf$QEp, digits = 2, pname = " ", add0 = TRUE, sep = TRUE, equal = TRUE)))))
dev.off()

system("dx upload forestplot_meta_T2D_age_*_after2011.png")

