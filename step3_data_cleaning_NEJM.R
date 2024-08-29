# Script 3 - merge genotype and phenotype data together and apply data cleaning and corrections
# 22 February 2024

system("dx cd /G6PD_diabetes_comps")

### Read in the genotype and phenotype data to merge together ###

system("dx download G6PD_genetic_data.txt")
gendata <- read.table("G6PD_genetic_data.txt", header = TRUE, stringsAsFactors = FALSE)

system("dx download G6PD_phenotype_data.txt")
phenodata <- read.csv("G6PD_phenotype_data.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")

alldata <- merge(gendata, phenodata, by.x = "ID", by.y = "Participant.ID")

colnames(alldata)[colnames(alldata) == "X23.154534419.G.A"] <- "23:154534419:G:A"
colnames(alldata)[colnames(alldata) == "X23.154536002.C.T"] <- "23:154536002:C:T"

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

# Generate binary variable to identify whether T2D was diagnosed after 2011
alldata$t2d_2011 <- 0
alldata[!is.na(alldata$t2d_date) & alldata$t2d_date >= as.Date("2011-1-1"),]$t2d_2011 <- 1

alldata$t2d_2011_self_GP <- 0
alldata[!is.na(alldata$t2d_date_GP_self_min) & alldata$t2d_date_GP_self_min >= as.Date("2011-1-1"),]$t2d_2011_self_GP <- 1

# Save this merged dataframe to take forward to analysis
write.table(alldata, file = "G6PD_geno_pheno_data_clean.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
system("dx upload G6PD_geno_pheno_data_clean.txt")

