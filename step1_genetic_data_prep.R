# Script 1 - extract genetic data for G6PD variants of interest from the UK Biobank WES data (after running PLINK)
# 3 January 2024

# Load required packages
library(tidyr)
library(stringr)

system("dx cd /G6PD_diabetes_comps")

##### Read in and prepare genotype data for two G6PD variants #####

# Read in PLINK G6PD variant genotyping results file
system("dx download /PLINK_genotypes/G6PD_diabetes_comps/genotype_files/var_plink_chrX_genotype.txt")

# Read in chromosome X-specific genotype file
mygeno <- read.table(paste0("var_plink_chrX_genotype.txt"), header = TRUE, stringsAsFactors = FALSE)
mygeno <- mygeno[, -which(names(mygeno) %in% c("FID", "PAT", "MAT", "SEX", "PHENOTYPE"))]
mygeno <- setNames(data.frame(t(mygeno[,-1])), mygeno[,1])
mygeno$variant <- rownames(mygeno)
mygeno <- mygeno[, c(ncol(mygeno), 1:(ncol(mygeno)-1))]

# Read in the missingness for each variant and add into the MAF dataframe
system("dx download /PLINK_genotypes/G6PD_diabetes_comps/missingness_files/var_plink_chrX.lmiss")
mymiss <- read.table("var_plink_chrX.lmiss", header = TRUE, stringsAsFactors = FALSE)

# Alter the variant ID in genotype dataframe to match MAF dataframes and create count_allele column
mygeno$variant <- str_remove_all(mygeno$variant, "X")
mygeno$variant <- str_replace_all(mygeno$variant, "[.]", ":")
mygeno <- mygeno %>% separate(variant, c("variant", "count_allele"), "_")

# Read in dataframe indicating those with the lowest p-value to add associated marker column
file_name <- "HbA1c_no_diag_diab"
system(paste0("dx download /variant_lists/G6PD/var_for_gene_masks_", file_name, ".txt"))
genedata <- read.table(paste0("var_for_gene_masks_", file_name, ".txt"), header = TRUE, stringsAsFactors = FALSE)
alldata <- merge(genedata[c("variant", "gene", "mask", "EffectAllele", "lowest_P")], mygeno, by = "variant", all.y = TRUE)

# Split up the gene column to be by gene name and ENST ID
alldata$gene_name <- sapply(strsplit(alldata$gene, "_"), "[", 1)
alldata$ENST_ID <- sapply(strsplit(alldata$gene, "_"), "[", 2)
alldata <- alldata[, -which(names(alldata) %in% c("gene"))]
colnames(alldata)[colnames(alldata) == "gene_name"] <- "gene"
alldata <- alldata[,c(1, (ncol(alldata)-1):ncol(alldata), 2:(ncol(alldata)-2))]

# Determine whether the allele corresponding to genotype counts is the same as the effect allele
# - if not - inverse the genotype counts to represent the effect allele
for (i in 1:nrow(alldata)) {
  if (!is.na(alldata$EffectAllele[i]) & !is.na(alldata$count_allele[i]) & (alldata$EffectAllele[i] != alldata$count_allele[i])) {
    alldata[i, grepl("[0-9]", colnames(alldata))] <- abs(alldata[i, grepl("[0-9]", colnames(alldata))]-2)
    alldata$count_allele[i] <- alldata$EffectAllele[i]
  }
}

# Transpose the above genotype dataframe
idcols <- colnames(alldata)[grepl("[0-9]", colnames(alldata))]
alldata <- alldata[c("variant", idcols)]
alldata <- unique(alldata)
varlist <- unique(alldata$variant)

alldata <- setNames(data.frame(t(alldata[,-1])), alldata[,1])
alldata$ID <- rownames(alldata)
# Keep only those rows corresponding to ID's genotypes
alldata <- alldata[grepl("[0-9]", alldata$ID),]
# Remove the 'X' from the ID names
alldata$ID <- str_remove_all(alldata$ID, "X")

# Save the genetic data
write.table(alldata, file = "G6PD_genetic_data.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
system("dx upload G6PD_genetic_data.txt")

