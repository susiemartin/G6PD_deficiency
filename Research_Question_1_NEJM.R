# Research Question 1
# Find the prevalence of two G6PD variants of interest in the UK Biobank
# 3 January 2024

# Load required packages
library(stringr)
library(ggplot2)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

system("dx cd /G6PD_diabetes_comps")

# Read in dataframe after applying exclusions for main analysis (2)
system("dx download G6PD_geno_pheno_data_clean.txt")
alldata <- read.csv("G6PD_geno_pheno_data_clean.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")

colnames(alldata)[colnames(alldata) == "X23.154536002.C.T"] <- "23:154536002:C:T"
colnames(alldata)[colnames(alldata) == "X23.154534419.G.A"] <- "23:154534419:G:A"

length(unique(alldata$ID))

# Read in the missingness for each variant and add into the MAF dataframe
system("dx download /PLINK_genotypes/G6PD_diabetes_comps/missingness_files/var_plink_chrX.lmiss")
mymiss <- read.table("var_plink_chrX.lmiss", header = TRUE, stringsAsFactors = FALSE)

# Calculate MAF for the IDs remaining in the cohort (e.g. those with pheno/ethnic data) - for self-reported ethnicity
ethnicities <- c("White", "Black", "Asian", "Mixed", "Other", "African", "Caribbean", "Other black", "South Asian", "East Asian", "Other Asian")
mafdata <- data.frame(matrix(NA, nrow = 2*length(ethnicities), ncol = 3))
colnames(mafdata) <- c("Variant", "Ethnicity", "MAF")
mafdata$Variant <- rep(c("23:154536002:C:T", "23:154534419:G:A"), each = length(ethnicities))
mafdata$Ethnicity <- ethnicities
for (i in 1:nrow(mafdata)) {
  subdata <- alldata[!is.na(alldata[mafdata$Variant[i]]),]
  if (mafdata$Ethnicity[i] %in% c("White", "Black", "Asian", "Mixed", "Other")) {
    subdata <- subdata[subdata$ethnic_group == mafdata$Ethnicity[i],]
  } else {
    subdata <- subdata[!is.na(subdata$sub_ethnic_group) & subdata$sub_ethnic_group == mafdata$Ethnicity[i],]
  }
  mafdata$MAF[i] <- sum(subdata[mafdata$Variant[i]])/(2*nrow(subdata))
}
mafdata <- merge(mafdata, mymiss[c("SNP", "F_MISS")], by.x = "Variant", by.y = "SNP", all.x = TRUE)
mafdata <- mafdata[order(mafdata$Variant, decreasing = TRUE),]

write.table(mafdata, file = "MAF_missingness_selfreport_pheno_cohort_only.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
system("dx upload MAF_missingness_selfreport_pheno_cohort_only.txt")

# Calculate MAF for the IDs remaining in the cohort (e.g. those with pheno/ethnic data) - for genetic-derived ancestry
ancestries <- c("EUR", "AFR", "SAS", "EAS", "AMR")
mafdata <- data.frame(matrix(NA, nrow = 2*length(ancestries), ncol = 3))
colnames(mafdata) <- c("Variant", "Ancestry", "MAF")
mafdata$Variant <- rep(c("23:154536002:C:T", "23:154534419:G:A"), each = length(ancestries))
mafdata$Ancestry <- ancestries
for (i in 1:nrow(mafdata)) {
  subdata <- alldata[!is.na(alldata[mafdata$Variant[i]]),]
  subdata <- subdata[!is.na(subdata$genetic_ancestry) & subdata$genetic_ancestry == mafdata$Ancestry[i],]
  mafdata$MAF[i] <- sum(subdata[mafdata$Variant[i]])/(2*nrow(subdata))
}
mafdata <- merge(mafdata, mymiss[c("SNP", "F_MISS")], by.x = "Variant", by.y = "SNP", all.x = TRUE)
mafdata <- mafdata[order(mafdata$Variant, decreasing = TRUE),]

write.table(mafdata, file = "MAF_missingness_geneticPCA_pheno_cohort_only.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
system("dx upload MAF_missingness_geneticPCA_pheno_cohort_only.txt")

##### Determine the prevalence of carriers and of G6PD deficiency and/or haemolytic anaemia diagnoses in carriers #####

## Plot barplots of prevalence of carriers per genotype (>0) across each ethnicity - for each G6PD variant
plotdata <- data.frame(matrix(NA, nrow = 2*length(ethnicities)*3, ncol = 10))
colnames(plotdata) <- c("Variant", "Ethnicity", "Sex", "Genotype", "Status", "N", "Prev", "SD", "LCI", "UCI")
plotdata$Variant <- rep(c("23:154536002:C:T", "23:154534419:G:A"), each = length(ethnicities)*3)
plotdata$Ethnicity <- rep(ethnicities, each = 3)
plotdata$Sex <- c("Male", rep("Female", times = 2))
plotdata$Genotype <- c(2,1,2)
plotdata$Status <- c("Hemizygote", "Heterozygote", "Homozygote")
for (i in 1:nrow(plotdata)) {
  subdata <- alldata[alldata$Sex == plotdata$Sex[i],]
  if (plotdata$Ethnicity[i] %in% c("White", "Black", "Asian", "Mixed", "Other")) {
    subdata <- subdata[subdata$ethnic_group == plotdata$Ethnicity[i],]
  } else {
    subdata <- subdata[!is.na(subdata$sub_ethnic_group) & subdata$sub_ethnic_group == plotdata$Ethnicity[i],]
  }
  subdata <- subdata[!is.na(subdata[plotdata$Variant[i]]),]
  subdata2 <- subdata[subdata[plotdata$Variant[i]] == plotdata$Genotype[i],]
  plotdata$N[i] <- prettyNum(length(unique(subdata2$ID)), big.mark = ",", scientific = FALSE)
  p <- length(unique(subdata2$ID))/length(unique(subdata$ID))
  plotdata$Prev[i] <- p*100
  plotdata$SD[i] <- sqrt((p*(1-p))/length(unique(subdata$ID)))*100
  plotdata$LCI[i] <- (p - (1.96*plotdata$SD[i]))*100
  plotdata$UCI[i] <- (p + (1.96*plotdata$SD[i]))*100
}

plotdata$Ethnicity <- factor(plotdata$Ethnicity, levels = ethnicities)
plotdata$Sex <- factor(plotdata$Sex, levels = c("Male", "Female"))
plotdata[plotdata$Prev == 0,]$Prev <- NA
plotdata <- plotdata[order(plotdata$Status),]
plotdata$label <- paste0(plotdata$Status, "\n N = ", plotdata$N)
plotdata$label <- factor(plotdata$label, levels = unique(plotdata$label))
plotdata$bar_label <- paste0(round(plotdata$Prev, 2), "%")
plotdata[is.na(plotdata$Prev),]$bar_label <- NA

# African variant (23:154536002:C:T)
png("barplot_carrier_23:154536002:C:T_main.png", 2000, 700)
ggplot(plotdata[plotdata$Variant == "23:154536002:C:T" & plotdata$Ethnicity %in% c("White", "Black", "Asian", "Mixed", "Other"),], aes(x = label, y = Prev, fill = Sex)) + 
  geom_bar(position = "stack", stat = "identity", lwd = 1.5) + geom_errorbar(aes(ymin = Prev - SD, ymax = Prev + SD), width = 0.4, size = 2) +
  facet_grid(~ Ethnicity, scales = "free", space = "free") + xlab("") + ylab("Prevalence of genotype (%)") +
  theme_bw() + theme(legend.position = "top", legend.title = element_blank()) + scale_fill_manual(values = cbPalette[6:7]) +
  geom_text(aes(y = Prev + SD, label = bar_label), size = 9, vjust = -0.5) + scale_y_continuous(limits = c(0,27), breaks = c(0,5,10,15,20,25)) +
  theme(axis.title = element_text(size = 30), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 30), strip.text = element_text(size = 30), legend.text = element_text(size = 30))
dev.off()

png("barplot_carrier_23:154536002:C:T_sub.png", 1500, 700)
ggplot(plotdata[plotdata$Variant == "23:154536002:C:T" & plotdata$Ethnicity %in% c("African", "Caribbean", "Other black"),], aes(x = label, y = Prev, fill = Sex)) + 
  geom_bar(position = "stack", stat = "identity", lwd = 1.5) + geom_errorbar(aes(ymin = Prev - SD, ymax = Prev + SD), width = 0.4, size = 2) +
  facet_grid(~ Ethnicity, scales = "free", space = "free") + xlab("") + ylab("Prevalence of genotype (%)") +
  theme_bw() + theme(legend.position = "top", legend.title = element_blank()) + scale_fill_manual(values = cbPalette[6:7]) +
  geom_text(aes(y = Prev + SD, label = bar_label), size = 9, vjust = -0.5) + scale_y_continuous(limits = c(0,28), breaks = c(0,5,10,15,20,25)) +
  theme(axis.title = element_text(size = 30), axis.text = element_text(size = 30), strip.text = element_text(size = 30), legend.text = element_text(size = 30))
dev.off()

# Asian variant (23:154534419:G:A)
png("barplot_carrier_23:154534419:G:A_main.png", 2000, 700)
ggplot(plotdata[plotdata$Variant == "23:154534419:G:A" & plotdata$Ethnicity %in% c("White", "Black", "Asian", "Mixed", "Other"),], aes(x = label, y = Prev, fill = Sex)) + 
  geom_bar(position = "stack", stat = "identity", lwd = 1.5) + geom_errorbar(aes(ymin = Prev - SD, ymax = Prev + SD), width = 0.4, size = 2) +
  facet_grid(~ Ethnicity, scales = "free", space = "free") + xlab("") + ylab("Prevalence of genotype (%)") +
  theme_bw() + theme(legend.position = "top", legend.title = element_blank()) + scale_fill_manual(values = cbPalette[6:7]) +
  geom_text(aes(y = Prev + SD, label = bar_label), size = 9, vjust = -0.4) + scale_y_continuous(limits = c(0,3), breaks = c(0,0.5,1,1.5,2,2.5)) +
  theme(axis.title = element_text(size = 30), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 30), strip.text = element_text(size = 30), legend.text = element_text(size = 30))
dev.off()

png("barplot_carrier_23:154534419:G:A_sub.png", 1500, 700)
ggplot(plotdata[plotdata$Variant == "23:154534419:G:A" & plotdata$Ethnicity %in% c("South Asian", "East Asian", "Other Asian"),], aes(x = label, y = Prev, fill = Sex)) + 
  geom_bar(position = "stack", stat = "identity", lwd = 1.5) + geom_errorbar(aes(ymin = Prev - SD, ymax = Prev + SD), width = 0.4, size = 2) +
  facet_grid(~ Ethnicity, scales = "free", space = "free") + xlab("") + ylab("Prevalence of genotype (%)") +
  theme_bw() + theme(legend.position = "top", legend.title = element_blank()) + scale_fill_manual(values = cbPalette[6:7]) +
  geom_text(aes(y = Prev + SD, label = bar_label), size = 9, vjust = -0.5) + scale_y_continuous(limits = c(0,2.2), breaks = c(0,0.5,1,1.5,2)) +
  theme(axis.title = element_text(size = 30), axis.text = element_text(size = 30), strip.text = element_text(size = 30), legend.text = element_text(size = 30))
dev.off()

system("dx upload barplot_carrier_*.png")

# Find the number of carriers of the Mediterranean variant in self-reported east Asians and other Asian backgrounds
aisdata <- data.frame(matrix(NA, nrow = 5*5, ncol = 5))
colnames(aisdata) <- c("Ethnicity", "Sex", "Genotype", "Total_N", "N_with_T2D")
aisdata$Ethnicity <- rep(c("Chinese", "Any other Asian background", "Bangladeshi", "Indian", "Pakistani"), each = 5)
aisdata$Sex <- c(rep("Male", 2), rep("Female", 3))
aisdata$Genotype <- c(0,2,0:2)
for (i in 1:nrow(aisdata)) {
  subdata <- alldata[alldata$ethnicity == aisdata$Ethnicity[i] & alldata$Sex == aisdata$Sex[i] & alldata$`23:154534419:G:A` == aisdata$Genotype[i],]
  aisdata$Total_N[i] <- length(unique(subdata$ID))
  aisdata$N_with_T2D[i] <- length(unique(subdata[subdata$t2d_anyage == 1,]$ID))
}
write.table(aisdata, file = "AIS_variant_sub_ethnicities.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
system("dx upload AIS_variant_sub_ethnicities.txt")

## Plot barplots of prevalence of diagnoses clustered by genotype - for corresponding ethnicity only for each G6PD variant
plotdata <- data.frame(matrix(NA, nrow = 2*length(ethnicities)*5, ncol = 10))
colnames(plotdata) <- c("Variant", "Ethnicity", "Sex", "Genotype", "Status", "N", "Percent_diag", "SD", "LCI", "UCI")
plotdata$Variant <- rep(c("23:154536002:C:T", "23:154534419:G:A"), each = length(ethnicities)*5)
plotdata$Ethnicity <- rep(ethnicities, each = 5)
plotdata$Sex <- c(rep("Male", times = 2), rep("Female", times = 3))
plotdata$Genotype <- c(0,2,0,1,2)
plotdata$Status <- c("No mutations", "Hemizygote", "No mutations", "Heterozygote", "Homozygote")
for (i in 1:nrow(plotdata)) {
  subdata <- alldata[alldata$Sex == plotdata$Sex[i],]
  if (plotdata$Ethnicity[i] %in% c("White", "Black", "Asian", "Mixed", "Other")) {
    subdata <- subdata[subdata$ethnic_group == plotdata$Ethnicity[i],]
  } else {
    subdata <- subdata[!is.na(subdata$sub_ethnic_group) & subdata$sub_ethnic_group == plotdata$Ethnicity[i],]
  }
  subdata <- subdata[!is.na(subdata[plotdata$Variant[i]]) & subdata[plotdata$Variant[i]] == plotdata$Genotype[i],]
  plotdata$N[i] <- nrow(subdata)
  p <- nrow(subdata[subdata$g6pd_or_haem == 1,])/plotdata$N[i]
  plotdata$Percent_diag[i] <- p*100
  plotdata$SD[i] <- sqrt((p*(1-p))/plotdata$N[i])*100
  plotdata$LCI[i] <- (p - (1.96*plotdata$SD[i]))*100
  plotdata$UCI[i] <- (p + (1.96*plotdata$SD[i]))*100
}

plotdata$Ethnicity <- factor(plotdata$Ethnicity, levels = ethnicities)
plotdata$Sex <- factor(plotdata$Sex, levels = c("Male", "Female"))
plotdata$Status <- factor(plotdata$Status, levels = unique(plotdata$Status))
plotdata <- plotdata[order(plotdata$Sex, plotdata$Status),]
plotdata$label <- paste0(plotdata$Status, "\n N = ", prettyNum(plotdata$N, big.mark = ",", scientific = FALSE))
plotdata$label <- factor(plotdata$label, levels = unique(plotdata$label))

# African variant (23:154536002:C:T)
png("barplot_carrier_diag_23:154536002:C:T_main.png", 1400, 800)
ggplot(plotdata[plotdata$Variant == "23:154536002:C:T" & plotdata$Ethnicity == "Black",], aes(x = label, y = Percent_diag, fill = Sex)) +
  geom_bar(stat = "identity") + geom_errorbar(aes(ymin = Percent_diag - SD, ymax = Percent_diag + SD), width = 0.4, size = 2) +
  facet_grid(~ Sex, scales = "free") + xlab("") + ylab("Percentage with haemolytic anaemia diagnosis (%)") +
  theme_bw() + theme(legend.position = "") + scale_fill_manual(values = cbPalette[6:7]) +
  theme(axis.title = element_text(size = 30), axis.text = element_text(size = 30), strip.text = element_text(size = 30))
dev.off()

png("barplot_carrier_diag_23:154536002:C:T_sub.png", 2000, 700)
ggplot(plotdata[plotdata$Variant == "23:154536002:C:T" & plotdata$Ethnicity %in% c("African", "Caribbean", "Other black"),], aes(x = label, y = Percent_diag, fill = Sex)) +
  geom_bar(stat = "identity") + geom_errorbar(aes(ymin = Percent_diag - SD, ymax = Percent_diag + SD), width = 0.4, size = 2) +
  facet_grid(~ Ethnicity, scales = "free") + xlab("") + ylab("Percentage with haemolytic anaemia diagnosis (%)") +
  theme_bw() + theme(legend.position = "top", legend.title = element_blank()) + scale_fill_manual(values = cbPalette[6:7]) +
  theme(axis.title = element_text(size = 30), axis.text.x = element_text(size = 25), axis.text.y = element_text(size = 30), strip.text = element_text(size = 30), legend.text = element_text(size = 30))
dev.off()

# Asian variant (23:154534419:G:A)
png("barplot_carrier_diag_23:154534419:G:A_main.png", 1400, 800)
ggplot(plotdata[plotdata$Variant == "23:154534419:G:A" & plotdata$Ethnicity == "Asian",], aes(x = label, y = Percent_diag, fill = Sex)) +
  geom_bar(stat = "identity") + geom_errorbar(aes(ymin = Percent_diag - SD, ymax = Percent_diag + SD), width = 0.4, size = 2) +
  facet_grid(~ Sex, scales = "free") + xlab("") + ylab("Percentage with haemolytic anaemia diagnosis (%)") +
  theme_bw() + theme(legend.position = "") + scale_fill_manual(values = cbPalette[6:7]) +
  theme(axis.title = element_text(size = 30), axis.text = element_text(size = 30), strip.text = element_text(size = 30))
dev.off()

png("barplot_carrier_diag_23:154534419:G:A_sub.png", 2000, 700)
ggplot(plotdata[plotdata$Variant == "23:154534419:G:A" & plotdata$Ethnicity %in% c("South Asian", "East Asian", "Other Asian"),], aes(x = label, y = Percent_diag, fill = Sex)) +
  geom_bar(stat = "identity") + geom_errorbar(aes(ymin = Percent_diag - SD, ymax = Percent_diag + SD), width = 0.4, size = 2) +
  facet_grid(~ Ethnicity, scales = "free") + xlab("") + ylab("Percentage with haemolytic anaemia diagnosis (%)") +
  theme_bw() + theme(legend.position = "top", legend.title = element_blank()) + scale_fill_manual(values = cbPalette[6:7]) +
  theme(axis.title = element_text(size = 30), axis.text.x = element_text(size = 25), axis.text.y = element_text(size = 30), strip.text = element_text(size = 30), legend.text = element_text(size = 30))
dev.off()

system("dx upload barplot_carrier_diag_*.png")

# Generate tables of genotypes versus number of diagnoses for each condition
myresults <- data.frame(matrix(NA, nrow = 2*length(ethnicities)*2*3, ncol = 8))
colnames(myresults) <- c("Variant", "Ethnicity", "Sex", "Genotype", "N", "G6PD_deficiency", "Haemolytic_anaemia", "Both")
myresults$Variant <- rep(c("23:154536002:C:T", "23:154534419:G:A"), each = length(ethnicities)*2*3)
myresults$Ethnicity <- rep(ethnicities, each = 2*3)
myresults$Sex <- rep(c("Male", "Female"), each = 3)
myresults$Genotype <- 0:2
for (i in 1:nrow(myresults)) {
  subdata <- alldata[alldata$Sex == myresults$Sex[i],]
  if (myresults$Ethnicity[i] %in% c("White", "Black", "Asian", "Mixed", "Other")) {
    subdata <- subdata[subdata$ethnic_group == myresults$Ethnicity[i],]
  } else {
    subdata <- subdata[!is.na(subdata$sub_ethnic_group) & subdata$sub_ethnic_group == myresults$Ethnicity[i],]
  }
  subdata <- subdata[!is.na(subdata[myresults$Variant[i]]) & subdata[myresults$Variant[i]] == myresults$Genotype[i],]
  myresults$N[i] <- length(unique(subdata$ID))
  myresults$G6PD_deficiency[i] <- length(unique(subdata[subdata$g6pd == 1,]$ID))
  myresults$Haemolytic_anaemia[i] <- length(unique(subdata[subdata$haem == 1,]$ID))
  myresults$Both[i] <- length(unique(subdata[subdata$g6pd_and_haem == 1,]$ID))
}
myresults <- myresults[!(myresults$Sex == "Male" & myresults$Genotype == 1),]

write.table(myresults, file = "genotype_g6pd_haem.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
system("dx upload genotype_g6pd_haem.txt")

# Identify proportion of diagnosed non-carriers (genotype = 0) have specifically G6PD deficiency diagnoses
subdata <- myresults[myresults$Variant == "23:154536002:C:T" & myresults$Genotype == 0 & myresults$Ethnicity %in% c("White", "Black", "Asian", "Mixed", "Other"),]
sum(subdata$G6PD_deficiency)/sum(subdata$Haemolytic_anaemia)*100

subdata <- myresults[myresults$Variant == "23:154534419:G:A" & myresults$Genotype == 0 & myresults$Ethnicity %in% c("White", "Black", "Asian", "Mixed", "Other"),]
sum(subdata$G6PD_deficiency)/sum(subdata$Haemolytic_anaemia)*100

