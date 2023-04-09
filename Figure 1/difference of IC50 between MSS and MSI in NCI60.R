# 1. load R packages --------------------------------------------------------------------

rm(list = ls())
options(stringsAsFactors = FALSE)

library(xlsx)
library(TCGAbiolinks)
library(ggnewscale)
# BiocManager::install("escape")
library(escape)
# BiocManager::install("dittoSeq")
library(dittoSeq)
# BiocManager::install("SingleCellExperiment")
# library(SingleCellExperiment)
# install.packages("Seurat")
library(Seurat)
# install.packages("GSVA")
library(clusterProfiler)
library(org.Hs.eg.db)
library(httr)
library(RCurl)
library(rjson)
library(tidyr)
library(GSEABase)
library(GSVA)
library(openxlsx)
library(limma)
library(TCGAbiolinks)
# install.packages("gcookbook")

library(gcookbook)
library(ggplot2)
# install.packages("ggthemes")
library(ggthemes)
library(reshape2)
library(DOSE)
# BiocManager::install("DESeq2")

library(DESeq2)
library(enrichplot)
library(GSA)
library(ggpubr)
library(GEOquery)
library(dplyr)
library(ComplexHeatmap)
library(enrichplot)
library(pRRophetic)

# 2. read data--------------------------------------------------------------------

# IC50 data
ic50 <- read.csv("/mnt/DATA/ytj/NCI-60/IC50.csv" , header = T)

# Cancer Chemotherapy National Service Center (NSC from DPT) drug number annotation data
# that means it contains only chemotherapy drugs
nsc_anno <- read.xlsx("/mnt/DATA/ytj/NCI-60/DTP_NCI60_ZSCORE.xlsx", colNames = F )


# ms data from GDSC mainly used for annotation of msi status
ms <- read.xlsx(
  "/mnt/DATA/ytj/GDSC/Cell_Lines_Details.xlsx")

# only keep tumor tissue
gdsc_cell_line <- read.csv("/mnt/DATA/ytj/GDSC/model_list_20211124.csv",
                           header = TRUE)
gdsc_cell_line1 <- gdsc_cell_line[,c(17 , 18 , 39)]
colnames(gdsc_cell_line1)[3] <- "COSMIC.identifier"
ms <- merge(gdsc_cell_line1 , ms , by.y = "COSMIC.identifier")
ms <- ms[ms$tissue_status%in%c("Metastasis","Tumour"),]

# 3. data preprocessing -----------------------------------------------------------

# annotating IC50 data with NSC data
nsc_anno1 <- nsc_anno[-(1:6),]
colnames(nsc_anno1) <- nsc_anno1[1,]
nsc_anno1 <- nsc_anno1[-1,]
nsc <- nsc_anno1[,1:2]
colnames(nsc)[1] <- "NSC"

ic501 <- merge(ic50 , nsc , by.y = "NSC")

# nsc data can not cover all drugs from ic50 data
unique(ic501$`Drug name`)
length(unique(ic50$NSC))
length(unique(nsc$NSC))

unique(ic501$CELL_NAME)
unique(ms$Sample.Name)
unique(ic501$CELL_NAME)[!unique(ic501$CELL_NAME)%in%unique(ms$Sample.Name)]

# some cell line symbols are not quite the same in gdsc and ic50
unique(ms$`Microsatellite.instability.Status.(MSI)`)
ms1 <- ms[ms$`Microsatellite.instability.Status.(MSI)`%in%c("MSS/MSI-L","MSI-H"),]
unique(ic501$CELL_NAME)[!unique(ic501$CELL_NAME)%in%unique(ms1$Sample.Name)]

ms1[,4] <- gsub("SF268","SF-268",ms1[,4])
ms1[,4] <- gsub("SF295","SF-295",ms1[,4])
ms1[,4] <- gsub("T47D","T-47D",ms1[,4])
ms1[,4] <- gsub("MDA-MB-231","MDA-MB-231/ATCC",ms1[,4])
ms1[,4] <- gsub("RXF393","RXF 393",ms1[,4])
ms1[,4] <- gsub("UO31","UO-31",ms1[,4])
ms1[,4] <- gsub("HT-29","HT29",ms1[,4])
ms1[,4] <- gsub("SW620","SW-620",ms1[,4])
ms1[,4] <- gsub("TK10","TK-10",ms1[,4])
ms1[,4] <- gsub("SF539","SF-539",ms1[,4])
ms1[,4] <- gsub("COLO-205","COLO 205",ms1[,4])
ms1[,4] <- gsub("A549","A549/ATCC",ms1[,4])
ms1[,4] <- gsub("HCC2998","HCC-2998",ms1[,4])
ms1[,4] <- gsub("IGROV-1","IGROV1",ms1[,4])
ms1[,4] <- gsub("LOXIMVI","LOX IMVI",ms1[,4])
ms1[,4] <- gsub("Hs-578-T","HS 578T",ms1[,4])
ms1[,4] <- gsub("HL-60","HL-60(TB)",ms1[,4])
ms1[,4] <- gsub("A172","A-172/H.Fine",ms1[,4])
ms1[,4] <- gsub("U251","U-251/H.Fine",ms1[,4])
ms1[,4] <- gsub("LN-229","LN-229/H.Fine",ms1[,4])
ms1[,4] <- gsub("T98G","T98G/H.Fine",ms1[,4])
ms1[,4] <- gsub("DMS-273","DMS 273",ms1[,4])
ms1[,4] <- gsub("DMS-114","DMS 114",ms1[,4])
ms1[,4] <- gsub("SW1088","SW 1088",ms1[,4])
ms1[,4] <- gsub("SW1783","SW 1783",ms1[,4])
ms1[,4] <- gsub("MCF7","MCF7/ATCC",ms1[,4])
ms1[,4] <- gsub("ES-2","ES-2 CNCR",ms1[,4])
ms1[,4] <- gsub("RKO","RKO Waf1",ms1[,4])
unique(ic501$CELL_NAME)[!unique(ic50$CELL_NAME)%in%unique(ms1$Sample.Name)]

ms2 <- ms1[,c(4,12,13)]
colnames(ms2)[1] <- "CELL_NAME"

# annotation msi status
ic502 <- merge(ic501, ms2 , by.y = "CELL_NAME")

# check log transformation
range(ic502$AVERAGE)

ic502_log <- ic502
ic502_log[,13] <- 10^ic502_log[,13]
range(ic502_log$AVERAGE)
ic50_msi <- ic502_log[ic502_log$`Microsatellite.instability.Status.(MSI)`%in%"MSI-H",]
ic50_mss <- ic502_log[ic502_log$`Microsatellite.instability.Status.(MSI)`%in%"MSS/MSI-L",]

# intersection of tumors
unique(ic50_msi$`Cancer.Type.(matching.TCGA.label)`)
unique(ic50_mss$`Cancer.Type.(matching.TCGA.label)`)
unique(ic502_log$`Cancer.Type.(matching.TCGA.label)`)
unique(ic50_msi$`Cancer.Type.(matching.TCGA.label)`)[unique(ic50_msi$`Cancer.Type.(matching.TCGA.label)`)%in%unique(ic50_mss$`Cancer.Type.(matching.TCGA.label)`)]
cancer_ic50 <- unique(ic50_msi$`Cancer.Type.(matching.TCGA.label)`)
ic50_mss <- ic50_mss[ic50_mss$`Cancer.Type.(matching.TCGA.label)`%in%cancer_ic50,]

# intersection of drug
length(unique(ic50_msi$`Drug name`))
length(unique(ic50_mss$`Drug name`))

drug_ic50 <- unique(ic50_msi$`Drug name`)

# take a list of how many samples correspond to each tumor and drug, and discard those with only one sample
length_msi <- data.frame()
length_mss <- data.frame()

drug_seq_predict <- c("Cytarabine","Fluorouracil","Teniposide","Epirubicin","Fludarabine","Cisplatin","Mitomycin","Bleomycin","Gemcitabine","Camptothecin","Topotecan","Irinotecan","Mitoxantrone","Etoposide","Vinblastine","Vinorelbine","Docetaxel","Paclitaxel","Actinomycin D","Carmustine","Vincristine","Chlorambucil","Lomustine","Melphalan","Carboplatin","Daunorubicin","Idarubicin")

for (i in 1:length(cancer_ic50)) {
  for (j in 1:length(drug_seq_predict)) {
    
    o <- ic50_msi[ic50_msi$`Cancer.Type.(matching.TCGA.label)`%in%cancer_ic50[i],]
    p <- o[o$`Drug name`%in%drug_seq_predict[j],]
    length_msi[i,j] <- nrow(p)
  }
}
rownames(length_msi) <- cancer_ic50
colnames(length_msi) <- drug_seq_predict
for (i in 1:length(cancer_ic50)) {
  for (j in 1:length(drug_seq_predict)) {
    
    m <- ic50_mss[ic50_mss$`Cancer.Type.(matching.TCGA.label)`%in%cancer_ic50[i],]
    n <- m[m$`Drug name`%in%drug_seq_predict[j],]
    length_mss[i,j] <- nrow(n)
    
    rownames(length_mss)[i] <- cancer_ic50[i]
    colnames(length_mss)[j] <- drug_seq_predict[j]
  }
}
rownames(length_mss) <- cancer_ic50
colnames(length_mss) <- drug_seq_predict


unique(as.character(length_msi[3,]))

ic50_mss1 <- ic50_mss[!ic50_mss$PANEL_NAME%in%"Lymphoma",]
ic50_msi1 <- ic50_msi[!ic50_msi$PANEL_NAME%in%"Lymphoma",]
cancer_ic501 <- cancer_ic50[-3]
length_mss1 <- length_mss[-3,]
length_msi1 <- length_msi[-3,]
for (i in 1:nrow(length_mss1)) {
  length_mss1 <- length_mss1[,!length_mss1[i,]<=1]
}
for (i in 1:nrow(length_msi1)) {
  length_msi1 <- length_msi1[,!length_msi1[i,]<=1]
}
drug_seq_predict1 <- colnames(length_mss1)
length_msi1 <- length_msi1[,drug_seq_predict1]

length(which(unique(colnames(length_msi1))%in%unique(colnames(length_msi1))))==ncol(length_msi1)

ic50_mss2 <- ic50_mss1[ic50_mss1$`Drug name`%in%drug_seq_predict1,]
ic50_msi2 <- ic50_msi1[ic50_msi1$`Drug name`%in%drug_seq_predict1,]


# 4. u test and heatmap -------------------------------------------------------------

wiltest <- data.frame()
wiltest_pval <- data.frame()

for (i in 1:length(cancer_ic501)) {
  for (j in 1:length(drug_seq_predict1)) {
    
    k <- ic50_mss2[ic50_mss2$`Cancer.Type.(matching.TCGA.label)`%in%cancer_ic501[i],]
    m <- k[k$`Drug name`%in%drug_seq_predict1[j],]
    l <- ic50_msi2[ic50_msi2$`Cancer.Type.(matching.TCGA.label)`%in%cancer_ic501[i],]
    n <- l[l$`Drug name`%in%drug_seq_predict1[j],]
    
    o <- wilcox.test(n$AVERAGE,m$AVERAGE , exact = F)
    wiltest[i,j] <- o$statistic
    wiltest_pval[i,j] <- o$p.value
  }
}
wiltest <- t(wiltest)
colnames(wiltest) <- cancer_ic501
rownames(wiltest) <- drug_seq_predict1

wiltest_pval <- t(wiltest_pval)
colnames(wiltest_pval) <- cancer_ic501
rownames(wiltest_pval) <- drug_seq_predict1
# remove groups with NA
ic50_msi2[ic50_msi2$`Drug name`%in%"Procarbazine"&ic50_msi2$PANEL_NAME%in%"Prostate Cancer",]

wiltest1 <- wiltest[!is.na(wiltest_pval[,1]),]
wiltest_pval1 <- wiltest_pval[!is.na(wiltest_pval[,1]),]


range(ic50_mss2$AVERAGE)
range(ic50_msi2$AVERAGE)
fold_change <- data.frame()
for (i in 1:length(cancer_ic501)) {
  for (j in 1:length(drug_seq_predict1)) {
    r <- ic50_mss2[ic50_mss2$`Cancer.Type.(matching.TCGA.label)`%in%cancer_ic501[i],]
    s <- r[r$`Drug name`%in%drug_seq_predict1[j],]
    t <- ic50_msi2[ic50_msi2$`Cancer.Type.(matching.TCGA.label)`%in%cancer_ic501[i],]
    u <- t[t$`Drug name`%in%drug_seq_predict1[j],]
    
    mss_mean <- mean(s$AVERAGE)
    msi_mean <- mean(u$AVERAGE)
    fold_change[i,j] <- msi_mean/mss_mean
    
    rownames(fold_change)[i] <- cancer_ic501[i]
    colnames(fold_change)[j] <- drug_seq_predict1[j]
  }
}
# log2 transformation of fold changes
fold_change <- t(log2(fold_change))
fold_change1 <- fold_change[!is.na(wiltest_pval[,1]),]


p_value = matrix(
  ifelse(wiltest_pval1>0.05," ",ifelse(wiltest_pval1>0.01,"*",ifelse(wiltest_pval1>0.001,"**",ifelse(wiltest_pval1>0.0001,"***","****")))),
  nrow(wiltest_pval1)
)
rownames(p_value) <- rownames(wiltest_pval1)
colnames(p_value) <- colnames(wiltest_pval1)


target <- data.frame(Targets = c("Antimetabolite","Antimetabolite (DNA & RNA)",
                                 "DNA replication","DNA replication(Anthracycline)",
                                 "DNA replication(Antimetabolite)",
                                 "DNA replication(DNA crosslinker)","DNA replication(DNA crosslinker)",
                                 "DNA replication(dsDNA break induction)","DNA replication(Pyrimidine antimetabolite)",
                                 "DNA replication(TOP1)","DNA replication(TOP1)","DNA replication(TOP1)",
                                 "DNA replication(TOP2)","DNA replication(TOP2)",
                                 "Mitosis(Microtubule destabiliser)","Mitosis(Microtubule destabiliser)",
                                 "Mitosis(Microtubule destabiliser)","RNA polymerase",
                                 "Unclassified","Unclassified", "Unclassified","Unclassified",
                                 "Unclassified","Unclassified","Unclassified","Unclassified"))

row_ha <- rowAnnotation(df = target
                        , col = list(Targets = c(
                          "DNA replication" = "#82D7FF",
                          "DNA replication(dsDNA break induction)" = "#82FF99",
                          "DNA replication(DNA crosslinker)" = "#CCFF7A",
                          "DNA replication(Pyrimidine antimetabolite)" = "#DBFEDC",
                          "DNA replication(TOP1)" = "#FFDFFC",
                          "DNA replication(TOP2)" = "#F9E2AE",
                          "DNA replication(Anthracycline)" = "#4FCAD2",
                          "DNA replication(Antimetabolite)" = "#94E3D7",
                          "RNA polymerase" = "#FBC7BD",
                          "Unclassified" = "#EFB7DC",
                          "Antimetabolite (DNA & RNA)" = "#3FB4F8",
                          "Antimetabolite" = "#889EFF",
                          "Mitosis(Microtubule destabiliser)" = "#F7C184")))

legendtilte <- paste0("* p < 0.05","\n\n",
                      "**P < 0.01","\n\n",
                      "***P < 0.001","\n\n",
                      "****P < 0.0001","\n\n",
                      "logFC")


rownames(fold_change1)[2] <- "5-Fluorouracil"
rownames(fold_change1)[18] <- "Dactinomycin"

rownames(p_value)[2] <- "5-Fluorouracil"
rownames(p_value)[18] <- "Dactinomycin"

rownames(wiltest_pval1)[2] <- "5-Fluorouracil"
rownames(wiltest_pval1)[18] <- "Dactinomycin"

drug_seq_predict1[2] <- "5-Fluorouracil"
drug_seq_predict1[18] <- "Dactinomycin"

range(fold_change1)
# col_fun = circlize::colorRamp2(c(-2,0,2), c('#A9F1DF', "white",'#FFBBBB'))
col_fun = circlize::colorRamp2(c(-1,0,1), c('#2AB7CA', "white",'#E73D2A'))
# col_fun = circlize::colorRamp2(c(-2,0,2), c('#BDF1FF', "white",'#F4ABF2'))
# col_fun = circlize::colorRamp2(c(-2,0,2), c('#B5ECDB', "white",'#F0C6C6'))

Heatmap(as.matrix(fold_change1),
        name = legendtilte,
        column_title_gp = gpar(fontsize = 16),
        row_title_gp = gpar(fontsize = 16),
        # color
        col = col_fun,
        # grid color
        rect_gp = gpar(col = "#F4F5EF"),
        row_title = "Drugs", row_title_side = "right", column_title = "Cancers", column_title_side = "bottom",
        cluster_rows = FALSE,
        row_order = drug_seq_predict1,
        # column_split = 6,
        cluster_columns = TRUE,
        right_annotation = row_ha,
        cell_fun = function(j,i,x,y,w,h,col){
          grid.text(p_value[i,j],x,y)
        }
)


