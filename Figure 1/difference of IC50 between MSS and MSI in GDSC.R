# 1. data and R packages--------------------------------------------------------------------
rm(list = ls())
options(stringsAsFactors = FALSE)

library(pheatmap)
library(scales)
library(tidyr)
library(reshape2)
library(ggtree)
library(aplot)
library(ggheatmap)
library(ggh4x)
library(ComplexHeatmap)
library(tidyverse)
library(randomcoloR)
library(ggplot2)
library(openxlsx)
library(ggpubr)
library(ggsci)

# Data of ln IC50 containing various drugs for various tumors
ln_ic50 <- read.xlsx(
  "/mnt/DATA/ytj/GDSC/GDSC2_fitted_dose_response_25Feb20.xlsx")
# save Cell_Lines_Details as Cell_Lines_Details2
ms <- read.xlsx(
  "/mnt/DATA/ytj/GDSC/Cell_Lines_Details2.xlsx")

# Merge ms into ln_ic50
all_data <- merge(ln_ic50,ms,
                  by.x = "COSMIC_ID",by.y = "COSMIC.identifier",
                  all.x = TRUE)

gdsc_cell_line <- read.csv("/mnt/DATA/ytj/GDSC/model_list_20211124.csv",
                           header = TRUE)

# 2. data preprocessing ----------------------------------------------------

# List all drugs in the data
a <- unique(all_data$DRUG_NAME)

# List commonly used chemotherapy drugs
b <- unique(c("Cisplatin","Oxaliplatin","Carmustine","Doxorubicin","Etoposide","Gemcitabine","Irinotecan",
              # 5-FU is not usually classified as DDR related drugs, it mainly affects DNA synthesis
              "Methotrexate","Mitoxantrone","5-Fluorouracil",
              "Nelarabine","Topotecan","Epirubicin","Camptothecin","Cyclophosphamide",
              "Pemetrexed","Teniposide","Temozolomide","Podophyllotoxin bromide",
              "Cisplatin","Oxaliplatin" ,"Carmustine" ,"Vincristine","Vinblastine",
              "Cytarabine","Fludarabine","Nelarabine","Gemcitabine","Vinorelbine"
              ,"Vinorelbine","Docetaxel","Paclitaxel","Topotecan","Irinotecan"
              ,"Mitoxantrone","Cytarabine","4-Fluorouracil","Fludarabine","Nelarabine"
              ,"Epirubicin","Dactinomycin","Camptothecin","Cyclophosphamide","Topotecan"
              ,"Teniposide","Temozolomide","Podophyllotoxin bromide","Doxorubicin","Adriamycin",
              "Bleomycin","Dacarbazine","Cyclophosphamide","bleomycin","Vincristine","Oncovin",
              "Dexamethasone","Etoposide","Prednisone","Hydroxydaunorubicin","BCNU",
              "Carmustine","VP-16","Procarbazine","Chlorambucil","Clarithromycin",
              "Pomalidomide","Dexamethasone","Methotrexate","Oncovin","Prednisone",
              "Procarbazine","Docetaxel","Taxotere","Thalidomide","Dexamethasone",
              "Carboplatin","Daunorubicin","ara-C","Tioguanine",
              "Dexamethasone","Platinum Agent","platinum agent","Platinum agent",
              "Rituximab","Ifosfamide","Etoposide","Thalidomide","Epirubicin",
              "Capecitabine","Prednisone","Mitomycin","Fludarabine","Leucovorin",
              "Folinic Acid","Folinic acid","folinic acid","Idarubicin","G-CSF",
              "Busulfan","Melphalan","Leucovorin","Vinorelbine","Pegylated","Liposomal",
              "Carboplatin","Leucovorin","Ifosfamide","Actinomycin D","Mesna","Novantrone",
              "Mitomycin","Methotrexate","Mechlorethamine","Procarbazine","Prednisone",
              "6-Mercaptopurine","Purinethol","MOPP","Revlimid","Bendamustine",
              "Trastuzumab","Carboplatin","Pertuzumab","Perjeta","Thalidomide","Lomustine",
              "Dacarbazine","Thalidomide"))
drug=intersect(a,b)

# Screen out chemotherapy drug data from the total data
all_chem_data <- all_data[all_data$DRUG_NAME%in%drug,]

unique(all_chem_data$DRUG_NAME)
unique(gdsc_cell_line$tissue_status)
gdsc_cell_line1 <- gdsc_cell_line[gdsc_cell_line$COSMIC_ID%in%all_chem_data$COSMIC_ID,]

unique(gdsc_cell_line1$tissue_status)
gdsc_cell_line2 <- gdsc_cell_line1[gdsc_cell_line1$tissue_status%in%c("Metastasis","Tumour"),]
all_chem_data <- all_chem_data[all_chem_data$COSMIC_ID%in%gdsc_cell_line2$COSMIC_ID,]

# only keep key columns
var_data <- all_chem_data[,c("COSMIC_ID","LN_IC50",
                             "Microsatellite.instability.Status.(MSI)",
                             "Cancer.Type.(matching.TCGA.label)",
                             "DRUG_NAME","PUTATIVE_TARGET")]

# devided data into MMSS/MSI-L group and MSI-H group
mss <- var_data[
  var_data$`Microsatellite.instability.Status.(MSI)`%in%"MSS/MSI-L",]
msi <- var_data[
  var_data$`Microsatellite.instability.Status.(MSI)`%in%"MSI-H",]

# Intersection of tumors and drugs
cancer_mss <- unique(mss$`Cancer.Type.(matching.TCGA.label)`)
cancer_msi <- unique(msi$`Cancer.Type.(matching.TCGA.label)`)
cancer_pre <- cancer_mss[cancer_mss%in%cancer_msi]
cancer <- cancer_pre[-1]

mss_pre <- mss[                 
  mss$`Cancer.Type.(matching.TCGA.label)`%in%cancer,]
msi_pre <- msi[
  msi$`Cancer.Type.(matching.TCGA.label)`%in%cancer,]
unique(mss_pre$`Cancer.Type.(matching.TCGA.label)`)
unique(msi_pre$`Cancer.Type.(matching.TCGA.label)`)
unique(mss_pre$DRUG_NAME)
unique(msi_pre$DRUG_NAME)
drug




# 3. using u test for difference analysis ---------------------------------
# Remove the ln transformation of IC50 before taking log2 transformation for fold change
mss_inter_rm <- mss_pre
msi_inter_rm <- msi_pre

# !!Notice: i didn't change column name when transformation!!
mss_inter_rm$LN_IC50 <- exp(mss_pre$LN_IC50)
msi_inter_rm$LN_IC50 <- exp(msi_pre$LN_IC50)

length_mss_rm <- data.frame()
length_msi_rm <- data.frame()

for (i in 1:length(cancer)) {
  for (j in 1:length(drug)) {
    
    m <- mss_inter_rm[mss_inter_rm$`Cancer.Type.(matching.TCGA.label)`%in%cancer[i],]
    n <- m[m$DRUG_NAME%in%drug[j],]
    length_mss_rm[i,j] <- nrow(n)
    
    rownames(length_mss_rm)[i] <- cancer[i]
    colnames(length_mss_rm)[j] <- drug[j]
    
    o <- msi_inter_rm[msi_inter_rm$`Cancer.Type.(matching.TCGA.label)`%in%cancer[i],]
    p <- o[o$DRUG_NAME%in%drug[j],]
    length_msi_rm[i,j] <- nrow(p)
    
    rownames(length_msi_rm)[i] <- cancer[i]
    colnames(length_msi_rm)[j] <- drug[j]
  }
}
# Remove the combination of tumors and drugs with only one sample 
cancer <- c("OV","PRAD","COAD/READ","ALL","STAD","UCEC")
mss_inter_rm <- mss_inter_rm[mss_inter_rm$`Cancer.Type.(matching.TCGA.label)`%in%cancer,]
msi_inter_rm <- msi_inter_rm[msi_inter_rm$`Cancer.Type.(matching.TCGA.label)`%in%cancer,]

wiltest <- data.frame()
wiltest_pval <- data.frame()
# Error in `[<-.data.frame`(`*tmp*`, i, j, value = "0.117647058823529") : 
# for(i in 1:length())   not  for(i in length())
for (i in 1:length(cancer)) {
  for (j in 1:length(drug)) {
    
    k <- mss_inter_rm[mss_inter_rm$`Cancer.Type.(matching.TCGA.label)`%in%cancer[i],]
    m <- k[k$DRUG_NAME%in%drug[j],]
    l <- msi_inter_rm[msi_inter_rm$`Cancer.Type.(matching.TCGA.label)`%in%cancer[i],]
    n <- l[l$DRUG_NAME%in%drug[j],]
    
    o <- wilcox.test(n$LN_IC50,m$LN_IC50)
    wiltest[i,j] <- o$statistic
    wiltest_pval[i,j] <- o$p.value
  }
}
wiltest <- t(wiltest)
colnames(wiltest) <- cancer
rownames(wiltest) <- drug

wiltest_pval <- t(wiltest_pval)
colnames(wiltest_pval) <- cancer
rownames(wiltest_pval) <- drug

# fold change is needed for heatmap
fold_change <- data.frame()
for (i in 1:length(cancer)) {
  for (j in 1:length(drug)) {
    r <- mss_inter_rm[mss_inter_rm$`Cancer.Type.(matching.TCGA.label)`%in%cancer[i],]
    s <- r[r$DRUG_NAME%in%drug[j],]
    t <- msi_inter_rm[msi_inter_rm$`Cancer.Type.(matching.TCGA.label)`%in%cancer[i],]
    u <- t[t$DRUG_NAME%in%drug[j],]
    
    mss_mean <- mean(s$LN_IC50)
    msi_mean <- mean(u$LN_IC50)
    fold_change[i,j] <- msi_mean/mss_mean
    
    rownames(fold_change)[i] <- cancer[i]
    colnames(fold_change)[j] <- drug[j]
  }
}
# log2 transformation for fold change
fold_change <- log2(fold_change)


# 4. heatmap --------------------------------------------------------------

p_value = matrix(
  ifelse(wiltest_pval>0.05," ",ifelse(wiltest_pval>0.01,"*",ifelse(wiltest_pval>0.001,"**",ifelse(wiltest_pval>0.0001,"***","****")))),
  nrow(wiltest_pval)
)
rownames(p_value) <- rownames(wiltest_pval)
colnames(p_value) <- colnames(wiltest_pval)

# target of row annotations
target <- data.frame(Targets = c("Antimetabolite","Antimetabolite (DNA & RNA)",
                                 "DNA replication","DNA replication","DNA replication(Anthracycline)",
                                 "DNA replication(Antimetabolite)","DNA replication(DNA alkylating agent)",
                                 "DNA replication(DNA alkylating agent)","DNA replication(DNA alkylating agent)",
                                 "DNA replication(DNA crosslinker)","DNA replication(Pyrimidine antimetabolite)",
                                 "DNA replication(TOP1)","DNA replication(TOP1)","DNA replication(TOP1)",
                                 "DNA replication(TOP2)","Mitosis(Microtubule destabiliser)",
                                 "Mitosis(Microtubule destabiliser)","Mitosis(Microtubule destabiliser)",
                                 "Mitosis(Microtubule destabiliser)","RNA polymerase",
                                 "Unclassified","Unclassified","Unclassified"))


row_ha <- rowAnnotation(df = target
                        , col = list(Targets = c(
                          "DNA replication" = "#82D7FF",
                          "DNA replication(DNA alkylating agent)" = "#82FF99",
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


# change COAD/READ to CRC
rownames(fold_change)[3] <- "CRC"

drug_seq <- c("Cytarabine","5-Fluorouracil","Nelarabine","Teniposide","Epirubicin","Fludarabine","Temozolomide","Cyclophosphamide","Oxaliplatin","Cisplatin","Gemcitabine","Topotecan","Irinotecan","Camptothecin","Mitoxantrone","Docetaxel","Paclitaxel","Vinblastine","Vinorelbine","Dactinomycin","Podophyllotoxin bromide","Carmustine","Vincristine")

fold_change <- fold_change[,drug_seq]
p_value <- p_value[drug_seq,]
wiltest_pval <- wiltest_pval[drug_seq,]

# write.xlsx(wiltest_pval , file = "~/泛癌化疗药敏项目/result/supplementary_table/p_value_GDSC.xlsx")
# write.xlsx(t(as.matrix(fold_change)) , file = "~/泛癌化疗药敏项目/result/supplementary_table/fold_change_GDSC.xlsx")

# color and extreme value
range(fold_change)
col_fun = circlize::colorRamp2(c(-2,0,2), c('#2AB7CA', "white",'#E73D2A'))
# col_fun = circlize::colorRamp2(c(-2,0,2), c('#C6B5F5', "white",'#43AEAD'))
# col_fun = circlize::colorRamp2(c(-2,0,2), c('#BDF1FF', "white",'#F4ABF2'))
# col_fun = circlize::colorRamp2(c(-2,0,2), c('#B5ECDB', "white",'#F0C6C6'))

Heatmap(t(as.matrix(fold_change)),
        # legend title
        name = legendtilte,
        # coloumn and row title font
        column_title_gp = gpar(fontsize = 16),
        row_title_gp = gpar(fontsize = 16),
        # heatmap color
        col = col_fun,
        # grid color
        rect_gp = gpar(col = "#F4F5EF"),
        row_title = "Drugs", row_title_side = "right", column_title = "Cancers", column_title_side = "bottom",
        # columns cluster
        cluster_rows = FALSE,
        cluster_columns = TRUE,
        # column_split = 6,
        row_order = drug_seq,
        right_annotation = row_ha,
        cell_fun = function(j,i,x,y,w,h,col){
          grid.text(p_value[i,j],x,y)
        }
)


# 5. boxplot --------------------------------------------------------------

# ln IC50 is needed for boxplot as the plot with IC50 value is too large to print completely
crc_mss <- mss_pre[mss_pre$`Cancer.Type.(matching.TCGA.label)`%in%"COAD/READ",]
crc_msi <- msi_pre[msi_pre$`Cancer.Type.(matching.TCGA.label)`%in%"COAD/READ",]
crc <- rbind(crc_msi,crc_mss)

# drugs with significance in CRC
crc_pval <- data.frame(wiltest_pval[,3])
all_drug <- rownames(crc_pval)
pval_drug <- all_drug[which(crc_pval<0.05)]
pval_drug

# only data of drugs with significance in CRC
crc_drug <- crc[crc$DRUG_NAME%in%pval_drug,]
crc_drug1 <- crc_drug
colnames(crc_drug1)[5] <- "Drugs"
colnames(crc_drug1)[2] <- "ln IC50"
colnames(crc_drug1)[3] <- "MSI Status"

# plot
col_fun = circlize::colorRamp2(c(-2,0,2), c('#A9F1DF', "white",'#FFBBBB'))
# col_fun = circlize::colorRamp2(c(-2,0,2), c('#C6B5F5', "white",'#43AEAD'))
# col_fun = circlize::colorRamp2(c(-2,0,2), c('#BDF1FF', "white",'#F4ABF2'))
# col_fun = circlize::colorRamp2(c(-2,0,2), c('#B5ECDB', "white",'#F0C6C6'))

class(crc_drug1$`MSI Status`)
crc_drug1$`MSI Status` <- factor(crc_drug1$`MSI Status`)
class(crc_drug1$`MSI Status`)
# For the sake of aesthetics, it is necessary to make the direction of the box diagram the same as the heat map
# so it is necessary to manually adjust the x-axis of the box diagram

mypal <- pal_npg("nrc", alpha = 0.7)(9)
mypal
show_col(mypal)
# Lancet
mypal <- pal_lancet("lanonc", alpha = 0.7)(9)
mypal
show_col(mypal)
# new england
mypal <- pal_nejm("default", alpha = 0.7)(9)
mypal
show_col(mypal)
# jama
mypal <- pal_jama("default", alpha = 0.7)(9)
mypal
show_col(mypal)

ggboxplot(crc_drug1 , x = "Drugs" , y = "ln IC50" , 
          title = "CRC",
          color = "MSI Status" , 
          # palette = c('#A9F1DF', '#FFBBBB'),
          palette = c('#C6B5F5', '#43AEAD'),#'#1F78B4', '#33A02C'
          # palette = c('#BDF1FF', '#F4ABF2'),
          # palette = c('#B5ECDB', '#F0C6C6'),
          order = drug_seq,
          add = "jitter",
          x.text.angle = 90)+
  stat_compare_means(
    # comparisons = list(c("MSI-H", "MSS/MSI-L")),
    aes(group = `MSI Status`),
    label = "p.signif",
    # paired = T,
    method = "wilcox.test"
  )
