# 1. load R packages ------------------------------------------------------

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
# 2. read data ------------------------------------------------------------
# Sys.setenv("VROOM_CONNECTION_SIZE" = 99999999999999999999 * 99999999999999999999)

# This is the dataset after drug treatment, where 0M of data is considered pre medication data
nci <- getGEO("GSE116436",
              # destdir = "/mnt/DATA/ytj/NCI-60",
              AnnotGPL = F , getGPL = F)

ddr <- getGmt("/mnt/DATA/ytj/GSEA/ddr_final.gmt")

# 3. data preprocessing ---------------------------------------------------

# phenotype data
pdata <- pData(nci[[1]])

pdata_sep<- separate(data = pdata , col = 1 ,
                     into = c("cell lines" , "DRUG_NAME" , "dose" , "observe time") ,
                     sep = "_")

colnames(pdata_drug)[1] <- "Sample.Name"

unique(pdata_drug[,1])[!unique(pdata_drug[,1])%in%ms1$Sample.Name]
pdata_drug[,1] <- gsub("HS578T","Hs-578-T",pdata_drug[,1])
pdata_drug[,1] <- gsub("HS-578T","Hs-578-T",pdata_drug[,1])
pdata_drug[,1] <- gsub("IGROV1","IGROV-1",pdata_drug[,1])
pdata_drug[,1] <- gsub("IGR-OV1","IGROV-1",pdata_drug[,1])
# pdata_drug[,1] <- gsub("COLO205","COLO-205",pdata_drug[,1])
pdata_drug[,1] <- gsub("T-47D","T47D",pdata_drug[,1])
unique(pdata_drug[,1])[!unique(pdata_drug[,1])%in%ms1$Sample.Name]

pdata_drug <- merge(pdata_drug , ms1 , by.y = "Sample.Name")
unique(pdata_drug[,50])
# There are many null values in the data of large cell lung cancer and T-cell non Hodgkin's lymphoma, and neither tumor is the focus of our research. Simply remove them
pdata_drug <- pdata_drug[!pdata_drug[,50]%in%NA,]
# remove LUSC and LGG for samll sample quantity
pdata_drug <- pdata_drug[!pdata_drug$`Cancer.Type.(matching.TCGA.label)`%in%c("LUSC","LGG"),]

# before and after drug treatment
sample_0 <- pdata_drug[pdata_drug$dose%in%"0nM",]
sample_1 <- pdata_drug[!pdata_drug$dose%in%"0nM",]

pdata_drug <- pdata_sep[pdata_sep$DRUG_NAME%in%c("cisplatin","paclitaxel","topotecan",
                                                 "5-Azacytidine",
                                                 "doxorubicin",
                                                 "gemcitibine"),]


# expression data
exp <- exprs(nci[[1]])
norm_exp <- exp
library(hgu133a2.db)
ls("package:hgu133a2.db")
anno <- toTable(hgu133a2SYMBOL)
head(anno)


# transformation
norm_exp1 <- as.data.frame(norm_exp)
norm_exp1[,7624] <- rownames(norm_exp1)
colnames(norm_exp1)[7624] <- "probe_id"
exp_sym <- merge(norm_exp1 , anno , by.y = "probe_id")
# median
exp_sym$median <- apply(exp_sym[,2:7624] , 1 , median)
exp_sym <- exp_sym[order(exp_sym$symbol , exp_sym$median , decreasing = T),]
exp_sym1 <- exp_sym[!duplicated(exp_sym$symbol),]
rownames(exp_sym1) <- exp_sym1[,7625]
exp_sym1 <- exp_sym1[,-c(1,7625,7626)]

exp_0 <- exp_sym1[,sample_0[,5]]
exp_1 <- exp_sym1[,sample_1[,5]]

# Separate the data before and after the use of different drugs for each type of tumor，Put it in one list
unique(sample_0$`Cancer.Type.(matching.TCGA.label)`)
unique(sample_1$`Cancer.Type.(matching.TCGA.label)`)
n_cancer <- unique(sample_0$`Cancer.Type.(matching.TCGA.label)`)
n_drug <- c("cisplatin","paclitaxel","topotecan",
            "5-Azacytidine",
            "doxorubicin",
            "gemcitibine")

list_ <- list() 
list0 <- list()

for (i in 1:length(n_cancer)) {
  
  can0 <- sample_0[sample_0$`Cancer.Type.(matching.TCGA.label)`%in%n_cancer[i],]
  
  for (j in 1:length(n_drug)) {
    
    dru0 <- can0[can0$DRUG_NAME%in%n_drug[j],]
    
    list_[[j]] <- exp_0[,dru0$geo_accession]
    
  }
  
  list0[(i*6-5):(i*6)] <- list_
  
}

n_drug_ <- c("Cisplatin","Paclitaxel","Topotecan",
             "5-Azacytidine",
             "Doxorubicin",
             "Gemcitibine")

name_list0 <- data.frame()
for (i in 1:length(n_cancer)) {
  for (j in 1:length(n_drug)) {
    
    name_list0[j,i]  <-  paste(n_cancer[i] , n_drug_[j] , sep = "_")
  }
}
name_list0_ <- as.vector(as.matrix(name_list0))
names(list0) <- name_list0_

list1 <- list()

for (k in 1:length(n_cancer)) {
  
  can1 <- sample_1[sample_1$`Cancer.Type.(matching.TCGA.label)`%in%n_cancer[k],]
  
  for (l in 1:length(n_drug)) {
    
    dru1 <- can1[can1$DRUG_NAME%in%n_drug[l],]
    
    list_[[l]] <- exp_1[,dru1$geo_accession]
    
  }
  list1[(k*6-5):(k*6)] <- list_
}

name_list1 <- data.frame()
for (k in 1:length(n_cancer)) {
  for (l in 1:length(n_drug)) {
    name_list1[l,k]  <-  paste(n_cancer[k] , n_drug_[l] , sep = "_")
  }
}
name_list1_ <- as.vector(as.matrix(name_list1))
names(list1) <- name_list1_

# 4. GSVA -----------------------------------------------------------------

list0_matrix <- lapply(list0, as.matrix)
list0_num <- list()
for (i in 1:length(list0_matrix)) {
  list0_num[[i]] <- apply(list0_matrix[[i]], 2, as.numeric)
  rownames(list0_num[[i]]) <- rownames(list0_matrix[[i]])
}
names(list0_num) <- names(list0_matrix)

list1_matrix <- lapply(list1, as.matrix)
list1_num <- list()
for (i in 1:length(list1_matrix)) {
  list1_num[[i]] <- apply(list1_matrix[[i]], 2, as.numeric)
  rownames(list1_num[[i]]) <- rownames(list1_matrix[[i]])
}
names(list1_num) <- names(list1_matrix)

list0_gsva <- list()
for (i in 1:length(list0_num)) {
  list0_gsva[[i]] <- gsva(list0_num[[i]] , ddr , method = "ssgsea" , parallel.sz = 45)
}
names(list0_gsva) <- names(list0_num)

list1_gsva <- list()
for (i in 1:length(list1_num)) {
  list1_gsva[[i]] <- gsva(list1_num[[i]] , ddr , method = "ssgsea" , parallel.sz = 45)
}
names(list1_gsva) <- names(list1_num)



# 5. difference analysis --------------------------------------------------

p_pan <- matrix(nrow = 72 , ncol = 9)
for (i in 1:length(list0_gsva)) {
  for (j in 1:length(names(ddr))) {
    
    o <- wilcox.test(list1_gsva[[i]][j,] , list0_gsva[[i]][j,])
    p_pan[i,j] <- o$p.value
    
  }}
colnames(p_pan) <- names(ddr)
rownames(p_pan) <- names(list0_gsva)

fc_pan <- matrix(nrow = 72 , ncol = 9)
for (i in 1:length(list0_gsva)) {
  for (j in 1:length(names(ddr))) {
    
    list1_mean <- mean(list1_gsva[[i]][j,])
    list0_mean <- mean(list0_gsva[[i]][j,])
    fc_pan[i,j]<- list1_mean/list0_mean
  }
}
colnames(fc_pan) <- names(ddr)
rownames(fc_pan) <- names(list0_gsva)


# 6. heatmap --------------------------------------------------------------

legendtilte <- paste0("* p < 0.05","\n\n",
                      "**P < 0.01","\n\n",
                      "***P < 0.001","\n\n",
                      "****P < 0.0001","\n\n",
                      "FC")


colnames(fc_pan)[which(colnames(fc_pan)=="REACTOME_BASE_EXCISION_REPAIR")] = "BER"
colnames(fc_pan)[which(colnames(fc_pan)=="REACTOME_SUMOYLATION_OF_DNA_DAMAGE_RESPONSE_AND_REPAIR_PROTEINS")] = "DDR"
colnames(fc_pan)[which(colnames(fc_pan)=="REACTOME_NUCLEOTIDE_EXCISION_REPAIR")] = "NER"
colnames(fc_pan)[which(colnames(fc_pan)=="REACTOME_MISMATCH_REPAIR")] = "MMR"
colnames(fc_pan)[which(colnames(fc_pan)=="GOMF_SINGLE_STRANDED_DNA_BINDING")] = "SSB"
colnames(fc_pan)[which(colnames(fc_pan)=="REACTOME_NONHOMOLOGOUS_END_JOINING_NHEJ")] = "NHEJ"
colnames(fc_pan)[which(colnames(fc_pan)=="KEGG_HOMOLOGOUS_RECOMBINATION")] = "HR"
colnames(fc_pan)[which(colnames(fc_pan)=="REACTOME_FANCONI_ANEMIA_PATHWAY")] = "FA"
colnames(fc_pan)[which(colnames(fc_pan)=="REACTOME_DOUBLE_STRAND_BREAK_REPAIR")] = "DSB"
p_value = matrix(
  ifelse(p_pan>0.05," ",ifelse(p_pan>0.01,"*",ifelse(p_pan>0.001,"**",ifelse(p_pan>0.0001,"***","****")))),
  nrow(p_pan)
)
rownames(p_value) <- rownames(fc_pan)
colnames(p_value) <- colnames(fc_pan)

# col_fun = circlize::colorRamp2(c(0.85,1,1.02), c('#2fa1dd', "white",'#f87669'))

drug_seq <- c(
  "KIRC_5-Azacytidine",
  "LUAD_5-Azacytidine",
  "BRCA_5-Azacytidine",
  "ALL_5-Azacytidine", 
  "COAD_5-Azacytidine",
  "PRAD_5-Azacytidine", 
  "OV_5-Azacytidine",
  "LAML_5-Azacytidine",
  "LCML_5-Azacytidine",
  "MM_5-Azacytidine",
  "GBM_5-Azacytidine",
  "SKCM_5-Azacytidine",
  
  "KIRC_Paclitaxel",
  "LUAD_Paclitaxel",
  "BRCA_Paclitaxel",
  "ALL_Paclitaxel", 
  "COAD_Paclitaxel",
  "PRAD_Paclitaxel", 
  "OV_Paclitaxel",
  "LAML_Paclitaxel",
  "LCML_Paclitaxel",
  "MM_Paclitaxel",
  "GBM_Paclitaxel",
  "SKCM_Paclitaxel",
  
  "KIRC_Cisplatin",
  "LUAD_Cisplatin",
  "BRCA_Cisplatin",
  "ALL_Cisplatin", 
  "COAD_Cisplatin",
  "PRAD_Cisplatin", 
  "OV_Cisplatin",
  "LAML_Cisplatin",
  "LCML_Cisplatin",
  "MM_Cisplatin",
  "GBM_Cisplatin",
  "SKCM_Cisplatin",
  
  "KIRC_Gemcitibine",
  "LUAD_Gemcitibine",
  "BRCA_Gemcitibine",
  "ALL_Gemcitibine", 
  "COAD_Gemcitibine",
  "PRAD_Gemcitibine", 
  "OV_Gemcitibine",
  "LAML_Gemcitibine",
  "LCML_Gemcitibine",
  "MM_Gemcitibine",
  "GBM_Gemcitibine",
  "SKCM_Gemcitibine",
  
  
  "KIRC_Topotecan",
  "LUAD_Topotecan",
  "BRCA_Topotecan",
  "ALL_Topotecan", 
  "COAD_Topotecan",
  "PRAD_Topotecan", 
  "OV_Topotecan",
  "LAML_Topotecan",
  "LCML_Topotecan",
  "MM_Topotecan",
  "GBM_Topotecan",
  "SKCM_Topotecan",
  
  "KIRC_Doxorubicin",
  "LUAD_Doxorubicin",
  "BRCA_Doxorubicin",
  "ALL_Doxorubicin", 
  "COAD_Doxorubicin",
  "PRAD_Doxorubicin", 
  "OV_Doxorubicin",
  "LAML_Doxorubicin",
  "LCML_Doxorubicin",
  "MM_Doxorubicin",
  "GBM_Doxorubicin",
  "SKCM_Doxorubicin"
)
# annotation on the right side of the plot
annotation_drugseq1 <- separate(as.data.frame(drug_seq) , col = 1 , into = c("Cancer" , "Drug") , sep = "_")
rownames(annotation_drugseq1) <- drug_seq
# annotation_drugseq[,3] <- drug_seq
# annotation_drugseq[,1] <- factor(annotation_drugseq[,1])
# annotation_drugseq[,2] <- factor(annotation_drugseq[,2])
row_annotation1 <- rowAnnotation(df = annotation_drugseq1
                                 , col = list(Cancer= c("KIRC"= "#94E3D7",
                                                        "LUAD"="#DBFEDC" ,
                                                        "BRCA"="#3FB4F8" ,
                                                        "ALL"="#889EFF" ,
                                                        "COAD"="#82D7FF" ,
                                                        "PRAD"="#F7C184" ,
                                                        "OV"="#F9E2AE" ,
                                                        "LAML"="#82FF99" ,
                                                        "LCML"="#CCFF7A",
                                                        "MM"="#FFDFFC",# #EFB7DC
                                                        "GBM"="#4FCAD2",
                                                        "SKCM"="#FBC7BD"),
                                              Drug = c("5-Azacytidine" = "#3FB4F8" , "Paclitaxel" = "#EFB7DC" ,
                                                       "Gemcitibine" = "#FFDFFC" , "Cisplatin" = "#94E3D7" ,
                                                       "Doxorubicin" = "#DBFEDC" , "Topotecan" = "#F7C184"))
)

# col_fun = circlize::colorRamp2(c(0.95,1,1.05), c('#A9F1DF', "white",'#FFBBBB'))
col_fun = circlize::colorRamp2(c(-2,0,2), c('#2AB7CA', "white",'#E73D2A'))
# col_fun = circlize::colorRamp2(c(-2,0,2), c('#BDF1FF', "white",'#F4ABF2'))
# col_fun = circlize::colorRamp2(c(-2,0,2), c('#B5ECDB', "white",'#F0C6C6'))
p_value <- p_value[drug_seq,]
p_pan <- p_pan[drug_seq,]


Heatmap(fc_pan[drug_seq,], 
        
        name = legendtilte,
       
        column_title_gp = gpar(fontsize = 16),
        row_title_gp = gpar(fontsize = 16),
        
        col = col_fun,
        
        rect_gp = gpar(col = "#F4F5EF"),
        row_title = "Cancers_Drugs", row_title_side = "right", column_title = "Pathways", column_title_side = "bottom",
        cluster_rows = FALSE,
        # column_split = 6,
        cluster_columns = TRUE,
        right_annotation = row_annotation1,
        row_order = drug_seq,
        width = ncol(fc_pan)*unit(8, "mm"), 
        height = nrow(fc_pan)*unit(5, "mm"),
        cell_fun = function(j,i,x,y,w,h,col){
          grid.text(p_value[i,j],x,y)
        })

# 7. violin plot ----------------------------------------------------------
# extract COAD data，devide into msi and mss
sample0_coad_mss <- sample_0[sample_0$`Cancer.Type.(matching.TCGA.label)`%in%"COAD"&sample_0$`Microsatellite.instability.Status.(MSI)`%in%"MSS/MSI-L",]
sample0_coad_msi <- sample_0[sample_0$`Cancer.Type.(matching.TCGA.label)`%in%"COAD"&sample_0$`Microsatellite.instability.Status.(MSI)`%in%"MSI-H",]
sample1_coad_mss <- sample_1[sample_1$`Cancer.Type.(matching.TCGA.label)`%in%"COAD"&sample_1$`Microsatellite.instability.Status.(MSI)`%in%"MSS/MSI-L",]
sample1_coad_msi <- sample_1[sample_1$`Cancer.Type.(matching.TCGA.label)`%in%"COAD"&sample_1$`Microsatellite.instability.Status.(MSI)`%in%"MSI-H",]

gsva0_coad_overall <- list0_gsva[25:30]
gsva1_coad_overall <- list1_gsva[25:30]
names(gsva0_coad_overall) <- n_drug_
names(gsva1_coad_overall) <- n_drug_

gsva0_coad_mss <- gsva0_coad_overall
gsva0_coad_msi <- gsva0_coad_overall
for (i in 1:length(gsva0_coad_overall)) {
  gsva0_coad_mss[[i]] <- gsva0_coad_overall[[i]][,colnames(gsva0_coad_overall[[i]])%in%sample0_coad_mss$geo_accession]
  gsva0_coad_msi[[i]] <- gsva0_coad_overall[[i]][,colnames(gsva0_coad_overall[[i]])%in%sample0_coad_msi$geo_accession]
}

gsva1_coad_mss <- gsva1_coad_overall
gsva1_coad_msi <- gsva1_coad_overall
for (i in 1:length(gsva1_coad_overall)) {
  gsva1_coad_mss[[i]] <- gsva1_coad_overall[[i]][,colnames(gsva1_coad_overall[[i]])%in%sample1_coad_mss$geo_accession]
  gsva1_coad_msi[[i]] <- gsva1_coad_overall[[i]][,colnames(gsva1_coad_overall[[i]])%in%sample1_coad_msi$geo_accession]
}
# combine control and treatment data as a list
dose_overal <- list(gsva0_coad_overall , gsva1_coad_overall)
dose_overall <- list(dose_overal[[1]][[1]] , dose_overal[[1]][[2]] , dose_overal[[1]][[3]] , dose_overal[[1]][[4]] , dose_overal[[1]][[5]] , dose_overal[[1]][[6]] , 
                     dose_overal[[2]][[1]] , dose_overal[[2]][[2]] , dose_overal[[2]][[3]] , dose_overal[[2]][[4]] , dose_overal[[2]][[5]] , dose_overal[[2]][[6]])
dose_mi <- list(gsva0_coad_msi , gsva1_coad_msi)
dose_msi <- list(dose_mi[[1]][[1]] , dose_mi[[1]][[2]] , dose_mi[[1]][[3]] , dose_mi[[1]][[4]] , dose_mi[[1]][[5]] , dose_mi[[1]][[6]] , 
                 dose_mi[[2]][[1]] , dose_mi[[2]][[2]] , dose_mi[[2]][[3]] , dose_mi[[2]][[4]] , dose_mi[[2]][[5]] , dose_mi[[2]][[6]])
dose_ms <- list(gsva0_coad_mss , gsva1_coad_mss)
dose_mss <- list(dose_ms[[1]][[1]] , dose_ms[[1]][[2]] , dose_ms[[1]][[3]] , dose_ms[[1]][[4]] , dose_ms[[1]][[5]] , dose_ms[[1]][[6]] , 
                 dose_ms[[2]][[1]] , dose_ms[[2]][[2]] , dose_ms[[2]][[3]] , dose_ms[[2]][[4]] , dose_ms[[2]][[5]] , dose_ms[[2]][[6]])

# adjust the order to the one of heatmap
ddr_seq <- c("REACTOME_BASE_EXCISION_REPAIR","REACTOME_MISMATCH_REPAIR","REACTOME_FANCONI_ANEMIA_PATHWAY",
             "KEGG_HOMOLOGOUS_RECOMBINATION","REACTOME_NUCLEOTIDE_EXCISION_REPAIR","GOMF_SINGLE_STRANDED_DNA_BINDING",
             "REACTOME_SUMOYLATION_OF_DNA_DAMAGE_RESPONSE_AND_REPAIR_PROTEINS","REACTOME_DOUBLE_STRAND_BREAK_REPAIR",
             "REACTOME_NONHOMOLOGOUS_END_JOINING_NHEJ")
for (i in 1:length(dose_overall)) {
  dose_overall[[i]] <- dose_overall[[i]][ddr_seq,]
  dose_msi[[i]] <- dose_msi[[i]][ddr_seq,]
  dose_mss[[i]] <- dose_mss[[i]][ddr_seq,]
}

for (i in 1:6) {
  names(dose_overall)[i] <- paste0("Control" , sep = "_" , n_drug_[i])
  names(dose_overall)[6+i] <- paste0("Treatment" , sep = "_" , n_drug_[i])
  
  names(dose_mss)[i] <- paste0("Control" , sep = "_" , n_drug_[i])
  names(dose_mss)[6+i] <- paste0("Treatment" , sep = "_" , n_drug_[i])
  
  names(dose_msi)[i] <- paste0("Control" , sep = "_" , n_drug_[i])
  names(dose_msi)[6+i] <- paste0("Treatment" , sep = "_" , n_drug_[i])
  
}



mss_list <- dose_mss
msi_list <- dose_msi
overall_list <- dose_overall

# convert to long data
mss_list_long <- list()
for (i in 1:length(mss_list)) {
  mss_pre <- as.data.frame(mss_list[[i]])
  mss_pre[,(ncol(mss_pre)+1)] <- rownames(mss_pre)
  colnames(mss_pre)[ncol(mss_pre)] <- "Pathways"
  mss_list_long[[i]] <- melt(mss_pre,
                             id.var = c("Pathways"),
                             variable.names = "Samples",
                             value.name = "ssGSEA ES")
}
names(mss_list_long) <- names(mss_list)

msi_list_long <- list()
for (i in 1:length(msi_list)) {
  msi_pre <- as.data.frame(msi_list[[i]])
  msi_pre[,(ncol(msi_pre)+1)] <- rownames(msi_pre)
  colnames(msi_pre)[ncol(msi_pre)] <- "Pathways"
  msi_list_long[[i]] <- melt(msi_pre,
                             id.var = c("Pathways"),
                             variable.names = "Samples",
                             value.name = "ssGSEA ES")
}
names(msi_list_long) <- names(msi_list)

overall_list_long <- list()
for (i in 1:length(overall_list)) {
  overall_pre <- as.data.frame(overall_list[[i]])
  overall_pre[,(ncol(overall_pre)+1)] <- rownames(overall_pre)
  colnames(overall_pre)[ncol(overall_pre)] <- "Pathways"
  overall_list_long[[i]] <- melt(overall_pre,
                                 id.var = c("Pathways"),
                                 variable.names = "Samples",
                                 value.name = "ssGSEA ES")
}
names(overall_list_long) <- names(overall_list)


# abbreviation
msi_list_sep <- list()
for (i in 1:length(msi_list_long)) {
  
  msi_pre <- msi_list_long[[i]]
  
  msi_pre$Pathways[which(msi_pre$Pathways=="REACTOME_BASE_EXCISION_REPAIR")] = "BER"
  msi_pre$Pathways[which(msi_pre$Pathways=="REACTOME_SUMOYLATION_OF_DNA_DAMAGE_RESPONSE_AND_REPAIR_PROTEINS")] = "DDR"
  msi_pre$Pathways[which(msi_pre$Pathways=="REACTOME_NUCLEOTIDE_EXCISION_REPAIR")] = "NER"
  msi_pre$Pathways[which(msi_pre$Pathways=="REACTOME_MISMATCH_REPAIR")] = "MMR"
  msi_pre$Pathways[which(msi_pre$Pathways=="GOMF_SINGLE_STRANDED_DNA_BINDING")] = "SSB"
  msi_pre$Pathways[which(msi_pre$Pathways=="REACTOME_NONHOMOLOGOUS_END_JOINING_NHEJ")] = "NHEJ"
  msi_pre$Pathways[which(msi_pre$Pathways=="KEGG_HOMOLOGOUS_RECOMBINATION")] = "HR"
  msi_pre$Pathways[which(msi_pre$Pathways=="REACTOME_FANCONI_ANEMIA_PATHWAY")] = "FA"
  msi_pre$Pathways[which(msi_pre$Pathways=="REACTOME_DOUBLE_STRAND_BREAK_REPAIR")] = "DSB"
  
  msi_pre[,4] <- names(msi_list_long)[[i]]
  msi_list_sep[[i]] <- separate(data = msi_pre , col = 4 ,
                                into = c("Groups" , "Drugs") ,
                                sep = "_")
  # colnames(msi_pre)[4] <- "Drugs"
  # msi_list_sep[[i]] <- msi_pre
}
names(msi_list_sep) <- names(msi_list_long)


mss_list_sep <- list()
for (i in 1:length(mss_list_long)) {
  
  mss_pre <- mss_list_long[[i]]
  
  mss_pre$Pathways[which(mss_pre$Pathways=="REACTOME_BASE_EXCISION_REPAIR")] = "BER"
  mss_pre$Pathways[which(mss_pre$Pathways=="REACTOME_SUMOYLATION_OF_DNA_DAMAGE_RESPONSE_AND_REPAIR_PROTEINS")] = "DDR"
  mss_pre$Pathways[which(mss_pre$Pathways=="REACTOME_NUCLEOTIDE_EXCISION_REPAIR")] = "NER"
  mss_pre$Pathways[which(mss_pre$Pathways=="REACTOME_MISMATCH_REPAIR")] = "MMR"
  mss_pre$Pathways[which(mss_pre$Pathways=="GOMF_SINGLE_STRANDED_DNA_BINDING")] = "SSB"
  mss_pre$Pathways[which(mss_pre$Pathways=="REACTOME_NONHOMOLOGOUS_END_JOINING_NHEJ")] = "NHEJ"
  mss_pre$Pathways[which(mss_pre$Pathways=="KEGG_HOMOLOGOUS_RECOMBINATION")] = "HR"
  mss_pre$Pathways[which(mss_pre$Pathways=="REACTOME_FANCONI_ANEMIA_PATHWAY")] = "FA"
  mss_pre$Pathways[which(mss_pre$Pathways=="REACTOME_DOUBLE_STRAND_BREAK_REPAIR")] = "DSB"
  
  mss_pre[,4] <- names(mss_list_long)[[i]]
  mss_list_sep[[i]] <- separate(data = mss_pre , col = 4 ,
                                into = c("Groups" , "Drugs") ,
                                sep = "_")
  # colnames(mss_pre)[4] <- "Drugs"
  # mss_list_sep[[i]] <- mss_pre
}
names(mss_list_sep) <- names(mss_list_long)



overall_list_sep <- list()
for (i in 1:length(overall_list_long)) {
  
  overall_pre <- overall_list_long[[i]]
  
  overall_pre$Pathways[which(overall_pre$Pathways=="REACTOME_BASE_EXCISION_REPAIR")] = "BER"
  overall_pre$Pathways[which(overall_pre$Pathways=="REACTOME_SUMOYLATION_OF_DNA_DAMAGE_RESPONSE_AND_REPAIR_PROTEINS")] = "DDR"
  overall_pre$Pathways[which(overall_pre$Pathways=="REACTOME_NUCLEOTIDE_EXCISION_REPAIR")] = "NER"
  overall_pre$Pathways[which(overall_pre$Pathways=="REACTOME_MISMATCH_REPAIR")] = "MMR"
  overall_pre$Pathways[which(overall_pre$Pathways=="GOMF_SINGLE_STRANDED_DNA_BINDING")] = "SSB"
  overall_pre$Pathways[which(overall_pre$Pathways=="REACTOME_NONHOMOLOGOUS_END_JOINING_NHEJ")] = "NHEJ"
  overall_pre$Pathways[which(overall_pre$Pathways=="KEGG_HOMOLOGOUS_RECOMBINATION")] = "HR"
  overall_pre$Pathways[which(overall_pre$Pathways=="REACTOME_FANCONI_ANEMIA_PATHWAY")] = "FA"
  overall_pre$Pathways[which(overall_pre$Pathways=="REACTOME_DOUBLE_STRAND_BREAK_REPAIR")] = "DSB"
  
  overall_pre[,4] <- names(overall_list_long)[[i]]
  overall_list_sep[[i]] <- separate(data = overall_pre , col = 4 ,
                                    into = c("Groups" , "Drugs") ,
                                    sep = "_")
  # colnames(overall_pre)[4] <- "Drugs"
  # overall_list_sep[[i]] <- overall_pre
}
names(overall_list_sep) <- names(overall_list_long)



# combine every lists
msi_data <- rbind(msi_list_sep[[1]],msi_list_sep[[2]],msi_list_sep[[3]],
                  msi_list_sep[[4]],msi_list_sep[[5]],msi_list_sep[[6]],
                  msi_list_sep[[7]],msi_list_sep[[8]],msi_list_sep[[9]],
                  msi_list_sep[[10]],msi_list_sep[[11]],msi_list_sep[[12]])

mss_data <- rbind(mss_list_sep[[1]],mss_list_sep[[2]],mss_list_sep[[3]],
                  mss_list_sep[[4]],mss_list_sep[[5]],mss_list_sep[[6]],
                  mss_list_sep[[7]],mss_list_sep[[8]],mss_list_sep[[9]],
                  mss_list_sep[[10]],mss_list_sep[[11]],mss_list_sep[[12]])

overall_data <- rbind(overall_list_sep[[1]],overall_list_sep[[2]],overall_list_sep[[3]],
                      overall_list_sep[[4]],overall_list_sep[[5]],overall_list_sep[[6]],
                      overall_list_sep[[7]],overall_list_sep[[8]],overall_list_sep[[9]],
                      overall_list_sep[[10]],overall_list_sep[[11]],overall_list_sep[[12]])


ggviolin(overall_data, x = "Pathways", y = "ssGSEA ES", fill = "Groups" , color = "Groups" ,
         add = c("boxplot"), add.params = list(fill="white") ,
         x.text.angle = 90  , title = "Overall",
         palette = c('#F3DBB5', '#277ABB'))+
  stat_compare_means(aes(group = Groups), label = "p.signif" , hide.ns = TRUE , method = "wilcox.test")+
  theme(plot.title = element_text(size = 20))+
  facet_wrap(~Drugs , ncol = 3) +
  theme(legend.position = 'top')

ggviolin(msi_data, x = "Pathways", y = "ssGSEA ES", fill = "Groups" , color = "Groups" ,
         add = c("boxplot"), add.params = list(fill="white") ,
         x.text.angle = 90 , title = "MSI-H",
         palette = c('#F3DBB5', '#277ABB')) +
  theme(plot.title = element_text(size = 20))+
  stat_compare_means(aes(group = Groups), label = "p.signif" , hide.ns = TRUE , method = "wilcox.test")+
  facet_wrap(~Drugs , ncol = 3) + 
  theme(legend.position = 'top')

ggviolin(mss_data, x = "Pathways", y = "ssGSEA ES", fill = "Groups" , color = "Groups" ,
         add = c("boxplot"), add.params = list(fill="white") ,
         x.text.angle = 90 , title = "MSS/MSI-L",
         palette = c('#F3DBB5', '#277ABB'))+
  stat_compare_means(aes(group = Groups), label = "p.signif" , hide.ns = TRUE , method = "wilcox.test")+
  theme(plot.title = element_text(size = 20))+
  facet_wrap(~Drugs , ncol = 3) +
  theme(legend.position = 'top')








