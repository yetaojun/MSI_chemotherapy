# 1. load R packages -----------------------------------------------------------
library(ggsci)
# install.packages("xlsx")
library(xlsx)
library(scales)
library(ggthemes)
library(GSVA)
library(GSEABase)
library(limma)
library(tidyverse)
library(survival)
# devtools::install_github("NightingaleHealth/ggforestplot")
library(ggforestplot)
# install.packages("forestplot")
library(forestplot)
# BiocManager::install("ggforest")
library(ggforest)
library(survminer)
library(cowplot)
# BiocManager::install(c('sva' , 'preprocessCore', 'ridge')) #pRRophetic dependency package(5 followings)
library(car)
library(genefilter)
library(ridge)
library(sva)
library(preprocessCore)
library(pRRophetic)
library(org.Hs.eg.db)
library(clusterProfiler)
library(TCGAbiolinks)
library(limma)
library(stats)
library(ComplexHeatmap)
# install.packages("KMunicate")
library(KMunicate)
library(reshape2)
library(ggsci)


# 2. read data -----------------------------------------------------------------

sur <- read.table("/mnt/DATA/ytj/Xena/GDC-PANCAN.survival.tsv" ,
                  header = TRUE , sep = "\t")

phe <- read.table("/mnt/DATA/ytj/Xena/GDC-PANCAN.basic_phenotype.tsv" ,
                  header = TRUE , sep = "\t")

exp <- read.table("/mnt/DATA/ytj/Xena/GDC-PANCAN.htseq_fpkm-uq.tsv" ,
                  header = TRUE , sep = "\t")
mantis <- read.csv("/mnt/DATA/ytj/Xena/Mantis score.csv",
                   header = T)

# using ddr gene sets from MSigDB creat ddr_final.gmt
gsea <- getGmt("/mnt/DATA/ytj/GSEA/ddr_final.gmt")
gsea_old <- getGmt("/mnt/DATA/ytj/GSEA/ddr.gmt")
# 3. data preprocessing ----------------------------------------------------------------

# change rowname from number to ensembl
exp1 <- exp
rownames(exp1) <- exp1[,1]
exp1 <- exp1[,-1]

# preparation for the transformation of standardization from fpkmuq to tpm
exp_fpkmuq <- 2^exp1
exp_fpkmuq1 <- as.data.frame(exp_fpkmuq-1)

# gene symbol annotation file
a <- read.table('/mnt/DATA/ytj/TCGA/gencode.v22.annotation.gene.probeMap',
                header = T)
# a can cover all ensembl id from tcga
m <- unique(a$id)
n <- unique(rownames(exp_fpkmuq1))

# match
symboltrans <- a[match(rownames(exp_fpkmuq1),a$id),1:2]
symboltrans$id[-which(n%in%symboltrans$id)]

# some gene symbol match more than one ensenbl id
count_symbol <- table(symboltrans$gene)
# 2000+
length(which(symboltrans$gene %in% names(count_symbol[count_symbol>=2])))

# remove ensenbl id with low expression level in one gene symbol
colnames(symboltrans)=c('probe_id','symbol')  

# make sure that there is no null
symboltrans <- symboltrans[symboltrans$symbol != '',]

symboltrans <- symboltrans[symboltrans$probe_id%in%rownames(exp_fpkmuq1),]
expr_exp_fpkmuq1 <- exp_fpkmuq1[symboltrans$probe_id,]

# median
length(which(rownames(expr_exp_fpkmuq1)==symboltrans$probe_id))
symboltrans$median <- apply(expr_exp_fpkmuq1,1,median) 
# order by median
symboltrans <- symboltrans[order(symboltrans$symbol,symboltrans$median,decreasing = T),]
# remove duplicated symbol and keep ensembl gene with maximal expression value
symboltrans <- symboltrans[!duplicated(symboltrans$symbol),]

# change rownames from ensembl id to symbol
expr_exp_fpkmuq1 <- expr_exp_fpkmuq1[symboltrans$probe_id,] 
rownames(expr_exp_fpkmuq1)=symboltrans$symbol

# adjust colnames of phenotype data to colnames of expression data
col_exp_fpkmuq1 <- colnames(expr_exp_fpkmuq1)

# When the object to be gsub has a wildcard character such as'. ', 
# fixed=TRUE needs to be used to forcibly cancel the wildcard
col_exp_fpkmuq1 <- gsub(".","-",col_exp_fpkmuq1,fixed = TRUE)
exp_fpkmuq2 <- expr_exp_fpkmuq1
colnames(exp_fpkmuq2) <- col_exp_fpkmuq1

# change into tpm
fpkmuqToTpm_crc <- function(fpkmuq)
{
  exp(log(fpkmuq) - log(sum(fpkmuq)) + log(1e6))
}

exp_tpm <- apply(exp_fpkmuq2,2,fpkmuqToTpm_crc)

log_tpm <- log2(exp_tpm+1)

# Because survival data is not required for drug prediction, the data preprocessing method and cox regression analysis are different
phe1 <- phe[phe$sample_type_id<10,]
phe1[,8] <- phe1[1]
phe2 <- separate(data = phe1 , col = 1 , into = letters[1:4] , sep = "-")
phe2[,12] <- unite(phe2[,1:3] , col = "X_PATIENT" , sep = "-")
which(duplicated(phe2[,12]))
phe2[phe2$X_PATIENT%in%phe2[duplicated(phe2[,12]),12],c(4,12)]
colnames(log_tpm) 
colnames(phe2)[11] <- "sample"


unique(phe2[phe2$program%in%"TARGET",8])
phe3 <- phe2[phe2[,5]%in%"TCGA",]


unique(phe3$d)
# Samples with decimal points are not classified into tumor types, so delete invalid data of this type
phe4 <- phe3[!phe3$d%in%c("01A.1" , "02A.1" , "01A.2" , "01B.1"),]
phe5 <- separate(data = phe4 , col = 4 , into = letters[1:4] , sep = "")
phe5 <- phe5[,-4]
sur5 <- phe5[order(phe5$d),]
unique(phe5$d)

patient <- unique(phe5$X_PATIENT)
phe6 <- data.frame()
for (i in 1:length(patient)) {
  f <- phe5[phe5$X_PATIENT%in%patient[i],]
  if("A"%in%f$d){
    phe6[i,1:11] <- f[f$d%in%"A",]  }
  else{if("B"%in%f$d){
    phe6[i,1:11] <- f[f$d%in%"B",]}
    else{if("C"%in%f$d){
      phe6[i,1:11] <- f[f$d%in%"C",]}
      else{if("D"%in%f$d){
        phe6[i,1:11] <- f[f$d%in%"D",]}
        else{if("E"%in%f$d){
          phe6[i,1:11] <- f[f$d%in%"E",]}
          else{if("F"%in%f$d){
            phe6[i,1:11] <- f[f$d%in%"F",]}
            else{if("G"%in%f$d){
              phe6[i,1:11] <- f[f$d%in%"G",]}
              else{if("H"%in%f$d){
                phe6[i,1:11] <- f[f$d%in%"H",]}
                else{if("R"%in%f$d){
                  phe6[i,1:11] <- f[f$d%in%"R",]}
                  else{}}}}}}}}}
}
unique(phe6$d)
phe7 <- phe6[,-(4:6)]
phe8 <- phe7
phe8[,9] <- ""
phe8[,2:9] <- phe7
colnames(phe8)[2:9] <- colnames(phe7)
phe8[,1] <- phe8[,8]
colnames(phe8)[1] <- "sample"
phe8 <- phe8[,-8]
which(duplicated(phe8[,1]))

phe9 <- phe8
phe9$project_id <- gsub("TCGA-COAD" , "TCGA-CRC" , phe9$project_id)
phe9$project_id <- gsub("TCGA-READ" , "TCGA-CRC" , phe9$project_id)
unique(phe9$project_id)

inter_syl <- intersect(colnames(log_tpm) , phe9[,1])
phe_pre <- phe9[phe9[,1]%in%inter_syl,]
log_tpm_pre <- log_tpm[,inter_syl]

unique(phe_pre$project_id)
cancer <- unique(phe_pre$project_id)

cancer_syl_pre <- list()
for (i in 1:length(cancer)) {
  cancer_syl_pre[[i]] <- phe_pre[phe_pre$project_id%in%cancer[i],1]
}

# intersection of drugs
pRRdrug <- c("A.443654","A.770041","ABT.263","ABT.888","AG.014699","AICAR","AKT.inhibitor.VIII","AMG.706","AP.24534","AS601245","ATRA","AUY922","Axitinib","AZ628","AZD.0530","AZD.2281","AZD6244","AZD6482","AZD7762","AZD8055","BAY.61.3606","Bexarotene","BI.2536","BIBW2992","Bicalutamide","BI.D1870","BIRB.0796","Bleomycin","BMS.509744","BMS.536924","BMS.708163","BMS.754807","Bortezomib","Bosutinib","Bryostatin.1","BX.795","Camptothecin","CCT007093","CCT018159","CEP.701","CGP.082996","CGP.60474","CHIR.99021","CI.1040","Cisplatin","CMK","Cyclopamine","Cytarabine","Dasatinib","DMOG","Docetaxel","Doxorubicin","EHT.1864","Elesclomol","Embelin","Epothilone.B","Erlotinib","Etoposide","FH535","FTI.277","GDC.0449","GDC0941","Gefitinib","Gemcitabine","GNF.2","GSK269962A","GSK.650394","GW.441756","GW843682X","Imatinib","IPA.3","JNJ.26854165","JNK.9L","JNK.Inhibitor.VIII","JW.7.52.1","KIN001.135","KU.55933","Lapatinib","Lenalidomide","LFM.A13","Metformin","Methotrexate","MG.132","Midostaurin","Mitomycin.C","MK.2206","MS.275","Nilotinib","NSC.87877","NU.7441","Nutlin.3a","NVP.BEZ235","NVP.TAE684","Obatoclax.Mesylate","OSI.906","PAC.1","Paclitaxel","Parthenolide","Pazopanib","PD.0325901","PD.0332991","PD.173074","PF.02341066","PF.4708671","PF.562271","PHA.665752","PLX4720","Pyrimethamine","QS11","Rapamycin","RDEA119","RO.3306","Roscovitine","Salubrinal","SB.216763","SB590885","Shikonin","SL.0101.1","Sorafenib","S.Trityl.L.cysteine","Sunitinib","Temsirolimus","Thapsigargin","Tipifarnib","TW.37","Vinblastine","Vinorelbine","Vorinostat","VX.680","VX.702","WH.4.023","WO2009093972","WZ.1.84","X17.AAG","X681640","XMD8.85","Z.LLNle.CHO","ZM.447439")
data(cgp2016ExprRma)
data(PANCANCER_IC_Tue_Aug_9_15_28_57_2016)
data(bortezomibData)
possibleDrugs2016 <- unique(drugData2016$Drug.name)
b <- unique(c("Cisplatin","Oxaliplatin","Carmustine","Doxorubicin","Etoposide","Gemcitabine","Irinotecan",
              "Methotrexate","Mitoxantrone","5-Fluorouracil",# 5-FU does not affect DNA damage repair, but affects DNA synthesis
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
drug <- possibleDrugs2016[possibleDrugs2016%in%b]

tpm_list <- list()
for (i in 1:length(cancer)) {
  tpm_list[[i]] <- log_tpm_pre[,colnames(log_tpm_pre)%in%cancer_syl_pre[[i]]]
}
names(tpm_list) <- cancer

tpm_list1 <- tpm_list
for (i in 1:length(cancer)) {
  tpm_list1[[i]] <- as.data.frame(t(tpm_list1[[i]]))
  tpm_list1[[i]][,(ncol(tpm_list1[[i]])+1)] <- rownames(tpm_list1[[i]])
  tpm_list1[[i]] <- t(separate(tpm_list1[[i]] ,
                               col = ncol(tpm_list1[[i]]) , into = letters[1:4] , sep = "-"))
}

tpm_list2 <- tpm_list1
for (i in 1:length(cancer)) {
  tpm_list2[[i]] <- as.data.frame(t(tpm_list1[[i]][-nrow(tpm_list1[[i]]),]))
  tpm_list2[[i]] <- unite(tpm_list2[[i]] , col = "cancer" , sep = "-" , (ncol(tpm_list2[[i]])-2):ncol(tpm_list2[[i]]))
  rownames(tpm_list2[[i]]) <- tpm_list2[[i]][,ncol(tpm_list2[[i]])]
  tpm_list2[[i]] <- t(tpm_list2[[i]][,-ncol(tpm_list2[[i]])])
  tpm_list2[[i]] <- apply(as.matrix(tpm_list2[[i]]), 2, as.numeric)
  rownames(tpm_list2[[i]]) <- rownames(tpm_list[[i]])
}



# check how many tumors in the mantis score file
unique(mantis$Cancer_Type)
# combine COAD and READ
mantis1 <- mantis
mantis1$Cancer_Type <- gsub("TCGA-COAD" , "TCGA-CRC" , mantis1$Cancer_Type)
mantis1$Cancer_Type <- gsub("TCGA-READ" , "TCGA-CRC" , mantis1$Cancer_Type)
unique(mantis1$Cancer_Type)
cancer%in%unique(mantis1$Cancer_Type)
# Distinguishing mss and msi with 0.4 as the threshold
mss <- mantis1[mantis1$MANTIS_Score<0.4,]
msi <- mantis1[mantis1$MANTIS_Score>=0.4,]
unique(msi$Cancer_Type)
unique(mss$Cancer_Type)
cancer1 <- unique(msi$Cancer_Type)
mss <- mss[mss$Cancer_Type%in%cancer1,]
unique(msi$Cancer_Type)
unique(mss$Cancer_Type)

ms <- rbind(mss,msi)
exp_ms <- tpm_list2[cancer1]
exp_ms1 <- list()
for (i in 1:length(cancer1)) {
  exp_ms1[[i]] <- exp_ms[[i]][,colnames(exp_ms[[i]])%in%ms$ID]
}
names(exp_ms1) <- names(exp_ms)

# this prediction analysis will take 10 hours
predict_ms <- list()
for (i in 1:length(exp_ms1)) {
  predictedPtype <- data.frame(1:ncol(exp_ms1[[i]]))
  for (j in 1:length(drug)) {
    predictedPtype[,j] <- pRRopheticPredict(
      testMatrix=exp_ms1[[i]],
      drug=drug[j],
      tissueType = "all",
      batchCorrect = "eb",
      selection=1,
      dataset = "cgp2016")
  }
  colnames(predictedPtype) <- drug[j]
  rownames(predictedPtype) <- colnames(exp_ms1[[i]])
  predict_ms[[i]] <- predictedPtype
}
names(predict_ms) <- names(exp_ms1)

# add names to it
predict_ms1 <- list()
for (i in 1:length(cancer1)) {
  predict_ms1[[i]] <- t(predict_ms[[i]])
  rownames(predict_ms1[[i]]) <- drug
}
names(predict_ms1) <- names(predict_ms)

predict_ms2 <- list()
for (i in 1:length(cancer1)) {
  predict_ms2[[(2*i-1)]] <- predict_ms1[[i]][,colnames(predict_ms1[[i]])%in%mss[,1]]
  names(predict_ms2)[(2*i-1)] <- paste(names(predict_ms1)[i] , "MSS" , sep = "-")
  predict_ms2[[(2*i)]] <- predict_ms1[[i]][,colnames(predict_ms1[[i]])%in%msi[,1]]
  names(predict_ms2)[(2*i)] <- paste(names(predict_ms1)[i] , "MSI" , sep = "-")
}

# as.matrix
predict_ms3 <- predict_ms2
predict_ms3[[16]] <- as.matrix(predict_ms3[[16]])
predict_ms3[[42]] <- as.matrix(predict_ms3[[42]])

# remove the groups with null
predict_ms3 <- predict_ms3[-c(9,10)]

# 4. using u test for difference analysis --------------------------------------------------------------

cancer2 <- cancer1[-5]

predict_ms4 <-sapply(predict_ms3, function(x)2^x)

utest <- data.frame()
utest_p <- data.frame()
for(i in 1:length(cancer2)){
  for (j in 1:length(drug)) {
    
    x <- predict_ms4[[(2*i-1)]][j,]
    y <- predict_ms4[[(2*i)]][j,]
    
    wil <- wilcox.test(x,y)
    utest[j,i] <- wil$statistic
    utest_p[j,i] <- wil$p.value
    
  }
}
colnames(utest) <- cancer2
colnames(utest_p) <- cancer2
rownames(utest) <- drug
rownames(utest_p) <- drug


fc <- data.frame()
for (i in 1:length(cancer2)) {
  for (j in 1:length(drug)) {
    
    x <- predict_ms4[[(2*i-1)]][j,]
    y <- predict_ms4[[(2*i)]][j,]
    
    xz <- mean(x , na.rm = T)
    yz <- mean(y , na.rm = T)
    # note: msi/mss not mss/msi
    fc[j,i] <- yz/xz
  }
}

colnames(fc) <- cancer2
rownames(fc) <- drug
fc <- as.matrix(log2(fc))


col <- as.data.frame(colnames(fc))
sep_col <- separate(col , into = c("TCGA" , "Cancer") , col = 1 , sep = "-")
colnames(fc) <- sep_col$Cancer

# 5. heatmap-----------------------------------------------------------------
legendtilte <- paste0("* p < 0.05","\n\n",
                      "**P < 0.01","\n\n",
                      "***P < 0.001","\n\n",
                      "****P < 0.0001","\n\n",
                      "logFC")

p_value = matrix(
  ifelse(utest_p>0.05," ",ifelse(utest_p>0.01,"*",ifelse(utest_p>0.001,"**",ifelse(utest_p>0.0001,"***","****")))),
  nrow(utest_p)
)
rownames(p_value) <- rownames(utest_p)
colnames(p_value) <- colnames(utest_p)


# target of row annoation
target <- data.frame(Targets = c("Antimetabolite","DNA replication(Antimetabolite)",
                                 "DNA replication(DNA alkylating agent)","DNA replication(DNA crosslinker)",
                                 "DNA replication(dsDNA break induction)","DNA replication(Pyrimidine antimetabolite)",
                                 "DNA replication(TOP1)","DNA replication(TOP2)","Mitosis(Microtubule destabiliser)",
                                 "Mitosis(Microtubule destabiliser)","Mitosis(Microtubule destabiliser)",
                                 "Mitosis(Microtubule destabiliser)"))


row_ha <- rowAnnotation(df = target
                        , col = list(Targets = c(
                          "DNA replication(Antimetabolite)" = "#4FCAD2",
                          "DNA replication(DNA alkylating agent)" = "#94E3D7",
                          "DNA replication(dsDNA break induction)" = "#FFDFFC",
                          "DNA replication(DNA crosslinker)" = "#DBFEDC",
                          "DNA replication(Pyrimidine antimetabolite)" = "#F9E2AE",
                          "DNA replication(TOP1)" = "#F7C184",
                          "DNA replication(TOP2)" = "#FBC7BD",
                          "Antimetabolite" = "#3FB4F8",
                          "Mitosis(Microtubule destabiliser)" = "#EFB7DC")))

drug_seq <- c("Cytarabine","Methotrexate","Temozolomide","Cisplatin","Bleomycin","Gemcitabine","Camptothecin","Etoposide","Vinblastine","Vinorelbine","Docetaxel","Paclitaxel")
fc <- fc[drug_seq,]
p_value <- p_value[drug_seq,]
utest_p <- utest_p[drug_seq,]

range(fc)
col_fun = circlize::colorRamp2(c(-0.4,0,0.4), c('#2AB7CA', "white",'#E73D2A'))


Heatmap(fc,
        name = legendtilte,
        column_title_gp = gpar(fontsize = 16),
        row_title_gp = gpar(fontsize = 16),
        col = col_fun,
        # grid color
        rect_gp = gpar(col = "#F4F5EF"),
        row_title = "Drugs", row_title_side = "right", column_title = "Cancers", column_title_side = "bottom",
        # 6 clusters
        cluster_rows = FALSE,
        row_order = drug_seq,
        # column_split = 3,
        cluster_columns = TRUE,
        right_annotation = row_ha,
        width = ncol(fc)*unit(8, "mm"), 
        height = nrow(fc)*unit(5, "mm"),
        cell_fun = function(j,i,x,y,w,h,col){
          grid.text(p_value[i,j],x,y)
        })
# 4. boxplot --------------------------------------------------------


# changing from wide data to long data

predict_mss <- as.data.frame(predict_ms4[[9]])
predict_mss[,(ncol(predict_mss)+1)] <- rownames(predict_mss)
colnames(predict_mss)[ncol(predict_mss)] <- "drugs"

predict_mss <- melt(
  predict_mss,
  id.vars = c("drugs"), 
  variable.name = "ID",
  value.name = "ln IC50"
)

predict_msi <- as.data.frame(predict_ms4[[10]])
predict_msi[,(ncol(predict_msi)+1)] <- rownames(predict_msi)
colnames(predict_msi)[ncol(predict_msi)] <- "drugs"

predict_msi <- melt(
  predict_msi,
  id.vars = c("drugs"), 
  variable.name = "ID",
  value.name = "ln IC50"
)
# Replace mantis score with MSI-H and MSS/MSI-L
mss3 <- mss
msi3 <- msi
mss3[,3] <- "MSS/MSI-L"
msi3[,3] <- "MSI-H"

# remove TCGA from cancer name
mss3 <- separate(mss3 , col = 2 , into = c("TCGA" , "Cancer"), sep = "-")
msi3 <- separate(msi3 , col = 2 , into = c("TCGA" , "Cancer"), sep = "-")

# log e
predict_mss[,3] <- log(predict_mss[,3])
predict_msi[,3] <- log(predict_msi[,3])

# Integrating IC50 data with phenotype data from MSS and MSI
crc_mss <- merge(predict_mss , mss3 , by.y = "ID")
crc_msi <- merge(predict_msi , msi3 , by.y = "ID")

crc <- rbind(crc_msi,crc_mss)

# only keep data of drug with significance
crc_drug1 <- crc
colnames(crc_drug1)[2] <- "Drugs"
colnames(crc_drug1)[6] <- "MSI Status"

# plot
col_fun = circlize::colorRamp2(c(-2,0,2), c('#A9F1DF', "white",'#FFBBBB'))
# col_fun = circlize::colorRamp2(c(-2,0,2), c('#C6B5F5', "white",'#43AEAD'))
# col_fun = circlize::colorRamp2(c(-2,0,2), c('#BDF1FF', "white",'#F4ABF2'))
# col_fun = circlize::colorRamp2(c(-2,0,2), c('#B5ECDB', "white",'#F0C6C6'))

class(crc_drug1$`MSI Status`)
crc_drug1$`MSI Status` <- factor(crc_drug1$`MSI Status`)
class(crc_drug1$`MSI Status`)

crc_drug1 <- crc_drug1[crc_drug1$Drugs%in%drug_seq,]
ggboxplot(crc_drug1 , x = "Drugs" , y = "ln IC50" , 
          title = "CRC",
          color = "MSI Status" , 
          # palette = c('#A9F1DF', '#FFBBBB'),
          palette = c('#C6B5F5', '#43AEAD'),
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

