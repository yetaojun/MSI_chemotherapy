# 1. load R packages -----------------------------------------------------------------
library(ggplot2)
library(xlsx)
library(limma)
library(ggsci)
library(edgeR)
library(statmod)
library(GSA)
library(clusterProfiler)
library(dittoSeq)
library(escape)
library(GSVA)
library(enrichplot)
library(plyr)
library(ggrepel)
library(RColorBrewer)
library(gridExtra)
library(ggplot2)
library(GEOquery)
library(tidyr)
library(openxlsx)
library(GSEABase)
library(DOSE)
library(org.Hs.eg.db)
library(DESeq2)
library(reshape2)
library(ggthemes)
library(pRRophetic)
library(IOBR) 

# 2. read data ------------------------------------------------------------

ddr <- getGmt("/mnt/DATA/ytj/GSEA/ddr_final.gmt")


ms <- read.xlsx(
  "/mnt/DATA/ytj/GDSC/Cell_Lines_Details2.xlsx")

gdsc_expr <- read.csv("/mnt/DATA/ytj/GDSC/rnaseq_tpm_20220118.csv",
                      header = TRUE
) 

# contains information on normal/tumor cell lines
gdsc_cell_line <- read.csv("/mnt/DATA/ytj/GDSC/model_list_20211124.csv",
                           header = TRUE)

gdsc_cell_line1 <- gdsc_cell_line[,c(17 , 18 , 39)]
colnames(gdsc_cell_line1)[3] <- "COSMIC.identifier"
ms <- merge(gdsc_cell_line1 , ms , by.y = "COSMIC.identifier")
# only keep metastasis and tumor samples
ms <- ms[ms$tissue_status%in%c("Metastasis","Tumour"),]

# data of ln IC50 containing various drugs for various tumors
ln_ic50 <- read.xlsx(
  "/mnt/DATA/ytj/GDSC/GDSC2_fitted_dose_response_25Feb20.xlsx")


# 3. data preprocessing ---------------------------------------------------

# adjust the rownames of gdsc_expr
gdsc_spl <- colnames(gdsc_expr)
gdsc_spl <- gsub(".","-",gdsc_spl,fixed = TRUE)
gdsc_expr1 <- gdsc_expr
colnames(gdsc_expr1) <- gdsc_spl

# gene SEPTIN4 duplicating
length(which(duplicated(gdsc_expr1[,1])))
# remove the one with low expression level
gdsc_expr1 <- gdsc_expr1[-2891,]    
genename_expr1 <- gdsc_expr1[,1]
gdsc_expr1 <- apply(gdsc_expr1[,-1],2,as.numeric) 

rownames(gdsc_expr1) <- genename_expr1 


# extract CRC data 
coad_ms_crc <- ms[ms$`Cancer.Type.(matching.TCGA.label)`%in%"COAD/READ",]

# devide into msi and mss
coad_mss_crc <- coad_ms_crc[coad_ms_crc$`Microsatellite.instability.Status.(MSI)`%in%"MSS/MSI-L",]
coad_msi_crc <- coad_ms_crc[coad_ms_crc$`Microsatellite.instability.Status.(MSI)`%in%"MSI-H",]

coad_mss_expr_crc <- gdsc_expr1[,colnames(gdsc_expr1)%in%coad_mss_crc$Sample.Name]
coad_msi_expr_crc <- gdsc_expr1[,colnames(gdsc_expr1)%in%coad_msi_crc$Sample.Name]

coad_mss_count_crc <- gdsc_count1[,colnames(gdsc_count1)%in%coad_mss_crc$Sample.Name]
coad_msi_count_crc <- gdsc_count1[,colnames(gdsc_count1)%in%coad_msi_crc$Sample.Name]

coad_expr_crc <- as.matrix(cbind(coad_mss_expr_crc, coad_msi_expr_crc))
coad_expr_crc <- coad_expr_crc[rownames(coad_expr_crc)[order(rownames(coad_expr_crc), decreasing = T)],]


# 4. ssGSEA and chemotherapy prediction -----------------------------------

coad_mss_matrix_crc <- as.matrix(coad_mss_expr_crc)
gdsc_mss_gsva_crc_ddr <- gsva(coad_mss_matrix_crc, ddr ,kcdf = "Gaussian" , method = "ssgsea", parallel.sz = 75)

coad_msi_matrix_crc <- as.matrix(coad_msi_expr_crc)
gdsc_msi_gsva_crc_ddr <- gsva(coad_msi_matrix_crc, ddr ,kcdf = "Gaussian" , method = "ssgsea", parallel.sz = 75)

expr_line_crc <- cbind(coad_msi_matrix_crc,coad_mss_matrix_crc)
gsva_line_crc <- gsva(expr_line_crc, ddr ,kcdf = "Gaussian" , method = "ssgsea", parallel.sz = 75)

pRRdrug <- c("A.443654","A.770041","ABT.263","ABT.888","AG.014699","AICAR","AKT.inhibitor.VIII","AMG.706","AP.24534","AS601245","ATRA","AUY922","Axitinib","AZ628","AZD.0530","AZD.2281","AZD6244","AZD6482","AZD7762","AZD8055","BAY.61.3606","Bexarotene","BI.2536","BIBW2992","Bicalutamide","BI.D1870","BIRB.0796","Bleomycin","BMS.509744","BMS.536924","BMS.708163","BMS.754807","Bortezomib","Bosutinib","Bryostatin.1","BX.795","Camptothecin","CCT007093","CCT018159","CEP.701","CGP.082996","CGP.60474","CHIR.99021","CI.1040","Cisplatin","CMK","Cyclopamine","Cytarabine","Dasatinib","DMOG","Docetaxel","Doxorubicin","EHT.1864","Elesclomol","Embelin","Epothilone.B","Erlotinib","Etoposide","FH535","FTI.277","GDC.0449","GDC0941","Gefitinib","Gemcitabine","GNF.2","GSK269962A","GSK.650394","GW.441756","GW843682X","Imatinib","IPA.3","JNJ.26854165","JNK.9L","JNK.Inhibitor.VIII","JW.7.52.1","KIN001.135","KU.55933","Lapatinib","Lenalidomide","LFM.A13","Metformin","Methotrexate","MG.132","Midostaurin","Mitomycin.C","MK.2206","MS.275","Nilotinib","NSC.87877","NU.7441","Nutlin.3a","NVP.BEZ235","NVP.TAE684","Obatoclax.Mesylate","OSI.906","PAC.1","Paclitaxel","Parthenolide","Pazopanib","PD.0325901","PD.0332991","PD.173074","PF.02341066","PF.4708671","PF.562271","PHA.665752","PLX4720","Pyrimethamine","QS11","Rapamycin","RDEA119","RO.3306","Roscovitine","Salubrinal","SB.216763","SB590885","Shikonin","SL.0101.1","Sorafenib","S.Trityl.L.cysteine","Sunitinib","Temsirolimus","Thapsigargin","Tipifarnib","TW.37","Vinblastine","Vinorelbine","Vorinostat","VX.680","VX.702","WH.4.023","WO2009093972","WZ.1.84","X17.AAG","X681640","XMD8.85","Z.LLNle.CHO","ZM.447439")
data(cgp2016ExprRma)
data(PANCANCER_IC_Tue_Aug_9_15_28_57_2016)
data(bortezomibData)
possibleDrugs2016 <- unique(drugData2016$Drug.name)
b <- unique(c("Cisplatin","Oxaliplatin","Carmustine","Doxorubicin","Etoposide","Gemcitabine","Irinotecan",
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
drug <- possibleDrugs2016[possibleDrugs2016%in%b]

predicted_gdsc <- data.frame(1:ncol(coad_expr_crc))
coad_expr_crc1 <- coad_expr_crc
colnames(coad_expr_crc1) <- paste0("s__",colnames(coad_expr_crc1))  
coad_expr_crc1 <- as.matrix(coad_expr_crc1) 
for (j in 1:length(drug)) {
  predicted_gdsc[,j] <- pRRopheticPredict(
    testMatrix = coad_expr_crc1, 
    drug=drug[j],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016")
}

# names
predicted_gdsc1 <- 2^predicted_gdsc
colnames(predicted_gdsc1) <- drug
rownames(predicted_gdsc1) <- colnames(coad_expr_crc)


ddr_drug <- c("Cisplatin","Doxorubicin","Etoposide","Gemcitabine","Irinotecan","Methotrexate",
              "Nelarabine","Oxaliplatin","Pemetrexed","Teniposide","Temozolomide","Cyclophosphamide",
              "Epirubicin","Fludarabine","Mitoxantrone","Topotecan","Camptothecin")

ic50_ddr <- ln_ic50[ln_ic50$DRUG_NAME%in%ddr_drug,]

# Transforming gsva data from wide data into long data
# overall
gsva_line1_crc <- as.data.frame(gsva_line_crc)
gsva_line1_crc[,(ncol(gsva_line_crc)+1)] <- rownames(gsva_line1_crc)
colnames(gsva_line1_crc)[(ncol(gsva_line_crc)+1)] <- "pathways"

gsva_line1_long_crc <- melt(
  gsva_line1_crc,
  id.vars = c("pathways"),
  variable.name = "cell lines",
  value.name = "ssgsea score"
)
colnames(ic50_ddr)[5] <- "cell lines"
line_crc <- merge(gsva_line1_long_crc,ic50_ddr,intersect(colnames(gsva_line1_long_crc),
                                                         colnames(ic50_ddr)))

line1_crc <- line_crc

line1_crc$pathways[which(line1_crc$pathways=="REACTOME_BASE_EXCISION_REPAIR")] = "BER"
line1_crc$pathways[which(line1_crc$pathways=="REACTOME_SUMOYLATION_OF_DNA_DAMAGE_RESPONSE_AND_REPAIR_PROTEINS")] = "DDR"
line1_crc$pathways[which(line1_crc$pathways=="REACTOME_NUCLEOTIDE_EXCISION_REPAIR")] = "NER"
line1_crc$pathways[which(line1_crc$pathways=="REACTOME_MISMATCH_REPAIR")] = "MMR"
line1_crc$pathways[which(line1_crc$pathways=="GOMF_SINGLE_STRANDED_DNA_BINDING")] = "SSB"
line1_crc$pathways[which(line1_crc$pathways=="REACTOME_NONHOMOLOGOUS_END_JOINING_NHEJ")] = "NHEJ"
line1_crc$pathways[which(line1_crc$pathways=="KEGG_HOMOLOGOUS_RECOMBINATION")] = "HR"
line1_crc$pathways[which(line1_crc$pathways=="REACTOME_FANCONI_ANEMIA_PATHWAY")] = "FA"
line1_crc$pathways[which(line1_crc$pathways=="REACTOME_DOUBLE_STRAND_BREAK_REPAIR")] = "DSB"

line1_crc[,22] <- "#E9E2EC"
colnames(line1_crc)[22] <- "Color"

line2_crc <- line1_crc
line2_crc[line2_crc$pathways%in%"BER"&line2_crc$DRUG_NAME%in%c("Camptothecin","Mitoxantrone","Irinotecan","Oxaliplatin","Topotecan"),]$Color <- "#63BCC9"#"#A7D676"
line2_crc[line2_crc$pathways%in%"DSB"&line2_crc$DRUG_NAME%in%c("Topotecan"),]$Color <- "#63BCC9"#"#63BCC9"
line2_crc[line2_crc$pathways%in%"MMR"&line2_crc$DRUG_NAME%in%c("Camptothecin","Irinotecan","Oxaliplatin","Epirubicin","Mitoxantrone","Oxaliplatin","Topotecan"),]$Color <- "#63BCC9"#"#F9E2AE"
line2_crc[line2_crc$pathways%in%"NER"&line2_crc$DRUG_NAME%in%c("Topotecan"),]$Color <- "#63BCC9"#"#FBC78D"
line2_crc[line2_crc$pathways%in%"NHEJ"&line2_crc$DRUG_NAME%in%c("Camptothecin","Cisplatin","Irinotecan","Oxaliplatin","Teniposide","Mitoxantrone","Oxaliplatin","Topotecan"),]$Color <- "#63BCC9"#"#ECB7CB"
line2_crc[line2_crc$pathways%in%"DDR"&line2_crc$DRUG_NAME%in%c("Camptothecin","Irinotecan","Oxaliplatin","Topotecan"),]$Color <- "#63BCC9"#"#85CBCC"
line3_crc <- line2_crc[!line2_crc$DRUG_NAME%in%c("Nelarabine","Temozolomide","Cyclophosphamide","Fludarabine","Gemcitabine"),]


# msi
gdsc_msi_gsva_crc_ddr1 <- as.data.frame(gdsc_msi_gsva_crc_ddr)
gdsc_msi_gsva_crc_ddr1[,(ncol(gdsc_msi_gsva_crc_ddr)+1)] <- rownames(gdsc_msi_gsva_crc_ddr1)
colnames(gdsc_msi_gsva_crc_ddr1)[(ncol(gdsc_msi_gsva_crc_ddr)+1)] <- "pathways"

gsva_msi1_long_crc <- melt(
  gdsc_msi_gsva_crc_ddr1,
  id.vars = c("pathways"),
  variable.name = "cell lines",
  value.name = "ssgsea score"
)
msi_crc <- merge(gsva_msi1_long_crc,ic50_ddr,intersect(colnames(gsva_msi1_long_crc),
                                                       colnames(ic50_ddr)))

msi1_crc <- msi_crc

msi1_crc$pathways[which(msi1_crc$pathways=="REACTOME_BASE_EXCISION_REPAIR")] = "BER"
msi1_crc$pathways[which(msi1_crc$pathways=="REACTOME_SUMOYLATION_OF_DNA_DAMAGE_RESPONSE_AND_REPAIR_PROTEINS")] = "DDR"
msi1_crc$pathways[which(msi1_crc$pathways=="REACTOME_NUCLEOTIDE_EXCISION_REPAIR")] = "NER"
msi1_crc$pathways[which(msi1_crc$pathways=="REACTOME_MISMATCH_REPAIR")] = "MMR"
msi1_crc$pathways[which(msi1_crc$pathways=="GOMF_SINGLE_STRANDED_DNA_BINDING")] = "SSB"
msi1_crc$pathways[which(msi1_crc$pathways=="REACTOME_NONHOMOLOGOUS_END_JOINING_NHEJ")] = "NHEJ"
msi1_crc$pathways[which(msi1_crc$pathways=="KEGG_HOMOLOGOUS_RECOMBINATION")] = "HR"
msi1_crc$pathways[which(msi1_crc$pathways=="REACTOME_FANCONI_ANEMIA_PATHWAY")] = "FA"
msi1_crc$pathways[which(msi1_crc$pathways=="REACTOME_DOUBLE_STRAND_BREAK_REPAIR")] = "DSB"

msi1_crc[,22] <- "#E9E2EC"
colnames(msi1_crc)[22] <- "Color"

msi2_crc <- msi1_crc
msi2_crc[msi2_crc$pathways%in%"BER"&msi2_crc$DRUG_NAME%in%c("Irinotecan","Topotecan"),]$Color <- "#63BCC9"#"#A7D676"
msi2_crc[msi2_crc$pathways%in%"NER"&msi2_crc$DRUG_NAME%in%c("Topotecan"),]$Color <- "#63BCC9"#"#FBC78D"
msi2_crc[msi2_crc$pathways%in%"NHEJ"&msi2_crc$DRUG_NAME%in%c("Camptothecin","Irinotecan","Oxaliplatin","Topotecan"),]$Color <- "#63BCC9"#"#ECB7CB"
msi2_crc[msi2_crc$pathways%in%"DDR"&msi2_crc$DRUG_NAME%in%c("Irinotecan","Topotecan"),]$Color <- "#63BCC9"#"#85CBCC"
msi3_crc <- msi2_crc[!msi2_crc$DRUG_NAME%in%c("Nelarabine","Temozolomide","Cyclophosphamide","Fludarabine","Gemcitabine"),]

# mss
gdsc_mss_gsva_crc_ddr1 <- as.data.frame(gdsc_mss_gsva_crc_ddr)
gdsc_mss_gsva_crc_ddr1[,(ncol(gdsc_mss_gsva_crc_ddr)+1)] <- rownames(gdsc_mss_gsva_crc_ddr1)
colnames(gdsc_mss_gsva_crc_ddr1)[(ncol(gdsc_mss_gsva_crc_ddr)+1)] <- "pathways"

gsva_mss1_long_crc <- melt(
  gdsc_mss_gsva_crc_ddr1,
  id.vars = c("pathways"),
  variable.name = "cell lines",
  value.name = "ssgsea score"
)
mss_crc <- merge(gsva_mss1_long_crc,ic50_ddr,intersect(colnames(gsva_mss1_long_crc),
                                                       colnames(ic50_ddr)))

mss1_crc <- mss_crc

mss1_crc$pathways[which(mss1_crc$pathways=="REACTOME_BASE_EXCISION_REPAIR")] = "BER"
mss1_crc$pathways[which(mss1_crc$pathways=="REACTOME_SUMOYLATION_OF_DNA_DAMAGE_RESPONSE_AND_REPAIR_PROTEINS")] = "DDR"
mss1_crc$pathways[which(mss1_crc$pathways=="REACTOME_NUCLEOTIDE_EXCISION_REPAIR")] = "NER"
mss1_crc$pathways[which(mss1_crc$pathways=="REACTOME_MISMATCH_REPAIR")] = "MMR"
mss1_crc$pathways[which(mss1_crc$pathways=="GOMF_SINGLE_STRANDED_DNA_BINDING")] = "SSB"
mss1_crc$pathways[which(mss1_crc$pathways=="REACTOME_NONHOMOLOGOUS_END_JOINING_NHEJ")] = "NHEJ"
mss1_crc$pathways[which(mss1_crc$pathways=="KEGG_HOMOLOGOUS_RECOMBINATION")] = "HR"
mss1_crc$pathways[which(mss1_crc$pathways=="REACTOME_FANCONI_ANEMIA_PATHWAY")] = "FA"
mss1_crc$pathways[which(mss1_crc$pathways=="REACTOME_DOUBLE_STRAND_BREAK_REPAIR")] = "DSB"

mss1_crc[,22] <- "#E9E2EC"
colnames(mss1_crc)[22] <- "Color"

mss2_crc <- mss1_crc
mss2_crc[mss2_crc$pathways%in%"BER"&mss2_crc$DRUG_NAME%in%c("Camptothecin","Irinotecan","Oxaliplatin"),]$Color <- "#63BCC9"#"#A7D676"
mss2_crc[mss2_crc$pathways%in%"DSB"&mss2_crc$DRUG_NAME%in%c("Camptothecin","Oxaliplatin","Topotecan"),]$Color <- "#63BCC9"#"#63BCC9"
mss2_crc[mss2_crc$pathways%in%"MMR"&mss2_crc$DRUG_NAME%in%c("Camptothecin","Oxaliplatin","Topotecan"),]$Color <- "#63BCC9"#"#F9E2AE"
mss2_crc[mss2_crc$pathways%in%"NER"&mss2_crc$DRUG_NAME%in%c("Oxaliplatin","Topotecan"),]$Color <- "#63BCC9"#"#FBC78D"
mss2_crc[mss2_crc$pathways%in%"SSB"&mss2_crc$DRUG_NAME%in%c("Topotecan"),]$Color <- "#63BCC9"#"#D78FB7"
mss2_crc[mss2_crc$pathways%in%"DDR"&mss2_crc$DRUG_NAME%in%c("Camptothecin","Oxaliplatin","Topotecan"),]$Color <- "#63BCC9"#"#85CBCC"
mss3_crc <- mss2_crc[!mss2_crc$DRUG_NAME%in%c("Nelarabine","Temozolomide","Cyclophosphamide","Fludarabine","Gemcitabine"),]


pathway <- unique(line3_crc$pathways)
drug_final <- unique(line3_crc$DRUG_NAME)

# 5. overall heatmap --------------------------------------------------------------



heatmap_data <- merge(pathway , drug_final)
colnames(heatmap_data)[1] <- "pathways"
colnames(heatmap_data)[2] <- "drugs"

for (i in 1:length(pathway)) {
  for (j in 1:length(drug_final)) {
    
    ssgsea <- line3_crc[line3_crc$pathways%in%pathway[i]&line3_crc$DRUG_NAME%in%drug_final[j],3]
    ic50 <- line3_crc[line3_crc$pathways%in%pathway[i]&line3_crc$DRUG_NAME%in%drug_final[j],18]
    
    a <- cor.test(ic50 , ssgsea , type = "spearman")
    heatmap_data[heatmap_data$pathways%in%pathway[i]&heatmap_data$drugs%in%drug_final[j],3] <- a$estimate
    heatmap_data[heatmap_data$pathways%in%pathway[i]&heatmap_data$drugs%in%drug_final[j],4] <- a$p.value
    
  }
}
colnames(heatmap_data)[3] <- "R"
colnames(heatmap_data)[4] <- "p value"

# Replace the p-value with*
heatmap_data[,5] <- 
  ifelse(heatmap_data[,4]>0.05," ",ifelse(heatmap_data[,4]>0.01,"*",ifelse(heatmap_data[,4]>0.001,"**",ifelse(heatmap_data[,4]>0.0001,"***","****"))))
colnames(heatmap_data)[5] <- "p symbol"

# 这次用ggplot画这个热图
ggplot(heatmap_data, aes(pathways, drugs)) + 
  geom_tile(aes(fill = R), colour = "white",size=1)+
  scale_fill_gradient2(low = "#2b8cbe",mid = "white",high = "#e41a1c")+
  geom_text(aes(label=`p symbol`),col ="black",size = 5)+
  theme_minimal()+# no backgroud
  theme(axis.title.x=element_blank(),# no title
        axis.ticks.x=element_blank(),# no x-axis
        axis.title.y=element_blank(),# no y-axis
        axis.text.x = element_text(angle = 45, hjust = 1),# x-axis text
        axis.text.y = element_text(size = 8))+# y-axis text
  # adjust legend
  labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n","*** p < 0.001","\n\n","**** p < 0.0001","\n\n","Correlation"))+
  ggtitle("Overall CRC")+
  theme(plot.title = element_text(size = 17))# title size




# 6. msi heatmap ----------------------------------------------------------

heatmap_data_msi <- merge(pathway , drug_final)
colnames(heatmap_data_msi)[1] <- "pathways"
colnames(heatmap_data_msi)[2] <- "drugs"

for (i in 1:length(pathway)) {
  for (j in 1:length(drug_final)) {
    
    ssgsea_msi <- msi3_crc[msi3_crc$pathways%in%pathway[i]&msi3_crc$DRUG_NAME%in%drug_final[j],3]
    ic50_msi <- msi3_crc[msi3_crc$pathways%in%pathway[i]&msi3_crc$DRUG_NAME%in%drug_final[j],18]
    
    a_msi <- cor.test(ic50_msi , ssgsea_msi , type = "spearman")
    heatmap_data_msi[heatmap_data_msi$pathways%in%pathway[i]&heatmap_data_msi$drugs%in%drug_final[j],3] <- a_msi$estimate
    heatmap_data_msi[heatmap_data_msi$pathways%in%pathway[i]&heatmap_data_msi$drugs%in%drug_final[j],4] <- a_msi$p.value
    
  }
}
colnames(heatmap_data_msi)[3] <- "R"
colnames(heatmap_data_msi)[4] <- "p value"

heatmap_data_msi[,5] <- 
  ifelse(heatmap_data_msi[,4]>0.05," ",ifelse(heatmap_data_msi[,4]>0.01,"*",ifelse(heatmap_data_msi[,4]>0.001,"**",ifelse(heatmap_data_msi[,4]>0.0001,"***","****"))))
colnames(heatmap_data_msi)[5] <- "p symbol"

write.xlsx(heatmap_data_msi , file = "~/泛癌化疗药敏项目/result/supplementary_table/GDSC_correlation_heatmap_msi.xlsx")

ggplot(heatmap_data_msi, aes(pathways, drugs)) + 
  geom_tile(aes(fill = R), colour = "white",size=1)+
  scale_fill_gradient2(low = "#2b8cbe",mid = "white",high = "#e41a1c")+
  geom_text(aes(label=`p symbol`),col ="black",size = 5)+
  theme_minimal()+
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8))+
  
  labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n","*** p < 0.001","\n\n","**** p < 0.0001","\n\n","Correlation"))+
  ggtitle("MSI-H")+
  theme(plot.title = element_text(size = 17))
# 7. mss heatmap ----------------------------------------------------------

heatmap_data_mss <- merge(pathway , drug_final)
colnames(heatmap_data_mss)[1] <- "pathways"
colnames(heatmap_data_mss)[2] <- "drugs"

for (i in 1:length(pathway)) {
  for (j in 1:length(drug_final)) {
    
    ssgsea_mss <- mss3_crc[mss3_crc$pathways%in%pathway[i]&mss3_crc$DRUG_NAME%in%drug_final[j],3]
    ic50_mss <- mss3_crc[mss3_crc$pathways%in%pathway[i]&mss3_crc$DRUG_NAME%in%drug_final[j],18]
    
    a_mss <- cor.test(ic50_mss , ssgsea_mss , type = "spearman")
    heatmap_data_mss[heatmap_data_mss$pathways%in%pathway[i]&heatmap_data_mss$drugs%in%drug_final[j],3] <- a_mss$estimate
    heatmap_data_mss[heatmap_data_mss$pathways%in%pathway[i]&heatmap_data_mss$drugs%in%drug_final[j],4] <- a_mss$p.value
    
  }
}
colnames(heatmap_data_mss)[3] <- "R"
colnames(heatmap_data_mss)[4] <- "p value"

heatmap_data_mss[,5] <- 
  ifelse(heatmap_data_mss[,4]>0.05," ",ifelse(heatmap_data_mss[,4]>0.01,"*",ifelse(heatmap_data_mss[,4]>0.001,"**",ifelse(heatmap_data_mss[,4]>0.0001,"***","****"))))
colnames(heatmap_data_mss)[5] <- "p symbol"

write.xlsx(heatmap_data_mss , file = "~/泛癌化疗药敏项目/result/supplementary_table/GDSC_correlation_heatmap_mss.xlsx")

ggplot(heatmap_data_mss, aes(pathways, drugs)) + 
  geom_tile(aes(fill = R), colour = "white",size=1)+
  scale_fill_gradient2(low = "#2b8cbe",mid = "white",high = "#e41a1c")+
  geom_text(aes(label=`p symbol`),col ="black",size = 5)+
  theme_minimal()+
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8))+
  
  labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n","*** p < 0.001","\n\n","**** p < 0.0001","\n\n","Correlation"))+
  ggtitle("MSS/MSI-L")+
  theme(plot.title = element_text(size = 17))
