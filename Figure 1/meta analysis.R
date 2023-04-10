# 1. load R packages -----------------------------------------------------------------

library(GEOquery)
library(BiocManager)
library(affy)
library(pRRophetic)
library(org.Hs.eg.db)
library(rtracklayer)
library(plyr)
library(GenomicFeatures)
library(xlsx)
library(openxlsx)
library(hgu133plus2.db)
# BiocManager::install("hgu133a.db")
library(hgu133a.db)
# BiocManager::install("hu6800.db")
library(hu6800.db)
library(hgu133a.db)
library(data.table)
library(readxl)
library(clusterProfiler)
library(limma)
library(stringr)
library(GSVA)
library(GSEABase)
library(meta)
library(metafor)
library(forestplot)
# Error in windowsFonts(myFont1 = windowsFont("Arial")) : 
#   could not find function "windowsFonts"
library(ggplot2)
# GSE185055---------------------------------------------------------------
# 2. read data -----------------------------------------------------------------

Sys.setenv("VROOM_CONNECTION_SIZE"=999999999999999 * 99999999999999999)

GSE185055_names <- list.files(path = "/mnt/DATA/ytj/GEO_meta/GSE185055_RAW")

GSE185055_list <- lapply(GSE185055_names, function(x){read.table(
  file.path("/mnt/DATA/ytj/GEO_meta/GSE185055_RAW/", x), header = T)})

# rbind all targets
GSE185055_matrix <- do.call(cbind, GSE185055_list)

GSE185055_matrix <- GSE185055_matrix[,-seq(3,ncol(GSE185055_matrix),2)]
# change colnames to GSM symbol
GSE185055_colnames <- substring(GSE185055_names, 1, 10)
colnames(GSE185055_matrix)[2:ncol(GSE185055_matrix)] <- GSE185055_colnames



GSE185055_sample <- read.xlsx(
  "/mnt/DATA/ytj/GEO_meta/GSE185055_msi.xlsx", sheet = 1, 
)

GSE185055_sample1 <- GSE185055_sample[,-1]
# Manually annotate MSI status
GSE185055_sample1[1,] <- gsub("cell line: HT29 cell line", "MSS/MSI-L", GSE185055_sample1[1,])
GSE185055_sample1[1,] <- gsub("cell line: HCT116 cell line", "MSI-H", GSE185055_sample1[1,])
GSE185055_sample1[1,] <- gsub("cell line: LS513 cell line", "MSS/MSI-L", GSE185055_sample1[1,])
GSE185055_sample1[1,] <- gsub("cell line: LS174T cell line", "MSI-H", GSE185055_sample1[1,])

# devide into msi and mss
GSE185055_sammsi <- GSE185055_sample1[1,GSE185055_sample1[1,]%in%"MSI-H"]
GSE185055_sammss <- GSE185055_sample1[1,GSE185055_sample1[1,]%in%"MSS/MSI-L"]


# 3. data preprocessing -------------------------------------------------------------

# chemotherapy drugs
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

# transformation from count to tpm for comparison
# gene length file from ensembl
glong <- read.csv("/home/stu045/泛癌化疗药敏项目/Data/Homo_sapiens.GRCh38.108.csv")
colnames(glong)[1] <- "Id"

# Place gene length in the last column of the expression dataframe 
GSE185055_matrix1 <- join(GSE185055_matrix, glong, by = "Id")

# transformation from ensembl id to gene symbol
# gene annotation file from ensembl
ensembl <- rtracklayer::import('/home/stu045/泛癌化疗药敏项目/Data/Homo_sapiens.GRCh38.108.chr.gtf')   
ensmbl <- as.data.frame(ensembl)
ensmbl1 <- ensmbl[,c(10, 12)]
colnames(ensmbl1)[1] <- "Id"

GSE185055_matrix2 <-  join(GSE185055_matrix1, ensmbl1, type = "left", match = "first")
# Removing duplicate genes by comparing the median
which(is.na(GSE185055_matrix2$gene_name))
GSE185055_matrix2 <- GSE185055_matrix2[which(!is.na(GSE185055_matrix2$gene_name)),]
GSE185055_matrix2[,(ncol(GSE185055_matrix2)+1)] <- apply(
  GSE185055_matrix2[,2:(ncol(GSE185055_matrix2)-1)],1,median)
GSE185055_matrix3 <- GSE185055_matrix2[order(
  GSE185055_matrix2$gene_name,GSE185055_matrix2$V196,decreasing = T),]
GSE185055_matrix3 <- GSE185055_matrix3[!duplicated(GSE185055_matrix3$gene_name),]
rownames(GSE185055_matrix3) <- GSE185055_matrix3$gene_name
GSE185055_matrix4 <- GSE185055_matrix3[,-c(
  1, ncol(GSE185055_matrix3), (ncol(GSE185055_matrix3)-1))]

colnames(GSE185055_matrix4) <- GSE185055_colnames

# transformation to tpm
GSE185055_length <- GSE185055_matrix4[,ncol(GSE185055_matrix4)]/1000
GSE185055_matrix5 <- GSE185055_matrix4[,1:(ncol(GSE185055_matrix4)-1)] / GSE185055_length
GSE185055_tpm <- t(t(GSE185055_matrix5)/colSums(GSE185055_matrix5) * 10^6)

GSE185055_msi <- GSE185055_tpm[,colnames(GSE185055_sammsi)]
GSE185055_mss <- GSE185055_tpm[,colnames(GSE185055_sammss)]


# 4. chemotherapy prediction -----------------------------------------------------------------

GSE185055_predmsi <- data.frame(1:ncol(GSE185055_msi))
GSE185055_predmss <- data.frame(1:ncol(GSE185055_mss))

for (i in 1:length(drug)) {
  GSE185055_predmsi[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE185055_msi,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016"
  )
  colnames(GSE185055_predmsi)[i] <- drug[i]
}

for (i in 1:length(drug)) {
  GSE185055_predmss[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE185055_mss,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016"
  )
  colnames(GSE185055_predmss)[i] <- drug[i]
}


# 5. create table for meta -------------------------------------------------------------



length_na <- function(x){length(na.omit(x))}

GSE185055_predmsi1 <- as.data.frame(apply(GSE185055_predmsi, 2, function(x)(2^x)))
GSE185055_predmsi1[(nrow(GSE185055_predmsi1)+1),] <- apply(
  GSE185055_predmsi, 2, length_na)
GSE185055_predmsi1[(nrow(GSE185055_predmsi1)+1),] <- apply(
  GSE185055_predmsi1[1:(nrow(GSE185055_predmsi1)-1),], 2, mean, na.rm = T)
GSE185055_predmsi1[(nrow(GSE185055_predmsi1)+1),] <- apply(
  GSE185055_predmsi1[1:(nrow(GSE185055_predmsi1)-2),], 2, var, na.rm = T)
rownames(GSE185055_predmsi1)[(nrow(GSE185055_predmsi1)-2)] <- "n1"
rownames(GSE185055_predmsi1)[(nrow(GSE185055_predmsi1)-1)] <- "mean1"
rownames(GSE185055_predmsi1)[nrow(GSE185055_predmsi1)] <- "sd1"

GSE185055_predmss1 <- as.data.frame(apply(GSE185055_predmss, 2, function(x)(2^x)))
GSE185055_predmss1[(nrow(GSE185055_predmss1)+1),] <- apply(
  GSE185055_predmss, 2, length_na)
GSE185055_predmss1[(nrow(GSE185055_predmss1)+1),] <- apply(
  GSE185055_predmss1[1:(nrow(GSE185055_predmss1)-1),], 2, mean, na.rm = T)
GSE185055_predmss1[(nrow(GSE185055_predmss1)+1),] <- apply(
  GSE185055_predmss1[1:(nrow(GSE185055_predmss1)-2),], 2, var, na.rm = T)
rownames(GSE185055_predmss1)[(nrow(GSE185055_predmss1)-2)] <- "n2"
rownames(GSE185055_predmss1)[(nrow(GSE185055_predmss1)-1)] <- "mean2"
rownames(GSE185055_predmss1)[nrow(GSE185055_predmss1)] <- "sd2"
# GSE156915 ---------------------------------------------------------------
# 2. read data -----------------------------------------------------------------

GSE156915_geo <- read.table(
  "/mnt/DATA/ytj/GEO_meta/GSE156915_series_matrix.txt.gz", 
  sep="\t", header=T, comment.char="!")


GSE156915_ms <- read.xlsx2(
  "/mnt/DATA/ytj/GEO_meta/GSE156915_msi.xlsx", sheetIndex = 1, header = F
)

# 3. data preprocessing ----------------------------------------------------------------

GSE156915_matrix <- GSE156915_geo
rownames(GSE156915_matrix) <- GSE156915_matrix[,1]
GSE156915_matrix <- GSE156915_matrix[,-1]

GSE156915_sammsi <- GSE156915_ms[2,GSE156915_ms[1,]%in%"msi: MSI"]
GSE156915_sammss <- GSE156915_ms[2,GSE156915_ms[1,]%in%"msi: MSS"]

GSE156915_msi <- as.matrix(GSE156915_matrix[,colnames(GSE156915_matrix)%in%GSE156915_sammsi])
GSE156915_mss <- as.matrix(GSE156915_matrix[,colnames(GSE156915_matrix)%in%GSE156915_sammss])

# 4. chemotherapy prediction -----------------------------------------------------------------

GSE156915_predmsi <- data.frame(1:ncol(GSE156915_msi))
GSE156915_predmss <- data.frame(1:ncol(GSE156915_mss))

for (i in 1:length(drug)) {
  GSE156915_predmsi[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE156915_msi,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016"
  )
  colnames(GSE156915_predmsi)[i] <- drug[i]
}

for (i in 1:length(drug)) {
  GSE156915_predmss[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE156915_mss,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016"
  )
  colnames(GSE156915_predmss)[i] <- drug[i]
}

# 5. create table for meta -------------------------------------------------------------

GSE156915_predmsi1 <- as.data.frame(apply(GSE156915_predmsi, 2, function(x)(2^x)))
GSE156915_predmsi1[(nrow(GSE156915_predmsi1)+1),] <- apply(
  GSE156915_predmsi, 2, length_na)
GSE156915_predmsi1[(nrow(GSE156915_predmsi1)+1),] <- apply(
  GSE156915_predmsi1[1:(nrow(GSE156915_predmsi1)-1),], 2, mean, na.rm = T)
GSE156915_predmsi1[(nrow(GSE156915_predmsi1)+1),] <- apply(
  GSE156915_predmsi1[1:(nrow(GSE156915_predmsi1)-2),], 2, var, na.rm = T)
rownames(GSE156915_predmsi1)[(nrow(GSE156915_predmsi1)-2)] <- "n1"
rownames(GSE156915_predmsi1)[(nrow(GSE156915_predmsi1)-1)] <- "mean1"
rownames(GSE156915_predmsi1)[nrow(GSE156915_predmsi1)] <- "sd1"

GSE156915_predmss1 <- as.data.frame(apply(GSE156915_predmss, 2, function(x)(2^x)))
GSE156915_predmss1[(nrow(GSE156915_predmss1)+1),] <- apply(
  GSE156915_predmss, 2, length_na)
GSE156915_predmss1[(nrow(GSE156915_predmss1)+1),] <- apply(
  GSE156915_predmss1[1:(nrow(GSE156915_predmss1)-1),], 2, mean, na.rm = T)
GSE156915_predmss1[(nrow(GSE156915_predmss1)+1),] <- apply(
  GSE156915_predmss1[1:(nrow(GSE156915_predmss1)-2),], 2, var, na.rm = T)
rownames(GSE156915_predmss1)[(nrow(GSE156915_predmss1)-2)] <- "n2"
rownames(GSE156915_predmss1)[(nrow(GSE156915_predmss1)-1)] <- "mean2"
rownames(GSE156915_predmss1)[nrow(GSE156915_predmss1)] <- "sd2"
# GSE146889 ---------------------------------------------------------------
# 样本量小于10个，chemotherapy prediction报警说会很大影响结果，如果这份GEO数据放在meta里不好看就丢弃
# 2. read data -----------------------------------------------------------------

# Manually extract MSI information in Excel
GSE146889_sample <- read.xlsx2("/mnt/DATA/ytj/GEO_meta/GSE146889_msi.xlsx", sheetIndex = 1, header = F)

# expression profile
GSE146889_geo <- read.table("/mnt/DATA/ytj/GEO_meta/GSE146889_GeneCount.tsv.gz",
                            header = T)

# 3. data preprocessing ----------------------------------------------------------------

# extract the rpkm data
GSE146889_matrix <- GSE146889_geo[,c(6,seq(9,ncol(GSE146889_geo), by = 2))]

# because there are duplicate genes, it is also necessary to extract those with higher expression levels and calculate the median first
# ensure there is not null
which(is.na(GSE146889_matrix$GeneName))
# median
GSE146889_matrix[,(ncol(GSE146889_matrix)+1)] <- apply(
  GSE146889_matrix[,2:(ncol(GSE146889_matrix)-1)],1,median)
# order by median 
GSE146889_matrix1 <- GSE146889_matrix[order(
  GSE146889_matrix$GeneName,GSE146889_matrix$V178,decreasing = T),]
# only keep symbol with the highest expression level
GSE146889_matrix1 <- GSE146889_matrix1[!duplicated(GSE146889_matrix1$GeneName),]
rownames(GSE146889_matrix1) <- GSE146889_matrix1$GeneName
GSE146889_matrix1 <- GSE146889_matrix1[,-c(
  1, ncol(GSE146889_matrix1))]

# remove the rpkm in the colnames for the combination with sample data
GSE146889_matrix2 <- GSE146889_matrix1
GSE146889_matrix2[(nrow(GSE146889_matrix2)+1),] <- colnames(GSE146889_matrix2)
# strsplit() for spliting strings
GSE146889_matrix2[nrow(GSE146889_matrix2),] <- strsplit(
  as.character(GSE146889_matrix2[nrow(GSE146889_matrix2),]), "_rpkm")
colnames(GSE146889_matrix2) <- GSE146889_matrix2[nrow(GSE146889_matrix2),]
GSE146889_matrix3 <- GSE146889_matrix2[-nrow(GSE146889_matrix2),]

# msi and mss CRC sample
GSE146889_sample1 <- GSE146889_sample[,-1]
GSE146889_samcrc <- GSE146889_sample1[,GSE146889_sample1[1,]%in%"Colorectal"]
GSE146889_samtum <- GSE146889_samcrc[,GSE146889_samcrc[3,]%in%"phenotype(normal/tumor): tumor"]
GSE146889_sammsi <- GSE146889_samtum[4, GSE146889_samtum[2,]%in%"phenotype (msi-h or mss): MSI"]
GSE146889_sammss <- GSE146889_samtum[4, GSE146889_samtum[2,]%in%"phenotype (msi-h or mss): MSS"]

# msi and mss expression
GSE146889_msirpkm <- apply(
  GSE146889_matrix3[,colnames(GSE146889_matrix3)%in%GSE146889_sammsi], 2, as.numeric)
GSE146889_mssrpkm <- apply(
  GSE146889_matrix3[,colnames(GSE146889_matrix3)%in%GSE146889_sammss], 2, as.numeric)
rownames(GSE146889_msirpkm) <- rownames(GSE146889_matrix3)
rownames(GSE146889_mssrpkm) <- rownames(GSE146889_matrix3)

# transformation of standardization from rpkm to tpm
rpkmTOtpm <- function(rpkm){
  exp(log(rpkm) - log(sum(rpkm)) + log(1e6))
} 
GSE146889_msi <- apply(GSE146889_msirpkm, 2, rpkmTOtpm)
GSE146889_mss <- apply(GSE146889_mssrpkm, 2, rpkmTOtpm)



# 4. chemotherapy prediction -----------------------------------------------------------------

GSE146889_predmsi <- data.frame(1:ncol(GSE146889_msi))
GSE146889_predmss <- data.frame(1:ncol(GSE146889_mss))

for (i in 1:length(drug)) {
  GSE146889_predmsi[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE146889_msi,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE146889_predmsi)[i] <- drug[i]
}

for (i in 1:length(drug)) {
  GSE146889_predmss[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE146889_mss,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE146889_predmss)[i] <- drug[i]
}

# 5. create table for meta -------------------------------------------------------------

GSE146889_predmsi1 <- as.data.frame(apply(GSE146889_predmsi, 2, function(x)(2^x)))
GSE146889_predmsi1[(nrow(GSE146889_predmsi1)+1),] <- apply(
  GSE146889_predmsi, 2, length_na)
GSE146889_predmsi1[(nrow(GSE146889_predmsi1)+1),] <- apply(
  GSE146889_predmsi1[1:(nrow(GSE146889_predmsi1)-1),], 2, mean, na.rm = T)
GSE146889_predmsi1[(nrow(GSE146889_predmsi1)+1),] <- apply(
  GSE146889_predmsi1[1:(nrow(GSE146889_predmsi1)-2),], 2, var, na.rm = T)
rownames(GSE146889_predmsi1)[(nrow(GSE146889_predmsi1)-2)] <- "n1"
rownames(GSE146889_predmsi1)[(nrow(GSE146889_predmsi1)-1)] <- "mean1"
rownames(GSE146889_predmsi1)[nrow(GSE146889_predmsi1)] <- "sd1"

GSE146889_predmss1 <- as.data.frame(apply(GSE146889_predmss, 2, function(x)(2^x)))
GSE146889_predmss1[(nrow(GSE146889_predmss1)+1),] <- apply(
  GSE146889_predmss, 2, length_na)
GSE146889_predmss1[(nrow(GSE146889_predmss1)+1),] <- apply(
  GSE146889_predmss1[1:(nrow(GSE146889_predmss1)-1),], 2, mean, na.rm = T)
GSE146889_predmss1[(nrow(GSE146889_predmss1)+1),] <- apply(
  GSE146889_predmss1[1:(nrow(GSE146889_predmss1)-2),], 2, var, na.rm = T)
rownames(GSE146889_predmss1)[(nrow(GSE146889_predmss1)-2)] <- "n2"
rownames(GSE146889_predmss1)[(nrow(GSE146889_predmss1)-1)] <- "mean2"
rownames(GSE146889_predmss1)[nrow(GSE146889_predmss1)] <- "sd2"
# GSE143985 ---------------------------------------------------------------
# 2. read data -----------------------------------------------------------------


GSE143985_geo <- read.table("/mnt/DATA/ytj/GEO_meta/GSE143985_series_matrix.txt.gz",
                            sep = "\t", comment.char = "!", header = T)

GSE143985_sample <- read.xlsx2("/mnt/DATA/ytj/GEO_meta/GSE143985_msi.xlsx",
                               sheetIndex = 1, header = F)

# 3. data preprocessing ----------------------------------------------------------------

range(GSE143985_geo[,2:ncol(GSE143985_geo)])

# change probe id 
at_ids <- toTable(hgu133plus2SYMBOL)
GSE143985_matrix <- GSE143985_geo
colnames(GSE143985_matrix)[1] <- "probe_id"
GSE143985_matrix1 <- merge(GSE143985_matrix, at_ids, by.x="probe_id", by.y="probe_id")

which(duplicated(GSE143985_matrix1[,ncol(GSE143985_matrix1)]))
# median for removing of duplicated gene id
GSE143985_matrix2 <- GSE143985_matrix1
GSE143985_matrix2[,(ncol(GSE143985_matrix2)+1)] <- apply(
  GSE143985_matrix2[,2:(ncol(GSE143985_matrix2)-1)], 1, median)
GSE143985_matrix2 <- GSE143985_matrix2[order(
  GSE143985_matrix2$symbol, GSE143985_matrix2$V94, decreasing = T),]
GSE143985_matrix2 <- GSE143985_matrix2[which(!duplicated(GSE143985_matrix2$symbol)),]

rownames(GSE143985_matrix2) <- GSE143985_matrix2$symbol
GSE143985_matrix2 <- GSE143985_matrix2[,-c(1, ncol(GSE143985_matrix2), (ncol(GSE143985_matrix2)-1))]

# log
GSE143985_matrix3 <- apply(GSE143985_matrix2, 2 , log2)

# msi and mss
GSE143985_sample1 <- GSE143985_sample[,-1]
GSE143985_sammsi <- GSE143985_sample1[1,GSE143985_sample1[2,]%in%"msi_status: MSI_H"]
GSE143985_sammss <- GSE143985_sample1[1,GSE143985_sample1[2,]%in%"msi_status: MSS"]
GSE143985_msi <- GSE143985_matrix3[,colnames(GSE143985_matrix3)%in%GSE143985_sammsi]
GSE143985_mss <- GSE143985_matrix3[,colnames(GSE143985_matrix3)%in%GSE143985_sammss]

# 4. chemotherapy prediction -----------------------------------------------------------------

GSE143985_predmsi <- data.frame(1:ncol(GSE143985_msi))
GSE143985_predmss <- data.frame(1:ncol(GSE143985_mss))

for (i in 1:length(drug)) {
  GSE143985_predmsi[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE143985_msi,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE143985_predmsi)[i] <- drug[i]
}

for (i in 1:length(drug)) {
  GSE143985_predmss[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE143985_mss,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE143985_predmss)[i] <- drug[i]
}

# 5. create table for meta -------------------------------------------------------------

GSE143985_predmsi1 <- as.data.frame(apply(GSE143985_predmsi, 2, function(x)(2^x)))
GSE143985_predmsi1[(nrow(GSE143985_predmsi1)+1),] <- apply(
  GSE143985_predmsi, 2, length_na)
GSE143985_predmsi1[(nrow(GSE143985_predmsi1)+1),] <- apply(
  GSE143985_predmsi1[1:(nrow(GSE143985_predmsi1)-1),], 2, mean, na.rm = T)
GSE143985_predmsi1[(nrow(GSE143985_predmsi1)+1),] <- apply(
  GSE143985_predmsi1[1:(nrow(GSE143985_predmsi1)-2),], 2, var, na.rm = T)
rownames(GSE143985_predmsi1)[(nrow(GSE143985_predmsi1)-2)] <- "n1"
rownames(GSE143985_predmsi1)[(nrow(GSE143985_predmsi1)-1)] <- "mean1"
rownames(GSE143985_predmsi1)[nrow(GSE143985_predmsi1)] <- "sd1"

GSE143985_predmss1 <- as.data.frame(apply(GSE143985_predmss, 2, function(x)(2^x)))
GSE143985_predmss1[(nrow(GSE143985_predmss1)+1),] <- apply(
  GSE143985_predmss, 2, length_na)
GSE143985_predmss1[(nrow(GSE143985_predmss1)+1),] <- apply(
  GSE143985_predmss1[1:(nrow(GSE143985_predmss1)-1),], 2, mean, na.rm = T)
GSE143985_predmss1[(nrow(GSE143985_predmss1)+1),] <- apply(
  GSE143985_predmss1[1:(nrow(GSE143985_predmss1)-2),], 2, var, na.rm = T)
rownames(GSE143985_predmss1)[(nrow(GSE143985_predmss1)-2)] <- "n2"
rownames(GSE143985_predmss1)[(nrow(GSE143985_predmss1)-1)] <- "mean2"
rownames(GSE143985_predmss1)[nrow(GSE143985_predmss1)] <- "sd2"
# GSE103340 ---------------------------------------------------------------
# 2. read data -----------------------------------------------------------------

GSE103340_sample <- read.xlsx2("/mnt/DATA/ytj/GEO_meta/GSE103340_msi.xlsx", sheetIndex = 1)

GSE103340_geo <- read.table("/mnt/DATA/ytj/GEO_meta/GSE103340_series_matrix.txt.gz",
                            header = T, comment.char = "!")

# 3. data preprocessing ----------------------------------------------------------------

# According to the series matrix file, the expression profile has been standardized by RMA and can be directly used for differential analysis

# download platform annotation file and extract two needed row from it manually
GPL13158 <- read.xlsx2("/mnt/DATA/ytj/GEO_meta/GPL13158.xlsx", sheetIndex = 1)

# change gene id
GSE103340_matrix <- merge(GSE103340_geo , GPL13158, by.x = "ID_REF", by.y = "ID")
which(duplicated(GSE103340_matrix[,ncol(GSE103340_matrix)]))

GSE103340_matrix[,(ncol(GSE103340_matrix)+1)] <- apply(
  GSE103340_matrix[,2:(ncol(GSE103340_matrix)-1)], 1, median)
GSE103340_matrix <- GSE103340_matrix[order(
  GSE103340_matrix$Gene.Symbol, GSE103340_matrix$V74, decreasing = T),]
GSE103340_matrix <- GSE103340_matrix[which(!duplicated(GSE103340_matrix$Gene.Symbol)),]
rownames(GSE103340_matrix) <- GSE103340_matrix$Gene.Symbol
GSE103340_matrix1 <- GSE103340_matrix[,-c(
  1,ncol(GSE103340_matrix),(ncol(GSE103340_matrix)-1))]


GSE103340_sample1 <- GSE103340_sample[,-1]
GSE103340_sammsi <- GSE103340_sample1[3,GSE103340_sample1[2,]%in%"msi: 1"]
GSE103340_sammss <- GSE103340_sample1[3,GSE103340_sample1[2,]%in%"msi: 0"]
GSE103340_msi <- as.matrix(GSE103340_matrix1[,colnames(GSE103340_matrix1)%in%GSE103340_sammsi])
GSE103340_mss <- as.matrix(GSE103340_matrix1[,colnames(GSE103340_matrix1)%in%GSE103340_sammss])


# 4. chemotherapy prediction -----------------------------------------------------------------


GSE103340_predmsi <- data.frame(1:ncol(GSE103340_msi))
GSE103340_predmss <- data.frame(1:ncol(GSE103340_mss))

for (i in 1:length(drug)) {
  GSE103340_predmsi[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE103340_msi,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE103340_predmsi)[i] <- drug[i]
}

for (i in 1:length(drug)) {
  GSE103340_predmss[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE103340_mss,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE103340_predmss)[i] <- drug[i]
}

# 5. create table for meta -------------------------------------------------------------

GSE103340_predmsi1 <- as.data.frame(apply(GSE103340_predmsi, 2, function(x)(2^x)))
GSE103340_predmsi1[(nrow(GSE103340_predmsi1)+1),] <- apply(
  GSE103340_predmsi, 2, length_na)
GSE103340_predmsi1[(nrow(GSE103340_predmsi1)+1),] <- apply(
  GSE103340_predmsi1[1:(nrow(GSE103340_predmsi1)-1),], 2, mean, na.rm = T)
GSE103340_predmsi1[(nrow(GSE103340_predmsi1)+1),] <- apply(
  GSE103340_predmsi1[1:(nrow(GSE103340_predmsi1)-2),], 2, var, na.rm = T)
rownames(GSE103340_predmsi1)[(nrow(GSE103340_predmsi1)-2)] <- "n1"
rownames(GSE103340_predmsi1)[(nrow(GSE103340_predmsi1)-1)] <- "mean1"
rownames(GSE103340_predmsi1)[nrow(GSE103340_predmsi1)] <- "sd1"

GSE103340_predmss1 <- as.data.frame(apply(GSE103340_predmss, 2, function(x)(2^x)))
GSE103340_predmss1[(nrow(GSE103340_predmss1)+1),] <- apply(
  GSE103340_predmss, 2, length_na)
GSE103340_predmss1[(nrow(GSE103340_predmss1)+1),] <- apply(
  GSE103340_predmss1[1:(nrow(GSE103340_predmss1)-1),], 2, mean, na.rm = T)
GSE103340_predmss1[(nrow(GSE103340_predmss1)+1),] <- apply(
  GSE103340_predmss1[1:(nrow(GSE103340_predmss1)-2),], 2, var, na.rm = T)
rownames(GSE103340_predmss1)[(nrow(GSE103340_predmss1)-2)] <- "n2"
rownames(GSE103340_predmss1)[(nrow(GSE103340_predmss1)-1)] <- "mean2"
rownames(GSE103340_predmss1)[nrow(GSE103340_predmss1)] <- "sd2"
# GSE75317 ----------------------------------------------------------------
# 2. read data -----------------------------------------------------------------

GSE75317_geo <- read.table("/mnt/DATA/ytj/GEO_meta/GSE75317-GPL570_series_matrix.txt.gz",
                           sep = "\t", header = T, comment.char = "!")
GSE75317_sample <- read.xlsx2("/mnt/DATA/ytj/GEO_meta/GSE75317-GPL570_msi.xlsx",
                              sheetIndex = 1)

# 3. data preprocessing ----------------------------------------------------------------


GSE75317_matrix <- merge(GSE75317_geo, at_ids, by.x="ID_REF", by.y="probe_id")
which(duplicated(GSE75317_matrix$symbol))

GSE75317_matrix1 <- GSE75317_matrix
GSE75317_matrix1[,(ncol(GSE75317_matrix1)+1)] <- apply(
  GSE75317_matrix1[,2:(ncol(GSE75317_matrix1)-1)], 1, median)

GSE75317_matrix1 <- GSE75317_matrix1[order(
  GSE75317_matrix1$symbol, GSE75317_matrix1$V62, decreasing = T),]
GSE75317_matrix1 <- GSE75317_matrix1[which(!duplicated(GSE75317_matrix1$symbol)),]
rownames(GSE75317_matrix1) <- GSE75317_matrix1[,(ncol(GSE75317_matrix1)-1)]
GSE75317_matrix1 <- GSE75317_matrix1[,-c(1, ncol(GSE75317_matrix1), (ncol(GSE75317_matrix1)-1))]

GSE75317_sample1 <- GSE75317_sample[,-1]
GSE75317_sammsi <- colnames(GSE75317_sample1[,GSE75317_sample1%in%"microsatellite status: MSI-H"])
GSE75317_sammss <- colnames(GSE75317_sample1[,GSE75317_sample1%in%"microsatellite status: MSS"])

GSE75317_msi <- as.matrix(GSE75317_matrix1[,colnames(GSE75317_matrix1)%in%GSE75317_sammsi])
GSE75317_mss <- as.matrix(GSE75317_matrix1[,colnames(GSE75317_matrix1)%in%GSE75317_sammss])


# 4. chemotherapy prediction -----------------------------------------------------------------


GSE75317_predmsi <- data.frame(1:ncol(GSE75317_msi))
GSE75317_predmss <- data.frame(1:ncol(GSE75317_mss))

for (i in 1:length(drug)) {
  GSE75317_predmsi[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE75317_msi,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE75317_predmsi)[i] <- drug[i]
}

for (i in 1:length(drug)) {
  GSE75317_predmss[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE75317_mss,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE75317_predmss)[i] <- drug[i]
}

# 5. create table for meta -------------------------------------------------------------

GSE75317_predmsi1 <- as.data.frame(apply(GSE75317_predmsi, 2, function(x)(2^x)))
GSE75317_predmsi1[(nrow(GSE75317_predmsi1)+1),] <- apply(
  GSE75317_predmsi, 2, length_na)
GSE75317_predmsi1[(nrow(GSE75317_predmsi1)+1),] <- apply(
  GSE75317_predmsi1[1:(nrow(GSE75317_predmsi1)-1),], 2, mean, na.rm = T)
GSE75317_predmsi1[(nrow(GSE75317_predmsi1)+1),] <- apply(
  GSE75317_predmsi1[1:(nrow(GSE75317_predmsi1)-2),], 2, var, na.rm = T)
rownames(GSE75317_predmsi1)[(nrow(GSE75317_predmsi1)-2)] <- "n1"
rownames(GSE75317_predmsi1)[(nrow(GSE75317_predmsi1)-1)] <- "mean1"
rownames(GSE75317_predmsi1)[nrow(GSE75317_predmsi1)] <- "sd1"

GSE75317_predmss1 <- as.data.frame(apply(GSE75317_predmss, 2, function(x)(2^x)))
GSE75317_predmss1[(nrow(GSE75317_predmss1)+1),] <- apply(
  GSE75317_predmss, 2, length_na)
GSE75317_predmss1[(nrow(GSE75317_predmss1)+1),] <- apply(
  GSE75317_predmss1[1:(nrow(GSE75317_predmss1)-1),], 2, mean, na.rm = T)
GSE75317_predmss1[(nrow(GSE75317_predmss1)+1),] <- apply(
  GSE75317_predmss1[1:(nrow(GSE75317_predmss1)-2),], 2, var, na.rm = T)
rownames(GSE75317_predmss1)[(nrow(GSE75317_predmss1)-2)] <- "n2"
rownames(GSE75317_predmss1)[(nrow(GSE75317_predmss1)-1)] <- "mean2"
rownames(GSE75317_predmss1)[nrow(GSE75317_predmss1)] <- "sd2"
# GSE75315 ----------------------------------------------------------------
# 2. read data -----------------------------------------------------------------

GSE75315_geo <- read.table("/mnt/DATA/ytj/GEO_meta/GSE75315_series_matrix.txt.gz",
                           sep = "\t", header = T, comment.char = "!")
GSE75315_sample <- read.xlsx2("/mnt/DATA/ytj/GEO_meta/GSE75315_msi.xlsx",
                              sheetIndex = 1)

# 3. data preprocessing ----------------------------------------------------------------

GPL5175 <- fread("/mnt/DATA/ytj/GEO_meta/GPL5175-3188.txt")
# Obtain gene annotation information for these 17565 probes
GSE75315_gene <- subset(
  GPL5175, GPL5175$ID%in%GSE75315_geo[,1], select=c(ID,gene_assignment))

GSE75315_gene_anno <- as.data.frame(
  str_split_fixed(GSE75315_gene$gene_assignment, " // ", n = 3))
GSE75315_gene_anno2 <- as.data.frame(
  cbind(GSE75315_gene$ID, GSE75315_gene_anno$V2))

colnames(GSE75315_gene_anno2) <- c("ID","symbol")
GSE75315_matrix <- merge(
  GSE75315_geo, GSE75315_gene_anno2, by.x = "ID_REF", by.y = "ID")

which(duplicated(GSE75315_matrix$symbol))

# median
GSE75315_matrix1 <- GSE75315_matrix
GSE75315_matrix1[,(ncol(GSE75315_matrix1)+1)] <- apply(
  GSE75315_matrix1[,2:(ncol(GSE75315_matrix1)-1)], 1, median)
GSE75315_matrix1 <- GSE75315_matrix1[order(
  GSE75315_matrix1$symbol, GSE75315_matrix1$V214, decreasing = T),]
GSE75315_matrix1 <- GSE75315_matrix1[which(!duplicated(GSE75315_matrix1$symbol)),]
rownames(GSE75315_matrix1) <- GSE75315_matrix1[,(ncol(GSE75315_matrix1)-1)]
GSE75315_matrix1 <- GSE75315_matrix1[,-c(1, ncol(GSE75315_matrix1), (ncol(GSE75315_matrix1)-1))]

GSE75315_sample1 <- GSE75315_sample[,-1]
GSE75315_sammsi <- colnames(GSE75315_sample1[,GSE75315_sample1%in%"microsatellite status: MSI-H"])
GSE75315_sammss <- colnames(GSE75315_sample1[,GSE75315_sample1%in%"microsatellite status: MSS"])

GSE75315_msi <- as.matrix(GSE75315_matrix1[,colnames(GSE75315_matrix1)%in%GSE75315_sammsi])
GSE75315_mss <- as.matrix(GSE75315_matrix1[,colnames(GSE75315_matrix1)%in%GSE75315_sammss])

# 4. chemotherapy prediction -----------------------------------------------------------------

GSE75315_predmsi <- data.frame(1:ncol(GSE75315_msi))
GSE75315_predmss <- data.frame(1:ncol(GSE75315_mss))

 
for (i in 1:length(drug)) {
  GSE75315_predmsi[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE75315_msi,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE75315_predmsi)[i] <- drug[i]
}

for (i in 1:length(drug)) {
  GSE75315_predmss[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE75315_mss,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE75315_predmss)[i] <- drug[i]
}

# 5. create table for meta -------------------------------------------------------------

GSE75315_predmsi1 <- as.data.frame(apply(GSE75315_predmsi, 2, function(x)(2^x)))
GSE75315_predmsi1[(nrow(GSE75315_predmsi1)+1),] <- apply(
  GSE75315_predmsi, 2, length_na)
GSE75315_predmsi1[(nrow(GSE75315_predmsi1)+1),] <- apply(
  GSE75315_predmsi1[1:(nrow(GSE75315_predmsi1)-1),], 2, mean, na.rm = T)
GSE75315_predmsi1[(nrow(GSE75315_predmsi1)+1),] <- apply(
  GSE75315_predmsi1[1:(nrow(GSE75315_predmsi1)-2),], 2, var, na.rm = T)
rownames(GSE75315_predmsi1)[(nrow(GSE75315_predmsi1)-2)] <- "n1"
rownames(GSE75315_predmsi1)[(nrow(GSE75315_predmsi1)-1)] <- "mean1"
rownames(GSE75315_predmsi1)[nrow(GSE75315_predmsi1)] <- "sd1"

GSE75315_predmss1 <- as.data.frame(apply(GSE75315_predmss, 2, function(x)(2^x)))
GSE75315_predmss1[(nrow(GSE75315_predmss1)+1),] <- apply(
  GSE75315_predmss, 2, length_na)
GSE75315_predmss1[(nrow(GSE75315_predmss1)+1),] <- apply(
  GSE75315_predmss1[1:(nrow(GSE75315_predmss1)-1),], 2, mean, na.rm = T)
GSE75315_predmss1[(nrow(GSE75315_predmss1)+1),] <- apply(
  GSE75315_predmss1[1:(nrow(GSE75315_predmss1)-2),], 2, var, na.rm = T)
rownames(GSE75315_predmss1)[(nrow(GSE75315_predmss1)-2)] <- "n2"
rownames(GSE75315_predmss1)[(nrow(GSE75315_predmss1)-1)] <- "mean2"
rownames(GSE75315_predmss1)[nrow(GSE75315_predmss1)] <- "sd2"
# GSE39084 ----------------------------------------------------------------
# 2. read data -----------------------------------------------------------------

GSE39084_geo <- read.table("/mnt/DATA/ytj/GEO_meta/GSE39084_series_matrix.txt",
                           header = T, comment.char = "!", sep = "\t")
GSE39084_sample <- read.xlsx2("/mnt/DATA/ytj/GEO_meta/GSE39084_msi.xlsx",
                              sheetIndex = 1)

# 3. data preprocessing ----------------------------------------------------------------

GSE39084_matrix <- merge(GSE39084_geo, at_ids, by.x = "ID_REF", by.y = "probe_id")

which(duplicated(GSE39084_matrix[,ncol(GSE39084_matrix)]))

# median
GSE39084_matrix1 <- GSE39084_matrix
GSE39084_matrix1[,(ncol(GSE39084_matrix1)+1)] <- apply(
  GSE39084_matrix1[,2:(ncol(GSE39084_matrix1)-1)], 1, median)
GSE39084_matrix1 <- GSE39084_matrix1[order(
  GSE39084_matrix1$symbol, GSE39084_matrix1$V73, decreasing = T),]
GSE39084_matrix1 <- GSE39084_matrix1[which(!duplicated(GSE39084_matrix1$symbol)),]

rownames(GSE39084_matrix1) <- GSE39084_matrix1$symbol
GSE39084_matrix1 <- GSE39084_matrix1[,-c(1, ncol(GSE39084_matrix1), (ncol(GSE39084_matrix1)-1))]

# msi and mss
GSE39084_sample1 <- GSE39084_sample[,-1]
GSE39084_sammsi <- colnames(
  GSE39084_sample1[,GSE39084_sample1[1,]%in%"msi.status [mismatch repair status]: High"])
GSE39084_sammss <- colnames(
  GSE39084_sample1[,GSE39084_sample1[1,]%in%c(
    "msi.status [mismatch repair status]: No", "msi.status [mismatch repair status]: Low")])
GSE39084_msi <- as.matrix(GSE39084_matrix1[,colnames(GSE39084_matrix1)%in%GSE39084_sammsi])
GSE39084_mss <- as.matrix(GSE39084_matrix1[,colnames(GSE39084_matrix1)%in%GSE39084_sammss])

# 4. chemotherapy prediction -----------------------------------------------------------------

GSE39084_predmsi <- data.frame(1:ncol(GSE39084_msi))
GSE39084_predmss <- data.frame(1:ncol(GSE39084_mss))

for (i in 1:length(drug)) {
  GSE39084_predmsi[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE39084_msi,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE39084_predmsi)[i] <- drug[i]
}

for (i in 1:length(drug)) {
  GSE39084_predmss[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE39084_mss,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE39084_predmss)[i] <- drug[i]
}

# 5. create table for meta -------------------------------------------------------------

GSE39084_predmsi1 <- as.data.frame(apply(GSE39084_predmsi, 2, function(x)(2^x)))
GSE39084_predmsi1[(nrow(GSE39084_predmsi1)+1),] <- apply(
  GSE39084_predmsi, 2, length_na)
GSE39084_predmsi1[(nrow(GSE39084_predmsi1)+1),] <- apply(
  GSE39084_predmsi1[1:(nrow(GSE39084_predmsi1)-1),], 2, mean, na.rm = T)
GSE39084_predmsi1[(nrow(GSE39084_predmsi1)+1),] <- apply(
  GSE39084_predmsi1[1:(nrow(GSE39084_predmsi1)-2),], 2, var, na.rm = T)
rownames(GSE39084_predmsi1)[(nrow(GSE39084_predmsi1)-2)] <- "n1"
rownames(GSE39084_predmsi1)[(nrow(GSE39084_predmsi1)-1)] <- "mean1"
rownames(GSE39084_predmsi1)[nrow(GSE39084_predmsi1)] <- "sd1"

GSE39084_predmss1 <- as.data.frame(apply(GSE39084_predmss, 2, function(x)(2^x)))
GSE39084_predmss1[(nrow(GSE39084_predmss1)+1),] <- apply(
  GSE39084_predmss, 2, length_na)
GSE39084_predmss1[(nrow(GSE39084_predmss1)+1),] <- apply(
  GSE39084_predmss1[1:(nrow(GSE39084_predmss1)-1),], 2, mean, na.rm = T)
GSE39084_predmss1[(nrow(GSE39084_predmss1)+1),] <- apply(
  GSE39084_predmss1[1:(nrow(GSE39084_predmss1)-2),], 2, var, na.rm = T)
rownames(GSE39084_predmss1)[(nrow(GSE39084_predmss1)-2)] <- "n2"
rownames(GSE39084_predmss1)[(nrow(GSE39084_predmss1)-1)] <- "mean2"
rownames(GSE39084_predmss1)[nrow(GSE39084_predmss1)] <- "sd2"
# GSE35896 ------------------------------------------------------------
# 2. read data -----------------------------------------------------------------

GSE35896_geo <- read.table("/mnt/DATA/ytj/GEO_meta/GSE35896_series_matrix.txt.gz",
                           header = T, comment.char = "!", sep = "\t")

GSE35896_sample <- read.xlsx2("/mnt/DATA/ytj/GEO_meta/GSE35896_msi.xlsx", sheetIndex = 1,
                              header = F)

# 3. data preprocessing ----------------------------------------------------------------

GSE35896_matrix <- merge(GSE35896_geo, at_ids, by.x = "ID_REF", by.y = "probe_id")

which(duplicated(GSE35896_matrix[,ncol(GSE35896_matrix)]))
# mdian
GSE35896_matrix1 <- GSE35896_matrix
GSE35896_matrix1[,(ncol(GSE35896_matrix1)+1)] <- apply(
  GSE35896_matrix1[,2:(ncol(GSE35896_matrix1)-1)], 1, median)
GSE35896_matrix1 <- GSE35896_matrix1[order(
  GSE35896_matrix1$symbol, GSE35896_matrix1$V65, decreasing = T),]
GSE35896_matrix1 <- GSE35896_matrix1[which(!duplicated(GSE35896_matrix1$symbol)),]

rownames(GSE35896_matrix1) <- GSE35896_matrix1$symbol
GSE35896_matrix1 <- GSE35896_matrix1[,-c(1, ncol(GSE35896_matrix1), (ncol(GSE35896_matrix1)-1))]

# msi and mss
GSE35896_sample1 <- GSE35896_sample[,-1]
GSE35896_sammsi <- GSE35896_sample1[2,GSE35896_sample1[1,]%in%"microsatellite.status: MSI"]
GSE35896_sammss <- GSE35896_sample1[2,GSE35896_sample1[1,]%in%"microsatellite.status: MSS"]
GSE35896_msi <- as.matrix(GSE35896_matrix1[,colnames(GSE35896_matrix1)%in%GSE35896_sammsi])
GSE35896_mss <- as.matrix(GSE35896_matrix1[,colnames(GSE35896_matrix1)%in%GSE35896_sammss])

# 4. chemotherapy prediction -----------------------------------------------------------------

GSE35896_predmsi <- data.frame(1:ncol(GSE35896_msi))
GSE35896_predmss <- data.frame(1:ncol(GSE35896_mss))

 
for (i in 1:length(drug)) {
  GSE35896_predmsi[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE35896_msi,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE35896_predmsi)[i] <- drug[i]
}

for (i in 1:length(drug)) {
  GSE35896_predmss[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE35896_mss,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE35896_predmss)[i] <- drug[i]
}

# 5. create table for meta -------------------------------------------------------------

GSE35896_predmsi1 <- as.data.frame(apply(GSE35896_predmsi, 2, function(x)(2^x)))
GSE35896_predmsi1[(nrow(GSE35896_predmsi1)+1),] <- apply(
  GSE35896_predmsi, 2, length_na)
GSE35896_predmsi1[(nrow(GSE35896_predmsi1)+1),] <- apply(
  GSE35896_predmsi1[1:(nrow(GSE35896_predmsi1)-1),], 2, mean, na.rm = T)
GSE35896_predmsi1[(nrow(GSE35896_predmsi1)+1),] <- apply(
  GSE35896_predmsi1[1:(nrow(GSE35896_predmsi1)-2),], 2, var, na.rm = T)
rownames(GSE35896_predmsi1)[(nrow(GSE35896_predmsi1)-2)] <- "n1"
rownames(GSE35896_predmsi1)[(nrow(GSE35896_predmsi1)-1)] <- "mean1"
rownames(GSE35896_predmsi1)[nrow(GSE35896_predmsi1)] <- "sd1"

GSE35896_predmss1 <- as.data.frame(apply(GSE35896_predmss, 2, function(x)(2^x)))
GSE35896_predmss1[(nrow(GSE35896_predmss1)+1),] <- apply(
  GSE35896_predmss, 2, length_na)
GSE35896_predmss1[(nrow(GSE35896_predmss1)+1),] <- apply(
  GSE35896_predmss1[1:(nrow(GSE35896_predmss1)-1),], 2, mean, na.rm = T)
GSE35896_predmss1[(nrow(GSE35896_predmss1)+1),] <- apply(
  GSE35896_predmss1[1:(nrow(GSE35896_predmss1)-2),], 2, var, na.rm = T)
rownames(GSE35896_predmss1)[(nrow(GSE35896_predmss1)-2)] <- "n2"
rownames(GSE35896_predmss1)[(nrow(GSE35896_predmss1)-1)] <- "mean2"
rownames(GSE35896_predmss1)[nrow(GSE35896_predmss1)] <- "sd2"
# GSE41258 ----------------------------------------------------------------
# 2. read data -----------------------------------------------------------------

GSE41258_geo <- read.table("/mnt/DATA/ytj/GEO_meta/GSE41258_series_matrix.txt.gz",
                           header = T, comment.char = "!", sep = "\t")

GSE41258_sample <- read.xlsx2("/mnt/DATA/ytj/GEO_meta/GSE41258_msi.xlsx", sheetIndex = 1,
                              header = F)

# 3. data preprocessing ----------------------------------------------------------------

# The expression spectrum should be standardized because it is not an integer, but the numerical value is a bit large and should not have been log transformed
# GPL96
at_ids96 <- toTable(hgu133aSYMBOL)
GSE41258_matrix <- merge(GSE41258_geo, at_ids, by.x = "ID_REF", by.y = "probe_id")

which(duplicated(GSE41258_matrix[,ncol(GSE41258_matrix)]))

GSE41258_matrix1 <- GSE41258_matrix
GSE41258_matrix1[,(ncol(GSE41258_matrix1)+1)] <- apply(
  GSE41258_matrix1[,2:(ncol(GSE41258_matrix1)-1)], 1, median)
GSE41258_matrix1 <- GSE41258_matrix1[order(
  GSE41258_matrix1$symbol, GSE41258_matrix1$V393, decreasing = T),]
GSE41258_matrix1 <- GSE41258_matrix1[which(!duplicated(GSE41258_matrix1$symbol)),]

rownames(GSE41258_matrix1) <- GSE41258_matrix1$symbol
GSE41258_matrix1 <- GSE41258_matrix1[,-c(1, ncol(GSE41258_matrix1), (ncol(GSE41258_matrix1)-1))]

GSE41258_matrix2 <- t(apply(GSE41258_matrix1, 1, log2))

# extract tumor sample
# devide into msi and mss
GSE41258_sample1 <- GSE41258_sample[,-1]
GSE41258_sample1 <- GSE41258_sample1[,!GSE41258_sample1[1,]%in%c(
  "tissue: Polyp", "tissue: Normal Liver", "tissue: Normal Lung", "tissue: Normal Colon",
  "tissue: Microadenoma", "tissue: Polyp, high grade","")]
GSE41258_sammsi <- GSE41258_sample1[3,GSE41258_sample1[2,]%in%"microsattelite instability: MSI-high"]
GSE41258_sammss <- GSE41258_sample1[3,GSE41258_sample1[2,]%in%c(
  "microsattelite instability: MSS","microsattelite instability: MSI-low")]
GSE41258_msi <- as.matrix(GSE41258_matrix2[,colnames(GSE41258_matrix2)%in%GSE41258_sammsi])
GSE41258_mss <- as.matrix(GSE41258_matrix2[,colnames(GSE41258_matrix2)%in%GSE41258_sammss])

# 4. chemotherapy prediction -----------------------------------------------------------------

GSE41258_predmsi <- data.frame(1:ncol(GSE41258_msi))
GSE41258_predmss <- data.frame(1:ncol(GSE41258_mss))

for (i in 1:length(drug)) {
  GSE41258_predmsi[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE41258_msi,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE41258_predmsi)[i] <- drug[i]
}

for (i in 1:length(drug)) {
  GSE41258_predmss[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE41258_mss,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE41258_predmss)[i] <- drug[i]
}

# 5. create table for meta -------------------------------------------------------------

GSE41258_predmsi1 <- as.data.frame(apply(GSE41258_predmsi, 2, function(x)(2^x)))
GSE41258_predmsi1[(nrow(GSE41258_predmsi1)+1),] <- apply(
  GSE41258_predmsi, 2, length_na)
GSE41258_predmsi1[(nrow(GSE41258_predmsi1)+1),] <- apply(
  GSE41258_predmsi1[1:(nrow(GSE41258_predmsi1)-1),], 2, mean, na.rm = T)
GSE41258_predmsi1[(nrow(GSE41258_predmsi1)+1),] <- apply(
  GSE41258_predmsi1[1:(nrow(GSE41258_predmsi1)-2),], 2, var, na.rm = T)
rownames(GSE41258_predmsi1)[(nrow(GSE41258_predmsi1)-2)] <- "n1"
rownames(GSE41258_predmsi1)[(nrow(GSE41258_predmsi1)-1)] <- "mean1"
rownames(GSE41258_predmsi1)[nrow(GSE41258_predmsi1)] <- "sd1"

GSE41258_predmss1 <- as.data.frame(apply(GSE41258_predmss, 2, function(x)(2^x)))
GSE41258_predmss1[(nrow(GSE41258_predmss1)+1),] <- apply(
  GSE41258_predmss, 2, length_na)
GSE41258_predmss1[(nrow(GSE41258_predmss1)+1),] <- apply(
  GSE41258_predmss1[1:(nrow(GSE41258_predmss1)-1),], 2, mean, na.rm = T)
GSE41258_predmss1[(nrow(GSE41258_predmss1)+1),] <- apply(
  GSE41258_predmss1[1:(nrow(GSE41258_predmss1)-2),], 2, var, na.rm = T)
rownames(GSE41258_predmss1)[(nrow(GSE41258_predmss1)-2)] <- "n2"
rownames(GSE41258_predmss1)[(nrow(GSE41258_predmss1)-1)] <- "mean2"
rownames(GSE41258_predmss1)[nrow(GSE41258_predmss1)] <- "sd2"
# GSE92921---------------------------------------------------------------
# 2. read data -----------------------------------------------------------------

GSE92921_geo <- read.table("/mnt/DATA/ytj/GEO_meta/GSE92921_series_matrix.txt.gz",
                           header = T, comment.char = "!", sep = "\t")

GSE92921_sample <- read.xlsx2("/mnt/DATA/ytj/GEO_meta/GSE92921_msi.xlsx", sheetIndex = 1,
                              header = F)

# 3. data preprocessing ----------------------------------------------------------------

# MAS5 standardization requires log2 transformation to perform difference analysis
# GPL570
GSE92921_matrix <- merge(GSE92921_geo, at_ids, by.x = "ID_REF", by.y = "probe_id")

which(duplicated(GSE92921_matrix[,ncol(GSE92921_matrix)]))
# median
GSE92921_matrix1 <- GSE92921_matrix
GSE92921_matrix1[,(ncol(GSE92921_matrix1)+1)] <- apply(
  GSE92921_matrix1[,2:(ncol(GSE92921_matrix1)-1)], 1, median)
GSE92921_matrix1 <- GSE92921_matrix1[order(
  GSE92921_matrix1$symbol, GSE92921_matrix1$V62, decreasing = T),]
GSE92921_matrix1 <- GSE92921_matrix1[which(!duplicated(GSE92921_matrix1$symbol)),]

rownames(GSE92921_matrix1) <- GSE92921_matrix1$symbol
GSE92921_matrix1 <- GSE92921_matrix1[,-c(1, ncol(GSE92921_matrix1), (ncol(GSE92921_matrix1)-1))]
GSE92921_matrix1 <- t(apply(GSE92921_matrix1, 1, log2))
# devided into msi and mss
GSE92921_sample1 <- GSE92921_sample[,-1]
GSE92921_sammsi <- GSE92921_sample1[1,GSE92921_sample1[2,]%in%"msi_status: MSI"]
GSE92921_sammss <- GSE92921_sample1[1,GSE92921_sample1[2,]%in%"msi_status: MSS"]
GSE92921_msi <- as.matrix(GSE92921_matrix1[,colnames(GSE92921_matrix1)%in%GSE92921_sammsi])
GSE92921_mss <- as.matrix(GSE92921_matrix1[,colnames(GSE92921_matrix1)%in%GSE92921_sammss])

# 4. chemotherapy prediction -----------------------------------------------------------------

GSE92921_predmsi <- data.frame(1:ncol(GSE92921_msi))
GSE92921_predmss <- data.frame(1:ncol(GSE92921_mss))

 
for (i in 1:length(drug)) {
  GSE92921_predmsi[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE92921_msi,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE92921_predmsi)[i] <- drug[i]
}

for (i in 1:length(drug)) {
  GSE92921_predmss[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE92921_mss,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE92921_predmss)[i] <- drug[i]
}
# 5. create table for meta -------------------------------------------------------------

GSE92921_predmsi1 <- as.data.frame(apply(GSE92921_predmsi, 2, function(x)(2^x)))
GSE92921_predmsi1[(nrow(GSE92921_predmsi1)+1),] <- apply(
  GSE92921_predmsi, 2, length_na)
GSE92921_predmsi1[(nrow(GSE92921_predmsi1)+1),] <- apply(
  GSE92921_predmsi1[1:(nrow(GSE92921_predmsi1)-1),], 2, mean, na.rm = T)
GSE92921_predmsi1[(nrow(GSE92921_predmsi1)+1),] <- apply(
  GSE92921_predmsi1[1:(nrow(GSE92921_predmsi1)-2),], 2, var, na.rm = T)
rownames(GSE92921_predmsi1)[(nrow(GSE92921_predmsi1)-2)] <- "n1"
rownames(GSE92921_predmsi1)[(nrow(GSE92921_predmsi1)-1)] <- "mean1"
rownames(GSE92921_predmsi1)[nrow(GSE92921_predmsi1)] <- "sd1"

GSE92921_predmss1 <- as.data.frame(apply(GSE92921_predmss, 2, function(x)(2^x)))
GSE92921_predmss1[(nrow(GSE92921_predmss1)+1),] <- apply(
  GSE92921_predmss, 2, length_na)
GSE92921_predmss1[(nrow(GSE92921_predmss1)+1),] <- apply(
  GSE92921_predmss1[1:(nrow(GSE92921_predmss1)-1),], 2, mean, na.rm = T)
GSE92921_predmss1[(nrow(GSE92921_predmss1)+1),] <- apply(
  GSE92921_predmss1[1:(nrow(GSE92921_predmss1)-2),], 2, var, na.rm = T)
rownames(GSE92921_predmss1)[(nrow(GSE92921_predmss1)-2)] <- "n2"
rownames(GSE92921_predmss1)[(nrow(GSE92921_predmss1)-1)] <- "mean2"
rownames(GSE92921_predmss1)[nrow(GSE92921_predmss1)] <- "sd2"
# GSE35566 ----------------------------------------------------------------
# 2. read data -----------------------------------------------------------------

GSE35566_geo <- read.table("/mnt/DATA/ytj/GEO_meta/GSE35566_series_matrix.txt.gz",
                           header = T, comment.char = "!", sep = "\t")

GSE35566_sample <- read.xlsx2("/mnt/DATA/ytj/GEO_meta/GSE35566_msi.xlsx", sheetIndex = 1,
                              header = F)

# 3. data preprocessing ----------------------------------------------------------------

# GPL570
GSE35566_matrix <- merge(GSE35566_geo, at_ids, by.x = "ID_REF", by.y = "probe_id")

which(duplicated(GSE35566_matrix[,ncol(GSE35566_matrix)]))
# median
GSE35566_matrix1 <- GSE35566_matrix
GSE35566_matrix1[,(ncol(GSE35566_matrix1)+1)] <- apply(
  GSE35566_matrix1[,2:(ncol(GSE35566_matrix1)-1)], 1, median)
GSE35566_matrix1 <- GSE35566_matrix1[order(
  GSE35566_matrix1$symbol, GSE35566_matrix1$V22, decreasing = T),]
GSE35566_matrix1 <- GSE35566_matrix1[which(!duplicated(GSE35566_matrix1$symbol)),]

rownames(GSE35566_matrix1) <- GSE35566_matrix1$symbol
GSE35566_matrix1 <- GSE35566_matrix1[,-c(1, ncol(GSE35566_matrix1), (ncol(GSE35566_matrix1)-1))]

# devided into msi and mss
GSE35566_sample1 <- GSE35566_sample[,-1]
GSE35566_sammsi <- GSE35566_sample1[1,GSE35566_sample1[2,]%in%"mss/msi status: MSI"]
GSE35566_sammss <- GSE35566_sample1[1,GSE35566_sample1[2,]%in%"mss/msi status: MSS"]
GSE35566_msi <- as.matrix(GSE35566_matrix1[,colnames(GSE35566_matrix1)%in%GSE35566_sammsi])
GSE35566_mss <- as.matrix(GSE35566_matrix1[,colnames(GSE35566_matrix1)%in%GSE35566_sammss])

# 4. chemotherapy prediction -----------------------------------------------------------------

GSE35566_predmsi <- data.frame(1:ncol(GSE35566_msi))
GSE35566_predmss <- data.frame(1:ncol(GSE35566_mss))

for (i in 1:length(drug)) {
  GSE35566_predmsi[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE35566_msi,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE35566_predmsi)[i] <- drug[i]
}

for (i in 1:length(drug)) {
  GSE35566_predmss[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE35566_mss,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE35566_predmss)[i] <- drug[i]
}


# 5. create table for meta -------------------------------------------------------------

GSE35566_predmsi1 <- as.data.frame(apply(GSE35566_predmsi, 2, function(x)(2^x)))
GSE35566_predmsi1[(nrow(GSE35566_predmsi1)+1),] <- apply(
  GSE35566_predmsi, 2, length_na)
GSE35566_predmsi1[(nrow(GSE35566_predmsi1)+1),] <- apply(
  GSE35566_predmsi1[1:(nrow(GSE35566_predmsi1)-1),], 2, mean, na.rm = T)
GSE35566_predmsi1[(nrow(GSE35566_predmsi1)+1),] <- apply(
  GSE35566_predmsi1[1:(nrow(GSE35566_predmsi1)-2),], 2, var, na.rm = T)
rownames(GSE35566_predmsi1)[(nrow(GSE35566_predmsi1)-2)] <- "n1"
rownames(GSE35566_predmsi1)[(nrow(GSE35566_predmsi1)-1)] <- "mean1"
rownames(GSE35566_predmsi1)[nrow(GSE35566_predmsi1)] <- "sd1"

GSE35566_predmss1 <- as.data.frame(apply(GSE35566_predmss, 2, function(x)(2^x)))
GSE35566_predmss1[(nrow(GSE35566_predmss1)+1),] <- apply(
  GSE35566_predmss, 2, length_na)
GSE35566_predmss1[(nrow(GSE35566_predmss1)+1),] <- apply(
  GSE35566_predmss1[1:(nrow(GSE35566_predmss1)-1),], 2, mean, na.rm = T)
GSE35566_predmss1[(nrow(GSE35566_predmss1)+1),] <- apply(
  GSE35566_predmss1[1:(nrow(GSE35566_predmss1)-2),], 2, var, na.rm = T)
rownames(GSE35566_predmss1)[(nrow(GSE35566_predmss1)-2)] <- "n2"
rownames(GSE35566_predmss1)[(nrow(GSE35566_predmss1)-1)] <- "mean2"
rownames(GSE35566_predmss1)[nrow(GSE35566_predmss1)] <- "sd2"
# GSE29638 ----------------------------------------------------------------
# 2. read data -----------------------------------------------------------------

GSE29638_geo <- read.table("/mnt/DATA/ytj/GEO_meta/GSE29638_series_matrix.txt.gz",
                           sep = "\t", header = T, comment.char = "!")
GSE29638_sample <- read.xlsx2("/mnt/DATA/ytj/GEO_meta/GSE29638_msi.xlsx",
                              sheetIndex = 1)

# 3. data preprocessing ----------------------------------------------------------------

GSE29638_gene <- subset(
  GPL5175, GPL5175$ID%in%GSE29638_geo[,1], select=c(ID,gene_assignment))

GSE29638_gene_anno <- as.data.frame(
  str_split_fixed(GSE29638_gene$gene_assignment, " // ", n=3))
GSE29638_gene_anno2 <- as.data.frame(
  cbind(GSE29638_gene$ID, GSE29638_gene_anno$V2))

colnames(GSE29638_gene_anno2) <- c("ID","symbol")
GSE29638_matrix <- merge(
  GSE29638_geo, GSE29638_gene_anno2, by.x = "ID_REF", by.y = "ID")

which(duplicated(GSE29638_matrix$symbol))
# median
GSE29638_matrix1 <- GSE29638_matrix
GSE29638_matrix1[,(ncol(GSE29638_matrix1)+1)] <- apply(
  GSE29638_matrix1[,2:(ncol(GSE29638_matrix1)-1)], 1, median)

GSE29638_matrix1 <- GSE29638_matrix1[order(
  GSE29638_matrix1$symbol, GSE29638_matrix1$V53, decreasing = T),]
GSE29638_matrix1 <- GSE29638_matrix1[which(!duplicated(GSE29638_matrix1$symbol)),]
rownames(GSE29638_matrix1) <- GSE29638_matrix1[,(ncol(GSE29638_matrix1)-1)]
GSE29638_matrix1 <- GSE29638_matrix1[,-c(1, ncol(GSE29638_matrix1), (ncol(GSE29638_matrix1)-1))]

GSE29638_sample1 <- GSE29638_sample[,-1]
GSE29638_sammsi <- colnames(GSE29638_sample1[,GSE29638_sample1%in%"msi-status: MSI-H"])
GSE29638_sammss <- colnames(GSE29638_sample1[,GSE29638_sample1%in%c(
  "msi-status: MSS","msi-status: MSI-L")])

GSE29638_msi <- as.matrix(GSE29638_matrix1[,colnames(GSE29638_matrix1)%in%GSE29638_sammsi])
GSE29638_mss <- as.matrix(GSE29638_matrix1[,colnames(GSE29638_matrix1)%in%GSE29638_sammss])


# 4. chemotherapy prediction -----------------------------------------------------------------


GSE29638_predmsi <- data.frame(1:ncol(GSE29638_msi))
GSE29638_predmss <- data.frame(1:ncol(GSE29638_mss))

 
for (i in 1:length(drug)) {
  GSE29638_predmsi[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE29638_msi,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE29638_predmsi)[i] <- drug[i]
}

for (i in 1:length(drug)) {
  GSE29638_predmss[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE29638_mss,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE29638_predmss)[i] <- drug[i]
}
# 5. create table for meta -------------------------------------------------------------

GSE29638_predmsi1 <- as.data.frame(apply(GSE29638_predmsi, 2, function(x)(2^x)))
GSE29638_predmsi1[(nrow(GSE29638_predmsi1)+1),] <- apply(
  GSE29638_predmsi, 2, length_na)
GSE29638_predmsi1[(nrow(GSE29638_predmsi1)+1),] <- apply(
  GSE29638_predmsi1[1:(nrow(GSE29638_predmsi1)-1),], 2, mean, na.rm = T)
GSE29638_predmsi1[(nrow(GSE29638_predmsi1)+1),] <- apply(
  GSE29638_predmsi1[1:(nrow(GSE29638_predmsi1)-2),], 2, var, na.rm = T)
rownames(GSE29638_predmsi1)[(nrow(GSE29638_predmsi1)-2)] <- "n1"
rownames(GSE29638_predmsi1)[(nrow(GSE29638_predmsi1)-1)] <- "mean1"
rownames(GSE29638_predmsi1)[nrow(GSE29638_predmsi1)] <- "sd1"

GSE29638_predmss1 <- as.data.frame(apply(GSE29638_predmss, 2, function(x)(2^x)))
GSE29638_predmss1[(nrow(GSE29638_predmss1)+1),] <- apply(
  GSE29638_predmss, 2, length_na)
GSE29638_predmss1[(nrow(GSE29638_predmss1)+1),] <- apply(
  GSE29638_predmss1[1:(nrow(GSE29638_predmss1)-1),], 2, mean, na.rm = T)
GSE29638_predmss1[(nrow(GSE29638_predmss1)+1),] <- apply(
  GSE29638_predmss1[1:(nrow(GSE29638_predmss1)-2),], 2, var, na.rm = T)
rownames(GSE29638_predmss1)[(nrow(GSE29638_predmss1)-2)] <- "n2"
rownames(GSE29638_predmss1)[(nrow(GSE29638_predmss1)-1)] <- "mean2"
rownames(GSE29638_predmss1)[nrow(GSE29638_predmss1)] <- "sd2"
# GSE30378 ----------------------------------------------------------------
# 2. read data -----------------------------------------------------------------

GSE30378_geo <- read.table("/mnt/DATA/ytj/GEO_meta/GSE30378_series_matrix.txt.gz",
                           sep = "\t", header = T, comment.char = "!")

GSE30378_sample <- read.xlsx2("/mnt/DATA/ytj/GEO_meta/GSE30378_msi.xlsx",
                              sheetIndex = 1)

# 3. data preprocessing ----------------------------------------------------------------

GPL5175 <- fread("/mnt/DATA/ytj/GEO_meta/GPL5175-3188.txt")
GSE30378_gene <- subset(
  GPL5175, GPL5175$ID%in%GSE30378_geo[,1], select=c(ID,gene_assignment))

GSE30378_gene_anno <- as.data.frame(
  str_split_fixed(GSE30378_gene$gene_assignment, " // ", n=3))
GSE30378_gene_anno2 <- as.data.frame(
  cbind(GSE30378_gene$ID, GSE30378_gene_anno$V2))

colnames(GSE30378_gene_anno2) <- c("ID","symbol")
GSE30378_matrix <- merge(
  GSE30378_geo, GSE30378_gene_anno2, by.x = "ID_REF", by.y = "ID")

which(duplicated(GSE30378_matrix$symbol))
GSE30378_matrix1 <- GSE30378_matrix
GSE30378_matrix1[,(ncol(GSE30378_matrix1)+1)] <- apply(
  GSE30378_matrix1[,2:(ncol(GSE30378_matrix1)-1)], 1, median)
GSE30378_matrix1 <- GSE30378_matrix1[order(
  GSE30378_matrix1$symbol, GSE30378_matrix1$V98, decreasing = T),]
GSE30378_matrix1 <- GSE30378_matrix1[which(!duplicated(GSE30378_matrix1$symbol)),]
rownames(GSE30378_matrix1) <- GSE30378_matrix1[,(ncol(GSE30378_matrix1)-1)]
GSE30378_matrix1 <- GSE30378_matrix1[,-c(1, ncol(GSE30378_matrix1), (ncol(GSE30378_matrix1)-1))]

GSE30378_sample1 <- GSE30378_sample[,-1]
GSE30378_sammsi <- colnames(GSE30378_sample1[,GSE30378_sample1%in%"msi-status: MSI-H"])
GSE30378_sammss <- colnames(GSE30378_sample1[,GSE30378_sample1%in%c("msi-status: MSI-L","msi-status: MSS")])

GSE30378_msi <- as.matrix(GSE30378_matrix1[,colnames(GSE30378_matrix1)%in%GSE30378_sammsi])
GSE30378_mss <- as.matrix(GSE30378_matrix1[,colnames(GSE30378_matrix1)%in%GSE30378_sammss])


# 4. chemotherapy prediction -----------------------------------------------------------------


GSE30378_predmsi <- data.frame(1:ncol(GSE30378_msi))
GSE30378_predmss <- data.frame(1:ncol(GSE30378_mss))

 
for (i in 1:length(drug)) {
  GSE30378_predmsi[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE30378_msi,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE30378_predmsi)[i] <- drug[i]
}

for (i in 1:length(drug)) {
  GSE30378_predmss[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE30378_mss,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE30378_predmss)[i] <- drug[i]
}
# 5. create table for meta -------------------------------------------------------------

GSE30378_predmsi1 <- as.data.frame(apply(GSE30378_predmsi, 2, function(x)(2^x)))
GSE30378_predmsi1[(nrow(GSE30378_predmsi1)+1),] <- apply(
  GSE30378_predmsi, 2, length_na)
GSE30378_predmsi1[(nrow(GSE30378_predmsi1)+1),] <- apply(
  GSE30378_predmsi1[1:(nrow(GSE30378_predmsi1)-1),], 2, mean, na.rm = T)
GSE30378_predmsi1[(nrow(GSE30378_predmsi1)+1),] <- apply(
  GSE30378_predmsi1[1:(nrow(GSE30378_predmsi1)-2),], 2, var, na.rm = T)
rownames(GSE30378_predmsi1)[(nrow(GSE30378_predmsi1)-2)] <- "n1"
rownames(GSE30378_predmsi1)[(nrow(GSE30378_predmsi1)-1)] <- "mean1"
rownames(GSE30378_predmsi1)[nrow(GSE30378_predmsi1)] <- "sd1"

GSE30378_predmss1 <- as.data.frame(apply(GSE30378_predmss, 2, function(x)(2^x)))
GSE30378_predmss1[(nrow(GSE30378_predmss1)+1),] <- apply(
  GSE30378_predmss, 2, length_na)
GSE30378_predmss1[(nrow(GSE30378_predmss1)+1),] <- apply(
  GSE30378_predmss1[1:(nrow(GSE30378_predmss1)-1),], 2, mean, na.rm = T)
GSE30378_predmss1[(nrow(GSE30378_predmss1)+1),] <- apply(
  GSE30378_predmss1[1:(nrow(GSE30378_predmss1)-2),], 2, var, na.rm = T)
rownames(GSE30378_predmss1)[(nrow(GSE30378_predmss1)-2)] <- "n2"
rownames(GSE30378_predmss1)[(nrow(GSE30378_predmss1)-1)] <- "mean2"
rownames(GSE30378_predmss1)[nrow(GSE30378_predmss1)] <- "sd2"
# GSE25071 ----------------------------------------------------------------
# 2. read data -----------------------------------------------------------------
GPL2986 <- getGEO(filename = "/mnt/DATA/ytj/GEO_meta/GPL2986_family.soft")

GSE25071_geo <- read.table("/mnt/DATA/ytj/GEO_meta/GSE25071_series_matrix.txt.gz",
                           sep = "\t", header = T, comment.char = "!")
GSE25071_sample <- read.xlsx2("/mnt/DATA/ytj/GEO_meta/GSE25071_msi.xlsx",
                              sheetIndex = 1)

# 3. data preprocessing ----------------------------------------------------------------

GPL2986_gene <- GPL2986@dataTable@table[,c(1,6)]

colnames(GPL2986_gene)[1] <- "ID_REF"
GSE25071_matrix <- join(GSE25071_geo, GPL2986_gene, type = "left", match = "first")

which(duplicated(GSE25071_matrix$`Gene Symbol`))
GSE25071_matrix1 <- GSE25071_matrix
GSE25071_matrix1[,(ncol(GSE25071_matrix1)+1)] <- apply(
  GSE25071_matrix1[,2:(ncol(GSE25071_matrix1)-1)], 1, median)
GSE25071_matrix1 <- GSE25071_matrix1[order(
  GSE25071_matrix1$`Gene Symbol`, GSE25071_matrix1$V53, decreasing = T),]
GSE25071_matrix1 <- GSE25071_matrix1[which(!duplicated(GSE25071_matrix1$`Gene Symbol`)),]
GSE25071_matrix1 <- GSE25071_matrix1[-which(is.na(GSE25071_matrix1$`Gene Symbol`)),]
rownames(GSE25071_matrix1) <- GSE25071_matrix1[,(ncol(GSE25071_matrix1)-1)]
GSE25071_matrix1 <- GSE25071_matrix1[,-c(1, ncol(GSE25071_matrix1), (ncol(GSE25071_matrix1)-1))]

GSE25071_sample1 <- GSE25071_sample[,-1]
GSE25071_sammsi <- colnames(GSE25071_sample1[,GSE25071_sample1%in%"microsatellite instability (msi) status: MSI-H"])
GSE25071_sammss <- colnames(GSE25071_sample1[,GSE25071_sample1%in%c(
  "microsatellite instability (msi) status: MSS","microsatellite instability (msi) status: MSI-L")])

GSE25071_msi <- as.matrix(GSE25071_matrix1[,colnames(GSE25071_matrix1)%in%GSE25071_sammsi])
GSE25071_mss <- as.matrix(GSE25071_matrix1[,colnames(GSE25071_matrix1)%in%GSE25071_sammss])


# 4. chemotherapy prediction -----------------------------------------------------------------


GSE25071_predmsi <- data.frame(1:ncol(GSE25071_msi))
GSE25071_predmss <- data.frame(1:ncol(GSE25071_mss))

 
for (i in 1:length(drug)) {
  GSE25071_predmsi[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE25071_msi,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE25071_predmsi)[i] <- drug[i]
}

for (i in 1:length(drug)) {
  GSE25071_predmss[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE25071_mss,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE25071_predmss)[i] <- drug[i]
}
# 5. create table for meta -------------------------------------------------------------

GSE25071_predmsi1 <- as.data.frame(apply(GSE25071_predmsi, 2, function(x)(2^x)))
GSE25071_predmsi1[(nrow(GSE25071_predmsi1)+1),] <- apply(
  GSE25071_predmsi, 2, length_na)
GSE25071_predmsi1[(nrow(GSE25071_predmsi1)+1),] <- apply(
  GSE25071_predmsi1[1:(nrow(GSE25071_predmsi1)-1),], 2, mean, na.rm = T)
GSE25071_predmsi1[(nrow(GSE25071_predmsi1)+1),] <- apply(
  GSE25071_predmsi1[1:(nrow(GSE25071_predmsi1)-2),], 2, var, na.rm = T)
rownames(GSE25071_predmsi1)[(nrow(GSE25071_predmsi1)-2)] <- "n1"
rownames(GSE25071_predmsi1)[(nrow(GSE25071_predmsi1)-1)] <- "mean1"
rownames(GSE25071_predmsi1)[nrow(GSE25071_predmsi1)] <- "sd1"

GSE25071_predmss1 <- as.data.frame(apply(GSE25071_predmss, 2, function(x)(2^x)))
GSE25071_predmss1[(nrow(GSE25071_predmss1)+1),] <- apply(
  GSE25071_predmss, 2, length_na)
GSE25071_predmss1[(nrow(GSE25071_predmss1)+1),] <- apply(
  GSE25071_predmss1[1:(nrow(GSE25071_predmss1)-1),], 2, mean, na.rm = T)
GSE25071_predmss1[(nrow(GSE25071_predmss1)+1),] <- apply(
  GSE25071_predmss1[1:(nrow(GSE25071_predmss1)-2),], 2, var, na.rm = T)
rownames(GSE25071_predmss1)[(nrow(GSE25071_predmss1)-2)] <- "n2"
rownames(GSE25071_predmss1)[(nrow(GSE25071_predmss1)-1)] <- "mean2"
rownames(GSE25071_predmss1)[nrow(GSE25071_predmss1)] <- "sd2"
# GSE24551------------------------------------------
# 2. read data -----------------------------------------------------------------

GSE24551_sample <- read.xlsx2("/mnt/DATA/ytj/GEO_meta/GSE24551-GPL11028_msi.xlsx", sheetIndex = 1)

GSE24551_geo <- read.table("/mnt/DATA/ytj/GEO_meta/GSE24551-GPL11028_series_matrix.txt.gz",
                           header = T, comment.char = "!", sep = "\t")

# 3. data preprocessing ----------------------------------------------------------------

GPL11028 <- read_excel("/mnt/DATA/ytj/GEO_meta/GPL11028.xlsx",sheet = 1)

GSE24551_matrix <- merge(GSE24551_geo , GPL11028, by.x = "ID_REF", by.y = "ID")
which(duplicated(GSE24551_matrix[,ncol(GSE24551_matrix)]))

GSE24551_matrix[,(ncol(GSE24551_matrix)+1)] <- apply(
  GSE24551_matrix[,2:(ncol(GSE24551_matrix)-1)], 1, median)
GSE24551_matrix <- GSE24551_matrix[order(
  GSE24551_matrix$gene_symbol, GSE24551_matrix$V163, decreasing = T),]
GSE24551_matrix <- GSE24551_matrix[which(!duplicated(GSE24551_matrix$gene_symbol)),]
# a missing gene name, directly delete that row
GSE24551_matrix <- GSE24551_matrix[-which(is.na(GSE24551_matrix$gene_symbol)),]
rownames(GSE24551_matrix) <- GSE24551_matrix$gene_symbol
GSE24551_matrix <- GSE24551_matrix[,-c(
  1,ncol(GSE24551_matrix),(ncol(GSE24551_matrix)-1))]
# Due to missing values in the expression matrix, the drug sensitivity prediction analysis reported an error
# Error in solve.default(crossprod(des), crossprod(des, y1)) : 
#   Lapack routine dgesv: system is exactly singular: U[1,1] = 0
GSE24551_matrix[which(is.na(GSE24551_matrix)),]
# Those genes are missing values in all samples, so we directly deleted those gene
GSE24551_matrix1 <- GSE24551_matrix[complete.cases(GSE24551_matrix),]
which(is.na(GSE24551_matrix1))


GSE24551_sample1 <- GSE24551_sample[,-1]
GSE24551_sammsi <- colnames(GSE24551_sample1[,GSE24551_sample1[1,]%in%"msi-status: MSI-H"])
GSE24551_sammss <- colnames(GSE24551_sample1[,GSE24551_sample1[1,]%in%c(
  "msi-status: MSI-L","msi-status: MSS")])
GSE24551_msi <- as.matrix(GSE24551_matrix1[,colnames(GSE24551_matrix1)%in%GSE24551_sammsi])
GSE24551_mss <- as.matrix(GSE24551_matrix1[,colnames(GSE24551_matrix1)%in%GSE24551_sammss])


# 4. chemotherapy prediction -----------------------------------------------------------------


GSE24551_predmsi <- data.frame(1:ncol(GSE24551_msi))
GSE24551_predmss <- data.frame(1:ncol(GSE24551_mss))

for (i in 1:length(drug)) {
  GSE24551_predmsi[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE24551_msi,
    drug = drug[1],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE24551_predmsi)[i] <- drug[i]
}

for (i in 1:length(drug)) {
  GSE24551_predmss[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE24551_mss,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE24551_predmss)[i] <- drug[i]
}

# 5. create table for meta -------------------------------------------------------------

GSE24551_predmsi1 <- as.data.frame(apply(GSE24551_predmsi, 2, function(x)(2^x)))
GSE24551_predmsi1[(nrow(GSE24551_predmsi1)+1),] <- apply(
  GSE24551_predmsi, 2, length_na)
GSE24551_predmsi1[(nrow(GSE24551_predmsi1)+1),] <- apply(
  GSE24551_predmsi1[1:(nrow(GSE24551_predmsi1)-1),], 2, mean, na.rm = T)
GSE24551_predmsi1[(nrow(GSE24551_predmsi1)+1),] <- apply(
  GSE24551_predmsi1[1:(nrow(GSE24551_predmsi1)-2),], 2, var, na.rm = T)
rownames(GSE24551_predmsi1)[(nrow(GSE24551_predmsi1)-2)] <- "n1"
rownames(GSE24551_predmsi1)[(nrow(GSE24551_predmsi1)-1)] <- "mean1"
rownames(GSE24551_predmsi1)[nrow(GSE24551_predmsi1)] <- "sd1"

GSE24551_predmss1 <- as.data.frame(apply(GSE24551_predmss, 2, function(x)(2^x)))
GSE24551_predmss1[(nrow(GSE24551_predmss1)+1),] <- apply(
  GSE24551_predmss, 2, length_na)
GSE24551_predmss1[(nrow(GSE24551_predmss1)+1),] <- apply(
  GSE24551_predmss1[1:(nrow(GSE24551_predmss1)-1),], 2, mean, na.rm = T)
GSE24551_predmss1[(nrow(GSE24551_predmss1)+1),] <- apply(
  GSE24551_predmss1[1:(nrow(GSE24551_predmss1)-2),], 2, var, na.rm = T)
rownames(GSE24551_predmss1)[(nrow(GSE24551_predmss1)-2)] <- "n2"
rownames(GSE24551_predmss1)[(nrow(GSE24551_predmss1)-1)] <- "mean2"
rownames(GSE24551_predmss1)[nrow(GSE24551_predmss1)] <- "sd2"
# GSE24550------------------------------------------
# 2. read data -----------------------------------------------------------------

GSE24550_sample <- read.xlsx2("/mnt/DATA/ytj/GEO_meta/GSE24550-GPL11028_msi.xlsx", sheetIndex = 1)

GSE24550_geo <- read.table("/mnt/DATA/ytj/GEO_meta/GSE24550-GPL11028_series_matrix.txt.gz",
                           header = T, comment.char = "!", sep = "\t")

# 3. data preprocessing ----------------------------------------------------------------

GSE24550_matrix <- merge(GSE24550_geo , GPL11028, by.x = "ID_REF", by.y = "ID")
which(duplicated(GSE24550_matrix[,ncol(GSE24550_matrix)]))

GSE24550_matrix[,(ncol(GSE24550_matrix)+1)] <- apply(
  GSE24550_matrix[,2:(ncol(GSE24550_matrix)-1)], 1, median)
GSE24550_matrix <- GSE24550_matrix[order(
  GSE24550_matrix$gene_symbol, GSE24550_matrix$V80, decreasing = T),]
GSE24550_matrix <- GSE24550_matrix[which(!duplicated(GSE24550_matrix$gene_symbol)),]
# A gene name is missing, directly delete that line
GSE24550_matrix <- GSE24550_matrix[-which(is.na(GSE24550_matrix$gene_symbol)),]
rownames(GSE24550_matrix) <- GSE24550_matrix$gene_symbol
GSE24550_matrix1 <- GSE24550_matrix[,-c(
  1,ncol(GSE24550_matrix),(ncol(GSE24550_matrix)-1))]
GSE24550_matrix1[which(is.na(GSE24550_matrix1)),]
GSE24550_matrix1 <- GSE24550_matrix1[which(!is.na(GSE24550_matrix1)),]

GSE24550_sample1 <- GSE24550_sample[,-1]
GSE24550_sammsi <- colnames(GSE24550_sample1[,GSE24550_sample1[1,]%in%"msi-status: MSI-H"])
GSE24550_sammss <- colnames(GSE24550_sample1[,GSE24550_sample1[1,]%in%c(
  "msi-status: MSI-L","msi-status: MSS")])
GSE24550_msi <- as.matrix(GSE24550_matrix1[,colnames(GSE24550_matrix1)%in%GSE24550_sammsi])
GSE24550_mss <- as.matrix(GSE24550_matrix1[,colnames(GSE24550_matrix1)%in%GSE24550_sammss])



# 4. chemotherapy prediction -----------------------------------------------------------------


GSE24550_predmsi <- data.frame(1:ncol(GSE24550_msi))
GSE24550_predmss <- data.frame(1:ncol(GSE24550_mss))

for (i in 1:length(drug)) {
  GSE24550_predmsi[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE24550_msi,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE24550_predmsi)[i] <- drug[i]
}

for (i in 1:length(drug)) {
  GSE24550_predmss[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE24550_mss,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE24550_predmss)[i] <- drug[i]
}

# 5. create table for meta -------------------------------------------------------------

GSE24550_predmsi1 <- as.data.frame(apply(GSE24550_predmsi, 2, function(x)(2^x)))
GSE24550_predmsi1[(nrow(GSE24550_predmsi1)+1),] <- apply(
  GSE24550_predmsi, 2, length_na)
GSE24550_predmsi1[(nrow(GSE24550_predmsi1)+1),] <- apply(
  GSE24550_predmsi1[1:(nrow(GSE24550_predmsi1)-1),], 2, mean, na.rm = T)
GSE24550_predmsi1[(nrow(GSE24550_predmsi1)+1),] <- apply(
  GSE24550_predmsi1[1:(nrow(GSE24550_predmsi1)-2),], 2, var, na.rm = T)
rownames(GSE24550_predmsi1)[(nrow(GSE24550_predmsi1)-2)] <- "n1"
rownames(GSE24550_predmsi1)[(nrow(GSE24550_predmsi1)-1)] <- "mean1"
rownames(GSE24550_predmsi1)[nrow(GSE24550_predmsi1)] <- "sd1"

GSE24550_predmss1 <- as.data.frame(apply(GSE24550_predmss, 2, function(x)(2^x)))
GSE24550_predmss1[(nrow(GSE24550_predmss1)+1),] <- apply(
  GSE24550_predmss, 2, length_na)
GSE24550_predmss1[(nrow(GSE24550_predmss1)+1),] <- apply(
  GSE24550_predmss1[1:(nrow(GSE24550_predmss1)-1),], 2, mean, na.rm = T)
GSE24550_predmss1[(nrow(GSE24550_predmss1)+1),] <- apply(
  GSE24550_predmss1[1:(nrow(GSE24550_predmss1)-2),], 2, var, na.rm = T)
rownames(GSE24550_predmss1)[(nrow(GSE24550_predmss1)-2)] <- "n2"
rownames(GSE24550_predmss1)[(nrow(GSE24550_predmss1)-1)] <- "mean2"
rownames(GSE24550_predmss1)[nrow(GSE24550_predmss1)] <- "sd2"
# GSE18088 ----------------------------------------------------------------
# 2. read data -----------------------------------------------------------------

GSE18088_geo <- read.table("/mnt/DATA/ytj/GEO_meta/GSE18088_series_matrix.txt.gz",
                           header = T, comment.char = "!", sep = "\t")

GSE18088_sample <- read.xlsx2("/mnt/DATA/ytj/GEO_meta/GSE18088_msi.xlsx", sheetIndex = 1,
                              header = F)

# 3. data preprocessing ----------------------------------------------------------------

# VSN standardization should allow for direct differential analysis, but the probe ID needs to be converted to GPL570
GSE18088_matrix <- merge(GSE18088_geo, at_ids, by.x = "ID_REF", by.y = "probe_id")

which(duplicated(GSE18088_matrix[,ncol(GSE18088_matrix)]))

GSE18088_matrix1 <- GSE18088_matrix
GSE18088_matrix1[,(ncol(GSE18088_matrix1)+1)] <- apply(
  GSE18088_matrix1[,2:(ncol(GSE18088_matrix1)-1)], 1, median)
GSE18088_matrix1 <- GSE18088_matrix1[order(
  GSE18088_matrix1$symbol, GSE18088_matrix1$V56, decreasing = T),]
GSE18088_matrix1 <- GSE18088_matrix1[which(!duplicated(GSE18088_matrix1$symbol)),]

rownames(GSE18088_matrix1) <- GSE18088_matrix1$symbol
GSE18088_matrix1 <- GSE18088_matrix1[,-c(1, ncol(GSE18088_matrix1), (ncol(GSE18088_matrix1)-1))]

GSE18088_sample1 <- GSE18088_sample[,-1]
GSE18088_sammsi <- GSE18088_sample1[1,GSE18088_sample1[2,]%in%"microsatellite status: MSI-high"]
GSE18088_sammss <- GSE18088_sample1[1,GSE18088_sample1[2,]%in%"microsatellite status: MSS"]
GSE18088_msi <- as.matrix(GSE18088_matrix1[,colnames(GSE18088_matrix1)%in%GSE18088_sammsi])
GSE18088_mss <- as.matrix(GSE18088_matrix1[,colnames(GSE18088_matrix1)%in%GSE18088_sammss])

# 4. chemotherapy prediction -----------------------------------------------------------------

GSE18088_predmsi <- data.frame(1:ncol(GSE18088_msi))
GSE18088_predmss <- data.frame(1:ncol(GSE18088_mss))

for (i in 1:length(drug)) {
  GSE18088_predmsi[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE18088_msi,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE18088_predmsi)[i] <- drug[i]
}

for (i in 1:length(drug)) {
  GSE18088_predmss[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE18088_mss,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE18088_predmss)[i] <- drug[i]
}

# 5. create table for meta -------------------------------------------------------------

GSE18088_predmsi1 <- as.data.frame(apply(GSE18088_predmsi, 2, function(x)(2^x)))
GSE18088_predmsi1[(nrow(GSE18088_predmsi1)+1),] <- apply(
  GSE18088_predmsi, 2, length_na)
GSE18088_predmsi1[(nrow(GSE18088_predmsi1)+1),] <- apply(
  GSE18088_predmsi1[1:(nrow(GSE18088_predmsi1)-1),], 2, mean, na.rm = T)
GSE18088_predmsi1[(nrow(GSE18088_predmsi1)+1),] <- apply(
  GSE18088_predmsi1[1:(nrow(GSE18088_predmsi1)-2),], 2, var, na.rm = T)
rownames(GSE18088_predmsi1)[(nrow(GSE18088_predmsi1)-2)] <- "n1"
rownames(GSE18088_predmsi1)[(nrow(GSE18088_predmsi1)-1)] <- "mean1"
rownames(GSE18088_predmsi1)[nrow(GSE18088_predmsi1)] <- "sd1"

GSE18088_predmss1 <- as.data.frame(apply(GSE18088_predmss, 2, function(x)(2^x)))
GSE18088_predmss1[(nrow(GSE18088_predmss1)+1),] <- apply(
  GSE18088_predmss, 2, length_na)
GSE18088_predmss1[(nrow(GSE18088_predmss1)+1),] <- apply(
  GSE18088_predmss1[1:(nrow(GSE18088_predmss1)-1),], 2, mean, na.rm = T)
GSE18088_predmss1[(nrow(GSE18088_predmss1)+1),] <- apply(
  GSE18088_predmss1[1:(nrow(GSE18088_predmss1)-2),], 2, var, na.rm = T)
rownames(GSE18088_predmss1)[(nrow(GSE18088_predmss1)-2)] <- "n2"
rownames(GSE18088_predmss1)[(nrow(GSE18088_predmss1)-1)] <- "mean2"
rownames(GSE18088_predmss1)[nrow(GSE18088_predmss1)] <- "sd2"
# GSE26682 ----------------------------------------------------------------
# 2. read data -----------------------------------------------------------------

GSE26682_geo <- read.table("/mnt/DATA/ytj/GEO_meta/GSE26682-GPL570_series_matrix.txt.gz",
                           header = T, comment.char = "!", sep = "\t")

GSE26682_sample <- read.xlsx2("/mnt/DATA/ytj/GEO_meta/GSE26682-GPL570_msi.xlsx", sheetIndex = 1,
                              header = F)

# 3. data preprocessing ----------------------------------------------------------------

GSE26682_matrix <- merge(GSE26682_geo, at_ids, by.x = "ID_REF", by.y = "probe_id")

which(duplicated(GSE26682_matrix[,ncol(GSE26682_matrix)]))

GSE26682_matrix1 <- GSE26682_matrix
GSE26682_matrix1[,(ncol(GSE26682_matrix1)+1)] <- apply(
  GSE26682_matrix1[,2:(ncol(GSE26682_matrix1)-1)], 1, median)
GSE26682_matrix1 <- GSE26682_matrix1[order(
  GSE26682_matrix1$symbol, GSE26682_matrix1$V179, decreasing = T),]
GSE26682_matrix1 <- GSE26682_matrix1[which(!duplicated(GSE26682_matrix1$symbol)),]

rownames(GSE26682_matrix1) <- GSE26682_matrix1$symbol
GSE26682_matrix1 <- GSE26682_matrix1[,-c(1, ncol(GSE26682_matrix1), (ncol(GSE26682_matrix1)-1))]

GSE26682_sample1 <- GSE26682_sample[,-1]
GSE26682_sammsi <- GSE26682_sample1[1,GSE26682_sample1[2,]%in%"microsatellite instability (msi) status: High [MSI-H]"]
GSE26682_sammss <- GSE26682_sample1[1,GSE26682_sample1[2,]%in%c(
  "microsatellite instability (msi) status: Stable [MSS]","microsatellite instability (msi) status: Low [MSI-L]")]
GSE26682_msi <- as.matrix(GSE26682_matrix1[,colnames(GSE26682_matrix1)%in%GSE26682_sammsi])
GSE26682_mss <- as.matrix(GSE26682_matrix1[,colnames(GSE26682_matrix1)%in%GSE26682_sammss])

# 4. chemotherapy prediction -----------------------------------------------------------------

GSE26682_predmsi <- data.frame(1:ncol(GSE26682_msi))
GSE26682_predmss <- data.frame(1:ncol(GSE26682_mss))

for (i in 1:length(drug)) {
  GSE26682_predmsi[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE26682_msi,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE26682_predmsi)[i] <- drug[i]
}

for (i in 1:length(drug)) {
  GSE26682_predmss[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE26682_mss,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE26682_predmss)[i] <- drug[i]
}

# 5. create table for meta -------------------------------------------------------------

GSE26682_predmsi1 <- as.data.frame(apply(GSE26682_predmsi, 2, function(x)(2^x)))
GSE26682_predmsi1[(nrow(GSE26682_predmsi1)+1),] <- apply(
  GSE26682_predmsi, 2, length_na)
GSE26682_predmsi1[(nrow(GSE26682_predmsi1)+1),] <- apply(
  GSE26682_predmsi1[1:(nrow(GSE26682_predmsi1)-1),], 2, mean, na.rm = T)
GSE26682_predmsi1[(nrow(GSE26682_predmsi1)+1),] <- apply(
  GSE26682_predmsi1[1:(nrow(GSE26682_predmsi1)-2),], 2, var, na.rm = T)
rownames(GSE26682_predmsi1)[(nrow(GSE26682_predmsi1)-2)] <- "n1"
rownames(GSE26682_predmsi1)[(nrow(GSE26682_predmsi1)-1)] <- "mean1"
rownames(GSE26682_predmsi1)[nrow(GSE26682_predmsi1)] <- "sd1"

GSE26682_predmss1 <- as.data.frame(apply(GSE26682_predmss, 2, function(x)(2^x)))
GSE26682_predmss1[(nrow(GSE26682_predmss1)+1),] <- apply(
  GSE26682_predmss, 2, length_na)
GSE26682_predmss1[(nrow(GSE26682_predmss1)+1),] <- apply(
  GSE26682_predmss1[1:(nrow(GSE26682_predmss1)-1),], 2, mean, na.rm = T)
GSE26682_predmss1[(nrow(GSE26682_predmss1)+1),] <- apply(
  GSE26682_predmss1[1:(nrow(GSE26682_predmss1)-2),], 2, var, na.rm = T)
rownames(GSE26682_predmss1)[(nrow(GSE26682_predmss1)-2)] <- "n2"
rownames(GSE26682_predmss1)[(nrow(GSE26682_predmss1)-1)] <- "mean2"
rownames(GSE26682_predmss1)[nrow(GSE26682_predmss1)] <- "sd2"
# GSE27544 ----------------------------------------------------------------
# 2. read data -----------------------------------------------------------------

GSE27544_sample <- read.xlsx2("/mnt/DATA/ytj/GEO_meta/GSE27544_msi.xlsx", sheetIndex = 1)

GSE27544_geo <- read.table("/mnt/DATA/ytj/GEO_meta/GSE27544_series_matrix.txt.gz",
                           header = T, comment.char = "!")

# 3. data preprocessing ----------------------------------------------------------------

GSE27544_matrix <- merge(GSE27544_geo , GPL13158, by.x = "ID_REF", by.y = "ID")
which(duplicated(GSE27544_matrix[,ncol(GSE27544_matrix)]))
GSE27544_matrix[,(ncol(GSE27544_matrix)+1)] <- apply(
  GSE27544_matrix[,2:(ncol(GSE27544_matrix)-1)], 1, median)
GSE27544_matrix <- GSE27544_matrix[order(
  GSE27544_matrix$Gene.Symbol, GSE27544_matrix$V25, decreasing = T),]
GSE27544_matrix <- GSE27544_matrix[which(!duplicated(GSE27544_matrix$Gene.Symbol)),]
rownames(GSE27544_matrix) <- GSE27544_matrix$Gene.Symbol
GSE27544_matrix1 <- GSE27544_matrix[,-c(
  1,ncol(GSE27544_matrix),(ncol(GSE27544_matrix)-1))]

GSE27544_sample1 <- GSE27544_sample[,-1]
GSE27544_sammsi <- colnames(
  GSE27544_sample1[,GSE27544_sample1[3,]%in%"microsatellite phenotype: high-level microsatellite instability"])
GSE27544_sammss <- colnames(
  GSE27544_sample1[,GSE27544_sample1[3,]%in%"microsatellite phenotype: microsatellite stability"])
GSE27544_msi <- as.matrix(GSE27544_matrix1[,colnames(GSE27544_matrix1)%in%GSE27544_sammsi])
GSE27544_mss <- as.matrix(GSE27544_matrix1[,colnames(GSE27544_matrix1)%in%GSE27544_sammss])


# 4. chemotherapy prediction -----------------------------------------------------------------


GSE27544_predmsi <- data.frame(1:ncol(GSE27544_msi))
GSE27544_predmss <- data.frame(1:ncol(GSE27544_mss))

for (i in 1:length(drug)) {
  GSE27544_predmsi[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE27544_msi,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE27544_predmsi)[i] <- drug[i]
}

for (i in 1:length(drug)) {
  GSE27544_predmss[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE27544_mss,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE27544_predmss)[i] <- drug[i]
}


# 5. create table for meta -------------------------------------------------------------

GSE27544_predmsi1 <- as.data.frame(apply(GSE27544_predmsi, 2, function(x)(2^x)))
GSE27544_predmsi1[(nrow(GSE27544_predmsi1)+1),] <- apply(
  GSE27544_predmsi, 2, length_na)
GSE27544_predmsi1[(nrow(GSE27544_predmsi1)+1),] <- apply(
  GSE27544_predmsi1[1:(nrow(GSE27544_predmsi1)-1),], 2, mean, na.rm = T)
GSE27544_predmsi1[(nrow(GSE27544_predmsi1)+1),] <- apply(
  GSE27544_predmsi1[1:(nrow(GSE27544_predmsi1)-2),], 2, var, na.rm = T)
rownames(GSE27544_predmsi1)[(nrow(GSE27544_predmsi1)-2)] <- "n1"
rownames(GSE27544_predmsi1)[(nrow(GSE27544_predmsi1)-1)] <- "mean1"
rownames(GSE27544_predmsi1)[nrow(GSE27544_predmsi1)] <- "sd1"

GSE27544_predmss1 <- as.data.frame(apply(GSE27544_predmss, 2, function(x)(2^x)))
GSE27544_predmss1[(nrow(GSE27544_predmss1)+1),] <- apply(
  GSE27544_predmss, 2, length_na)
GSE27544_predmss1[(nrow(GSE27544_predmss1)+1),] <- apply(
  GSE27544_predmss1[1:(nrow(GSE27544_predmss1)-1),], 2, mean, na.rm = T)
GSE27544_predmss1[(nrow(GSE27544_predmss1)+1),] <- apply(
  GSE27544_predmss1[1:(nrow(GSE27544_predmss1)-2),], 2, var, na.rm = T)
rownames(GSE27544_predmss1)[(nrow(GSE27544_predmss1)-2)] <- "n2"
rownames(GSE27544_predmss1)[(nrow(GSE27544_predmss1)-1)] <- "mean2"
rownames(GSE27544_predmss1)[nrow(GSE27544_predmss1)] <- "sd2"
# GSE24795 ----------------------------------------------------------------
# 2. read data -----------------------------------------------------------------

GSE24795_geo <- read.table("/mnt/DATA/ytj/GEO_meta/GSE24795_series_matrix.txt.gz",
                           header = T, comment.char = "!", sep = "\t")

GSE24795_sample <- read.xlsx2("/mnt/DATA/ytj/GEO_meta/GSE24795_msi.xlsx", sheetIndex = 1,
                              header = F)

# 3. data preprocessing ----------------------------------------------------------------

GSE24795_matrix <- merge(GSE24795_geo, at_ids, by.x = "ID_REF", by.y = "probe_id")

which(duplicated(GSE24795_matrix[,ncol(GSE24795_matrix)]))

GSE24795_matrix1 <- GSE24795_matrix
GSE24795_matrix1[,(ncol(GSE24795_matrix1)+1)] <- apply(
  GSE24795_matrix1[,2:(ncol(GSE24795_matrix1)-1)], 1, median)
GSE24795_matrix1 <- GSE24795_matrix1[order(
  GSE24795_matrix1$symbol, GSE24795_matrix1$V33, decreasing = T),]
GSE24795_matrix1 <- GSE24795_matrix1[which(!duplicated(GSE24795_matrix1$symbol)),]

rownames(GSE24795_matrix1) <- GSE24795_matrix1$symbol
GSE24795_matrix1 <- GSE24795_matrix1[,-c(1, ncol(GSE24795_matrix1), (ncol(GSE24795_matrix1)-1))]

GSE24795_sample1 <- GSE24795_sample[,-1]
GSE24795_sammsi <- GSE24795_sample1[1,GSE24795_sample1[2,]%in%"genotype/variation: Replication error positive (RER+/MSI+)"]
GSE24795_sammss <- GSE24795_sample1[1,GSE24795_sample1[2,]%in%"genotype/variation: Replication error negative (RER-/MSI-)"]
GSE24795_msi <- as.matrix(GSE24795_matrix1[,colnames(GSE24795_matrix1)%in%GSE24795_sammsi])
GSE24795_mss <- as.matrix(GSE24795_matrix1[,colnames(GSE24795_matrix1)%in%GSE24795_sammss])

# 4. chemotherapy prediction -----------------------------------------------------------------

GSE24795_predmsi <- data.frame(1:ncol(GSE24795_msi))
GSE24795_predmss <- data.frame(1:ncol(GSE24795_mss))

for (i in 1:length(drug)) {
  GSE24795_predmsi[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE24795_msi,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE24795_predmsi)[i] <- drug[i]
}

for (i in 1:length(drug)) {
  GSE24795_predmss[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE24795_mss,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE24795_predmss)[i] <- drug[i]
}


# 5. create table for meta -------------------------------------------------------------

GSE24795_predmsi1 <- as.data.frame(apply(GSE24795_predmsi, 2, function(x)(2^x)))
GSE24795_predmsi1[(nrow(GSE24795_predmsi1)+1),] <- apply(
  GSE24795_predmsi, 2, length_na)
GSE24795_predmsi1[(nrow(GSE24795_predmsi1)+1),] <- apply(
  GSE24795_predmsi1[1:(nrow(GSE24795_predmsi1)-1),], 2, mean, na.rm = T)
GSE24795_predmsi1[(nrow(GSE24795_predmsi1)+1),] <- apply(
  GSE24795_predmsi1[1:(nrow(GSE24795_predmsi1)-2),], 2, var, na.rm = T)
rownames(GSE24795_predmsi1)[(nrow(GSE24795_predmsi1)-2)] <- "n1"
rownames(GSE24795_predmsi1)[(nrow(GSE24795_predmsi1)-1)] <- "mean1"
rownames(GSE24795_predmsi1)[nrow(GSE24795_predmsi1)] <- "sd1"

GSE24795_predmss1 <- as.data.frame(apply(GSE24795_predmss, 2, function(x)(2^x)))
GSE24795_predmss1[(nrow(GSE24795_predmss1)+1),] <- apply(
  GSE24795_predmss, 2, length_na)
GSE24795_predmss1[(nrow(GSE24795_predmss1)+1),] <- apply(
  GSE24795_predmss1[1:(nrow(GSE24795_predmss1)-1),], 2, mean, na.rm = T)
GSE24795_predmss1[(nrow(GSE24795_predmss1)+1),] <- apply(
  GSE24795_predmss1[1:(nrow(GSE24795_predmss1)-2),], 2, var, na.rm = T)
rownames(GSE24795_predmss1)[(nrow(GSE24795_predmss1)-2)] <- "n2"
rownames(GSE24795_predmss1)[(nrow(GSE24795_predmss1)-1)] <- "mean2"
rownames(GSE24795_predmss1)[nrow(GSE24795_predmss1)] <- "sd2"
# GSE4459 ---------------------------------------------------
# 2. read data -----------------------------------------------------------------

GSE4459_geo <- read.table("/mnt/DATA/ytj/GEO_meta/GSE4459_series_matrix.txt.gz",
                          header = T, comment.char = "!", sep = "\t")

GSE4459_sample <- read.xlsx2("/mnt/DATA/ytj/GEO_meta/GSE4459_msi.xlsx", sheetIndex = 1,
                             header = F)
# 3. data preprocessing ----------------------------------------------------------------

GSE4459_matrix <- merge(GSE4459_geo, at_ids, by.x = "ID_REF", by.y = "probe_id")

which(duplicated(GSE4459_matrix[,ncol(GSE4459_matrix)]))

GSE4459_matrix1 <- GSE4459_matrix
GSE4459_matrix1[,(ncol(GSE4459_matrix1)+1)] <- apply(
  GSE4459_matrix1[,2:(ncol(GSE4459_matrix1)-1)], 1, median)
GSE4459_matrix1 <- GSE4459_matrix1[order(
  GSE4459_matrix1$symbol, GSE4459_matrix1$V92, decreasing = T),]
GSE4459_matrix1 <- GSE4459_matrix1[which(!duplicated(GSE4459_matrix1$symbol)),]

rownames(GSE4459_matrix1) <- GSE4459_matrix1$symbol
GSE4459_matrix1 <- GSE4459_matrix1[,-c(1, ncol(GSE4459_matrix1), (ncol(GSE4459_matrix1)-1))]
GSE4459_matrix1 <- t(apply(GSE4459_matrix1+1, 1, log2))
# Because there is a 0 in this expression spectrum, log2 becomes negative infinity, and the drug sensitivity analysis later reported an error
# Error: BiocParallel errors
# element index: 1, 2
# first error: missing value where TRUE/FALSE needed


GSE4459_sample1 <- GSE4459_sample[,-1]
GSE4459_sammsi <- GSE4459_sample1[1,GSE4459_sample1[2,]%in%c(
  "MSI, proximal","MSI, distal","MSI, unknown")]
GSE4459_sammss <- GSE4459_sample1[1,GSE4459_sample1[2,]%in%c(
  "MSS, distal","MSS, proximal")]
GSE4459_msi <- as.matrix(GSE4459_matrix1[,colnames(GSE4459_matrix1)%in%GSE4459_sammsi])
GSE4459_mss <- as.matrix(GSE4459_matrix1[,colnames(GSE4459_matrix1)%in%GSE4459_sammss])

# 4. chemotherapy prediction -----------------------------------------------------------------

GSE4459_predmsi <- data.frame(1:ncol(GSE4459_msi))
GSE4459_predmss <- data.frame(1:ncol(GSE4459_mss))

 
for (i in 1:length(drug)) {
  GSE4459_predmsi[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE4459_msi,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE4459_predmsi)[i] <- drug[i]
}

for (i in 1:length(drug)) {
  GSE4459_predmss[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE4459_mss,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE4459_predmss)[i] <- drug[i]
}

# 5. create table for meta -------------------------------------------------------------

GSE4459_predmsi1 <- as.data.frame(apply(GSE4459_predmsi, 2, function(x)(2^x)))
GSE4459_predmsi1[(nrow(GSE4459_predmsi1)+1),] <- apply(
  GSE4459_predmsi, 2, length_na)
GSE4459_predmsi1[(nrow(GSE4459_predmsi1)+1),] <- apply(
  GSE4459_predmsi1[1:(nrow(GSE4459_predmsi1)-1),], 2, mean, na.rm = T)
GSE4459_predmsi1[(nrow(GSE4459_predmsi1)+1),] <- apply(
  GSE4459_predmsi1[1:(nrow(GSE4459_predmsi1)-2),], 2, var, na.rm = T)
rownames(GSE4459_predmsi1)[(nrow(GSE4459_predmsi1)-2)] <- "n1"
rownames(GSE4459_predmsi1)[(nrow(GSE4459_predmsi1)-1)] <- "mean1"
rownames(GSE4459_predmsi1)[nrow(GSE4459_predmsi1)] <- "sd1"

GSE4459_predmss1 <- as.data.frame(apply(GSE4459_predmss, 2, function(x)(2^x)))
GSE4459_predmss1[(nrow(GSE4459_predmss1)+1),] <- apply(
  GSE4459_predmss, 2, length_na)
GSE4459_predmss1[(nrow(GSE4459_predmss1)+1),] <- apply(
  GSE4459_predmss1[1:(nrow(GSE4459_predmss1)-1),], 2, mean, na.rm = T)
GSE4459_predmss1[(nrow(GSE4459_predmss1)+1),] <- apply(
  GSE4459_predmss1[1:(nrow(GSE4459_predmss1)-2),], 2, var, na.rm = T)
rownames(GSE4459_predmss1)[(nrow(GSE4459_predmss1)-2)] <- "n2"
rownames(GSE4459_predmss1)[(nrow(GSE4459_predmss1)-1)] <- "mean2"
rownames(GSE4459_predmss1)[nrow(GSE4459_predmss1)] <- "sd2"
# GSE11543 ----------------------------------------------------------------
# 2. read data -----------------------------------------------------------------

GSE11543_geo <- read.table("/mnt/DATA/ytj/GEO_meta/GSE11543_series_matrix.txt.gz",
                           header = T, comment.char = "!", sep = "\t")

GSE11543_sample <- read.xlsx2("/mnt/DATA/ytj/GEO_meta/GSE11543_msi.xlsx", sheetIndex = 1,
                              header = F)
# 3. data preprocessing ----------------------------------------------------------------

at_ids_hu6800 <- toTable(hu6800SYMBOL)
GSE11543_matrix <- merge(GSE11543_geo, at_ids_hu6800, by.x = "ID_REF", by.y = "probe_id")

which(duplicated(GSE11543_matrix[,ncol(GSE11543_matrix)]))

GSE11543_matrix1 <- GSE11543_matrix
GSE11543_matrix1[,(ncol(GSE11543_matrix1)+1)] <- apply(
  GSE11543_matrix1[,2:(ncol(GSE11543_matrix1)-1)], 1, median)
GSE11543_matrix1 <- GSE11543_matrix1[order(
  GSE11543_matrix1$symbol, GSE11543_matrix1$V54, decreasing = T),]
GSE11543_matrix1 <- GSE11543_matrix1[which(!duplicated(GSE11543_matrix1$symbol)),]

rownames(GSE11543_matrix1) <- GSE11543_matrix1$symbol
GSE11543_matrix1 <- GSE11543_matrix1[,-c(1, ncol(GSE11543_matrix1), (ncol(GSE11543_matrix1)-1))]
GSE11543_matrix1 <- t(apply(GSE11543_matrix1, 1, log2))

GSE11543_sample1 <- GSE11543_sample[,-1]
GSE11543_sample1[1,] <- as.data.frame(
  strsplit(as.character(GSE11543_sample1[1,]), "colorectal tumor"))[1,]
GSE11543_sammsi <- GSE11543_sample1[2,GSE11543_sample1[1,]%in%"microsatellite instable "]
GSE11543_sammss <- GSE11543_sample1[2,GSE11543_sample1[1,]%in%"microsatellite stable "]
GSE11543_msi <- as.matrix(GSE11543_matrix1[,colnames(GSE11543_matrix1)%in%GSE11543_sammsi])
GSE11543_mss <- as.matrix(GSE11543_matrix1[,colnames(GSE11543_matrix1)%in%GSE11543_sammss])

# 4. chemotherapy prediction -----------------------------------------------------------------

GSE11543_predmsi <- data.frame(1:ncol(GSE11543_msi))
GSE11543_predmss <- data.frame(1:ncol(GSE11543_mss))

 
for (i in 1:length(drug)) {
  GSE11543_predmsi[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE11543_msi,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE11543_predmsi)[i] <- drug[i]
}

for (i in 1:length(drug)) {
  GSE11543_predmss[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE11543_mss,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE11543_predmss)[i] <- drug[i]
}

# 5. create table for meta -------------------------------------------------------------

GSE11543_predmsi1 <- as.data.frame(apply(GSE11543_predmsi, 2, function(x)(2^x)))
GSE11543_predmsi1[(nrow(GSE11543_predmsi1)+1),] <- apply(
  GSE11543_predmsi, 2, length_na)
GSE11543_predmsi1[(nrow(GSE11543_predmsi1)+1),] <- apply(
  GSE11543_predmsi1[1:(nrow(GSE11543_predmsi1)-1),], 2, mean, na.rm = T)
GSE11543_predmsi1[(nrow(GSE11543_predmsi1)+1),] <- apply(
  GSE11543_predmsi1[1:(nrow(GSE11543_predmsi1)-2),], 2, var, na.rm = T)
rownames(GSE11543_predmsi1)[(nrow(GSE11543_predmsi1)-2)] <- "n1"
rownames(GSE11543_predmsi1)[(nrow(GSE11543_predmsi1)-1)] <- "mean1"
rownames(GSE11543_predmsi1)[nrow(GSE11543_predmsi1)] <- "sd1"

GSE11543_predmss1 <- as.data.frame(apply(GSE11543_predmss, 2, function(x)(2^x)))
GSE11543_predmss1[(nrow(GSE11543_predmss1)+1),] <- apply(
  GSE11543_predmss, 2, length_na)
GSE11543_predmss1[(nrow(GSE11543_predmss1)+1),] <- apply(
  GSE11543_predmss1[1:(nrow(GSE11543_predmss1)-1),], 2, mean, na.rm = T)
GSE11543_predmss1[(nrow(GSE11543_predmss1)+1),] <- apply(
  GSE11543_predmss1[1:(nrow(GSE11543_predmss1)-2),], 2, var, na.rm = T)
rownames(GSE11543_predmss1)[(nrow(GSE11543_predmss1)-2)] <- "n2"
rownames(GSE11543_predmss1)[(nrow(GSE11543_predmss1)-1)] <- "mean2"
rownames(GSE11543_predmss1)[nrow(GSE11543_predmss1)] <- "sd2"
# GSE13294 ----------------------------------------------------------------
# 2. read data -----------------------------------------------------------------

GSE13294_geo <- read.table("/mnt/DATA/ytj/GEO_meta/GSE13294_series_matrix.txt.gz",
                           header = T, comment.char = "!", sep = "\t")

GSE13294_sample <- read.xlsx2("/mnt/DATA/ytj/GEO_meta/GSE13294_msi.xlsx", sheetIndex = 1,
                              header = F)
# 3. data preprocessing ----------------------------------------------------------------

# MAS5 strandardization after log transformation
GSE13294_matrix <- merge(GSE13294_geo, at_ids, by.x = "ID_REF", by.y = "probe_id")

which(duplicated(GSE13294_matrix[,ncol(GSE13294_matrix)]))

GSE13294_matrix1 <- GSE13294_matrix
GSE13294_matrix1[,(ncol(GSE13294_matrix1)+1)] <- apply(
  GSE13294_matrix1[,2:(ncol(GSE13294_matrix1)-1)], 1, median)
GSE13294_matrix1 <- GSE13294_matrix1[order(
  GSE13294_matrix1$symbol, GSE13294_matrix1$V158, decreasing = T),]
GSE13294_matrix1 <- GSE13294_matrix1[which(!duplicated(GSE13294_matrix1$symbol)),]

rownames(GSE13294_matrix1) <- GSE13294_matrix1$symbol
GSE13294_matrix1 <- GSE13294_matrix1[,-c(1, ncol(GSE13294_matrix1), (ncol(GSE13294_matrix1)-1))]

GSE13294_sample1 <- GSE13294_sample[,-1]
GSE13294_sammsi <- GSE13294_sample1[1,GSE13294_sample1[2,]%in%"MSI"]
GSE13294_sammss <- GSE13294_sample1[1,GSE13294_sample1[2,]%in%"MSS"]
GSE13294_msi <- as.matrix(GSE13294_matrix1[,colnames(GSE13294_matrix1)%in%GSE13294_sammsi])
GSE13294_mss <- as.matrix(GSE13294_matrix1[,colnames(GSE13294_matrix1)%in%GSE13294_sammss])

# 4. chemotherapy prediction -----------------------------------------------------------------

GSE13294_predmsi <- data.frame(1:ncol(GSE13294_msi))
GSE13294_predmss <- data.frame(1:ncol(GSE13294_mss))

 
for (i in 1:length(drug)) {
  GSE13294_predmsi[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE13294_msi,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE13294_predmsi)[i] <- drug[i]
}

for (i in 1:length(drug)) {
  GSE13294_predmss[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE13294_mss,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE13294_predmss)[i] <- drug[i]
}
# 5. create table for meta -------------------------------------------------------------

GSE13294_predmsi1 <- as.data.frame(apply(GSE13294_predmsi, 2, function(x)(2^x)))
GSE13294_predmsi1[(nrow(GSE13294_predmsi1)+1),] <- apply(
  GSE13294_predmsi, 2, length_na)
GSE13294_predmsi1[(nrow(GSE13294_predmsi1)+1),] <- apply(
  GSE13294_predmsi1[1:(nrow(GSE13294_predmsi1)-1),], 2, mean, na.rm = T)
GSE13294_predmsi1[(nrow(GSE13294_predmsi1)+1),] <- apply(
  GSE13294_predmsi1[1:(nrow(GSE13294_predmsi1)-2),], 2, var, na.rm = T)
rownames(GSE13294_predmsi1)[(nrow(GSE13294_predmsi1)-2)] <- "n1"
rownames(GSE13294_predmsi1)[(nrow(GSE13294_predmsi1)-1)] <- "mean1"
rownames(GSE13294_predmsi1)[nrow(GSE13294_predmsi1)] <- "sd1"

GSE13294_predmss1 <- as.data.frame(apply(GSE13294_predmss, 2, function(x)(2^x)))
GSE13294_predmss1[(nrow(GSE13294_predmss1)+1),] <- apply(
  GSE13294_predmss, 2, length_na)
GSE13294_predmss1[(nrow(GSE13294_predmss1)+1),] <- apply(
  GSE13294_predmss1[1:(nrow(GSE13294_predmss1)-1),], 2, mean, na.rm = T)
GSE13294_predmss1[(nrow(GSE13294_predmss1)+1),] <- apply(
  GSE13294_predmss1[1:(nrow(GSE13294_predmss1)-2),], 2, var, na.rm = T)
rownames(GSE13294_predmss1)[(nrow(GSE13294_predmss1)-2)] <- "n2"
rownames(GSE13294_predmss1)[(nrow(GSE13294_predmss1)-1)] <- "mean2"
rownames(GSE13294_predmss1)[nrow(GSE13294_predmss1)] <- "sd2"
# GSE13067 ----------------------------------------------------------------
# 2. read data -----------------------------------------------------------------

GSE13067_geo <- read.table("/mnt/DATA/ytj/GEO_meta/GSE13067_series_matrix.txt.gz",
                           header = T, comment.char = "!", sep = "\t")

GSE13067_sample <- read.xlsx2("/mnt/DATA/ytj/GEO_meta/GSE13067_msi.xlsx", sheetIndex = 1,
                              header = F)
# 3. data preprocessing ----------------------------------------------------------------

GSE13067_matrix <- merge(GSE13067_geo, at_ids, by.x = "ID_REF", by.y = "probe_id")

which(duplicated(GSE13067_matrix[,ncol(GSE13067_matrix)]))

GSE13067_matrix1 <- GSE13067_matrix
GSE13067_matrix1[,(ncol(GSE13067_matrix1)+1)] <- apply(
  GSE13067_matrix1[,2:(ncol(GSE13067_matrix1)-1)], 1, median)
GSE13067_matrix1 <- GSE13067_matrix1[order(
  GSE13067_matrix1$symbol, GSE13067_matrix1$V77, decreasing = T),]
GSE13067_matrix1 <- GSE13067_matrix1[which(!duplicated(GSE13067_matrix1$symbol)),]

rownames(GSE13067_matrix1) <- GSE13067_matrix1$symbol
GSE13067_matrix1 <- GSE13067_matrix1[,-c(1, ncol(GSE13067_matrix1), (ncol(GSE13067_matrix1)-1))]

GSE13067_sample1 <- GSE13067_sample[,-1]
GSE13067_sammsi <- GSE13067_sample1[1,GSE13067_sample1[2,]%in%"MSI"]
GSE13067_sammss <- GSE13067_sample1[1,GSE13067_sample1[2,]%in%"MSS"]
GSE13067_msi <- as.matrix(GSE13067_matrix1[,colnames(GSE13067_matrix1)%in%GSE13067_sammsi])
GSE13067_mss <- as.matrix(GSE13067_matrix1[,colnames(GSE13067_matrix1)%in%GSE13067_sammss])

# 4. chemotherapy prediction -----------------------------------------------------------------

GSE13067_predmsi <- data.frame(1:ncol(GSE13067_msi))
GSE13067_predmss <- data.frame(1:ncol(GSE13067_mss))

 
for (i in 1:length(drug)) {
  GSE13067_predmsi[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE13067_msi,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE13067_predmsi)[i] <- drug[i]
}

for (i in 1:length(drug)) {
  GSE13067_predmss[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE13067_mss,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE13067_predmss)[i] <- drug[i]
}


# 5. create table for meta -------------------------------------------------------------

GSE13067_predmsi1 <- as.data.frame(apply(GSE13067_predmsi, 2, function(x)(2^x)))
GSE13067_predmsi1[(nrow(GSE13067_predmsi1)+1),] <- apply(
  GSE13067_predmsi, 2, length_na)
GSE13067_predmsi1[(nrow(GSE13067_predmsi1)+1),] <- apply(
  GSE13067_predmsi1[1:(nrow(GSE13067_predmsi1)-1),], 2, mean, na.rm = T)
GSE13067_predmsi1[(nrow(GSE13067_predmsi1)+1),] <- apply(
  GSE13067_predmsi1[1:(nrow(GSE13067_predmsi1)-2),], 2, var, na.rm = T)
rownames(GSE13067_predmsi1)[(nrow(GSE13067_predmsi1)-2)] <- "n1"
rownames(GSE13067_predmsi1)[(nrow(GSE13067_predmsi1)-1)] <- "mean1"
rownames(GSE13067_predmsi1)[nrow(GSE13067_predmsi1)] <- "sd1"

GSE13067_predmss1 <- as.data.frame(apply(GSE13067_predmss, 2, function(x)(2^x)))
GSE13067_predmss1[(nrow(GSE13067_predmss1)+1),] <- apply(
  GSE13067_predmss, 2, length_na)
GSE13067_predmss1[(nrow(GSE13067_predmss1)+1),] <- apply(
  GSE13067_predmss1[1:(nrow(GSE13067_predmss1)-1),], 2, mean, na.rm = T)
GSE13067_predmss1[(nrow(GSE13067_predmss1)+1),] <- apply(
  GSE13067_predmss1[1:(nrow(GSE13067_predmss1)-2),], 2, var, na.rm = T)
rownames(GSE13067_predmss1)[(nrow(GSE13067_predmss1)-2)] <- "n2"
rownames(GSE13067_predmss1)[(nrow(GSE13067_predmss1)-1)] <- "mean2"
rownames(GSE13067_predmss1)[nrow(GSE13067_predmss1)] <- "sd2"
# GSE4559 -----------------------------------------------------------------
# 2. read data -----------------------------------------------------------------

GSE4559_geo <- read.table("/mnt/DATA/ytj/GEO_meta/GSE4459_series_matrix.txt.gz",
                          header = T, comment.char = "!", sep = "\t")

GSE4559_sample <- read.xlsx2("/mnt/DATA/ytj/GEO_meta/GSE4459_msi.xlsx", sheetIndex = 1,
                             header = F)
# 3. data preprocessing ----------------------------------------------------------------

GSE4559_matrix <- merge(GSE4559_geo, at_ids, by.x = "ID_REF", by.y = "probe_id")

which(duplicated(GSE4559_matrix[,ncol(GSE4559_matrix)]))

GSE4559_matrix1 <- GSE4559_matrix
GSE4559_matrix1[,(ncol(GSE4559_matrix1)+1)] <- apply(
  GSE4559_matrix1[,2:(ncol(GSE4559_matrix1)-1)], 1, median)
GSE4559_matrix1 <- GSE4559_matrix1[order(
  GSE4559_matrix1$symbol, GSE4559_matrix1$V92, decreasing = T),]
GSE4559_matrix1 <- GSE4559_matrix1[which(!duplicated(GSE4559_matrix1$symbol)),]

rownames(GSE4559_matrix1) <- GSE4559_matrix1$symbol
GSE4559_matrix1 <- GSE4559_matrix1[,-c(1, ncol(GSE4559_matrix1), (ncol(GSE4559_matrix1)-1))]

GSE4559_sample1 <- GSE4559_sample[,-1]
GSE4559_sammsi <- GSE4559_sample1[1,GSE4559_sample1[2,]%in%c(
  "MSI, proximal","MSI, distal","MSI, unknown")]
GSE4559_sammss <- GSE4559_sample1[1,GSE4559_sample1[2,]%in%c(
  "MSS, distal","MSS, proximal")]
GSE4559_msi <- as.matrix(GSE4559_matrix1[,colnames(GSE4559_matrix1)%in%GSE4559_sammsi])
GSE4559_mss <- as.matrix(GSE4559_matrix1[,colnames(GSE4559_matrix1)%in%GSE4559_sammss])

# 4. chemotherapy prediction -----------------------------------------------------------------

GSE4559_predmsi <- data.frame(1:ncol(GSE4559_msi))
GSE4559_predmss <- data.frame(1:ncol(GSE4559_mss))

 
for (i in 1:length(drug)) {
  GSE4559_predmsi[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE4559_msi,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE4559_predmsi)[i] <- drug[i]
}

for (i in 1:length(drug)) {
  GSE4559_predmss[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE4559_mss,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE4559_predmss)[i] <- drug[i]
}

# 5. create table for meta -------------------------------------------------------------

GSE4559_predmsi1 <- as.data.frame(apply(GSE4559_predmsi, 2, function(x)(2^x)))
GSE4559_predmsi1[(nrow(GSE4559_predmsi1)+1),] <- apply(
  GSE4559_predmsi, 2, length_na)
GSE4559_predmsi1[(nrow(GSE4559_predmsi1)+1),] <- apply(
  GSE4559_predmsi1[1:(nrow(GSE4559_predmsi1)-1),], 2, mean, na.rm = T)
GSE4559_predmsi1[(nrow(GSE4559_predmsi1)+1),] <- apply(
  GSE4559_predmsi1[1:(nrow(GSE4559_predmsi1)-2),], 2, var, na.rm = T)
rownames(GSE4559_predmsi1)[(nrow(GSE4559_predmsi1)-2)] <- "n1"
rownames(GSE4559_predmsi1)[(nrow(GSE4559_predmsi1)-1)] <- "mean1"
rownames(GSE4559_predmsi1)[nrow(GSE4559_predmsi1)] <- "sd1"

GSE4559_predmss1 <- as.data.frame(apply(GSE4559_predmss, 2, function(x)(2^x)))
GSE4559_predmss1[(nrow(GSE4559_predmss1)+1),] <- apply(
  GSE4559_predmss, 2, length_na)
GSE4559_predmss1[(nrow(GSE4559_predmss1)+1),] <- apply(
  GSE4559_predmss1[1:(nrow(GSE4559_predmss1)-1),], 2, mean, na.rm = T)
GSE4559_predmss1[(nrow(GSE4559_predmss1)+1),] <- apply(
  GSE4559_predmss1[1:(nrow(GSE4559_predmss1)-2),], 2, var, na.rm = T)
rownames(GSE4559_predmss1)[(nrow(GSE4559_predmss1)-2)] <- "n2"
rownames(GSE4559_predmss1)[(nrow(GSE4559_predmss1)-1)] <- "mean2"
rownames(GSE4559_predmss1)[nrow(GSE4559_predmss1)] <- "sd2"
# GSE2138 -----------------------------------------------------------------
# 2. read data -----------------------------------------------------------------

GSE2138_geo <- read.table("/mnt/DATA/ytj/GEO_meta/GSE2138-GPL96_series_matrix.txt.gz",
                          header = T, comment.char = "!", sep = "\t")

GSE2138_sample <- read.xlsx2("/mnt/DATA/ytj/GEO_meta/GSE2138-GPL96_msi.xlsx", sheetIndex = 1,
                             header = F)
# 3. data preprocessing ----------------------------------------------------------------

at_ids_hgu133a <- toTable(hgu133aSYMBOL)
GSE2138_matrix <- merge(GSE2138_geo, at_ids_hgu133a, by.x = "ID_REF", by.y = "probe_id")

which(duplicated(GSE2138_matrix[,ncol(GSE2138_matrix)]))

GSE2138_matrix1 <- GSE2138_matrix
GSE2138_matrix1[,(ncol(GSE2138_matrix1)+1)] <- apply(
  GSE2138_matrix1[,2:(ncol(GSE2138_matrix1)-1)], 1, median)
GSE2138_matrix1 <- GSE2138_matrix1[order(
  GSE2138_matrix1$symbol, GSE2138_matrix1$V23, decreasing = T),]
GSE2138_matrix1 <- GSE2138_matrix1[which(!duplicated(GSE2138_matrix1$symbol)),]

rownames(GSE2138_matrix1) <- GSE2138_matrix1$symbol
GSE2138_matrix1 <- GSE2138_matrix1[,-c(1, ncol(GSE2138_matrix1), (ncol(GSE2138_matrix1)-1))]

GSE2138_matrix1 <- t(apply(GSE2138_matrix1+1, 1, log2))
range(GSE2138_matrix1)

GSE2138_sample1 <- GSE2138_sample[,-1]
GSE2138_sammsi <- GSE2138_sample1[1,GSE2138_sample1[2,]%in%"msi status: Positive"]
GSE2138_sammss <- GSE2138_sample1[1,GSE2138_sample1[2,]%in%"msi status: Negative"]
GSE2138_msi <- as.matrix(GSE2138_matrix1[,colnames(GSE2138_matrix1)%in%GSE2138_sammsi])
GSE2138_mss <- as.matrix(GSE2138_matrix1[,colnames(GSE2138_matrix1)%in%GSE2138_sammss])

# 4. chemotherapy prediction -----------------------------------------------------------------

GSE2138_predmsi <- data.frame(1:ncol(GSE2138_msi))
GSE2138_predmss <- data.frame(1:ncol(GSE2138_mss))

 
for (i in 1:length(drug)) {
  GSE2138_predmsi[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE2138_msi,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE2138_predmsi)[i] <- drug[i]
}

for (i in 1:length(drug)) {
  GSE2138_predmss[,i] <- pRRophetic::pRRopheticPredict(
    testMatrix = GSE2138_mss,
    drug = drug[i],
    tissueType = "all",
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016",
    minNumSamples = 1
  )
  colnames(GSE2138_predmss)[i] <- drug[i]
}



# 5. create table for meta -------------------------------------------------------------

GSE2138_predmsi1 <- as.data.frame(apply(GSE2138_predmsi, 2, function(x)(2^x)))
GSE2138_predmsi1[(nrow(GSE2138_predmsi1)+1),] <- apply(
  GSE2138_predmsi, 2, length_na)
GSE2138_predmsi1[(nrow(GSE2138_predmsi1)+1),] <- apply(
  GSE2138_predmsi1[1:(nrow(GSE2138_predmsi1)-1),], 2, mean, na.rm = T)
GSE2138_predmsi1[(nrow(GSE2138_predmsi1)+1),] <- apply(
  GSE2138_predmsi1[1:(nrow(GSE2138_predmsi1)-2),], 2, var, na.rm = T)
rownames(GSE2138_predmsi1)[(nrow(GSE2138_predmsi1)-2)] <- "n1"
rownames(GSE2138_predmsi1)[(nrow(GSE2138_predmsi1)-1)] <- "mean1"
rownames(GSE2138_predmsi1)[nrow(GSE2138_predmsi1)] <- "sd1"

GSE2138_predmss1 <- as.data.frame(apply(GSE2138_predmss, 2, function(x)(2^x)))
GSE2138_predmss1[(nrow(GSE2138_predmss1)+1),] <- apply(
  GSE2138_predmss, 2, length_na)
GSE2138_predmss1[(nrow(GSE2138_predmss1)+1),] <- apply(
  GSE2138_predmss1[1:(nrow(GSE2138_predmss1)-1),], 2, mean, na.rm = T)
GSE2138_predmss1[(nrow(GSE2138_predmss1)+1),] <- apply(
  GSE2138_predmss1[1:(nrow(GSE2138_predmss1)-2),], 2, var, na.rm = T)
rownames(GSE2138_predmss1)[(nrow(GSE2138_predmss1)-2)] <- "n2"
rownames(GSE2138_predmss1)[(nrow(GSE2138_predmss1)-1)] <- "mean2"
rownames(GSE2138_predmss1)[nrow(GSE2138_predmss1)] <- "sd2"

# overall meta -----------------------------------------------------------------

samname <- c("GSE103340", "GSE11543", "GSE13067", "GSE13294", "GSE143985",
             "GSE146889", "GSE156915", "GSE18088", "GSE185055", "GSE2138",
             "GSE24550", "GSE24551", "GSE24795", "GSE25071", "GSE26682",
             "GSE27544", "GSE29638", "GSE30378", "GSE35566", "GSE35896",
             "GSE39084", "GSE41258", "GSE4459", "GSE75315",
             "GSE92921")


con_meta <- list()
predict_metadata <- list()
predict_table <- list()
for (i in 1:length(drug)) {
  
  msi_n <- c(GSE103340_predmsi1["n1",i], GSE11543_predmsi1["n1",i], GSE13067_predmsi1["n1",i], GSE13294_predmsi1["n1",i], GSE143985_predmsi1["n1",i],
             GSE146889_predmsi1["n1",i], GSE156915_predmsi1["n1",i], GSE18088_predmsi1["n1",i], GSE185055_predmsi1["n1",i], GSE2138_predmsi1["n1",i],
             GSE24550_predmsi1["n1",i], GSE24551_predmsi1["n1",i], GSE24795_predmsi1["n1",i], GSE25071_predmsi1["n1",i], GSE26682_predmsi1["n1",i],
             GSE27544_predmsi1["n1",i], GSE29638_predmsi1["n1",i], GSE30378_predmsi1["n1",i], GSE35566_predmsi1["n1",i], GSE35896_predmsi1["n1",i], 
             GSE39084_predmsi1["n1",i], GSE41258_predmsi1["n1",i], GSE4459_predmsi1["n1",i], GSE75315_predmsi1["n1",i],
             GSE92921_predmsi1["n1",i])
  msi_mean <- c(GSE103340_predmsi1["mean1",i], GSE11543_predmsi1["mean1",i], GSE13067_predmsi1["mean1",i], GSE13294_predmsi1["mean1",i], GSE143985_predmsi1["mean1",i],
                GSE146889_predmsi1["mean1",i], GSE156915_predmsi1["mean1",i], GSE18088_predmsi1["mean1",i], GSE185055_predmsi1["mean1",i], GSE2138_predmsi1["mean1",i],
                GSE24550_predmsi1["mean1",i], GSE24551_predmsi1["mean1",i], GSE24795_predmsi1["mean1",i], GSE25071_predmsi1["mean1",i], GSE26682_predmsi1["mean1",i],
                GSE27544_predmsi1["mean1",i], GSE29638_predmsi1["mean1",i], GSE30378_predmsi1["mean1",i], GSE35566_predmsi1["mean1",i], GSE35896_predmsi1["mean1",i], 
                GSE39084_predmsi1["mean1",i], GSE41258_predmsi1["mean1",i], GSE4459_predmsi1["mean1",i], GSE75315_predmsi1["mean1",i],
                GSE92921_predmsi1["mean1",i])
  msi_sd <- c(GSE103340_predmsi1["sd1",i], GSE11543_predmsi1["sd1",i], GSE13067_predmsi1["sd1",i], GSE13294_predmsi1["sd1",i], GSE143985_predmsi1["sd1",i],
              GSE146889_predmsi1["sd1",i], GSE156915_predmsi1["sd1",i], GSE18088_predmsi1["sd1",i], GSE185055_predmsi1["sd1",i], GSE2138_predmsi1["sd1",i],
              GSE24550_predmsi1["sd1",i], GSE24551_predmsi1["sd1",i], GSE24795_predmsi1["sd1",i], GSE25071_predmsi1["sd1",i], GSE26682_predmsi1["sd1",i],
              GSE27544_predmsi1["sd1",i], GSE29638_predmsi1["sd1",i], GSE30378_predmsi1["sd1",i], GSE35566_predmsi1["sd1",i], GSE35896_predmsi1["sd1",i], 
              GSE39084_predmsi1["sd1",i], GSE41258_predmsi1["sd1",i], GSE4459_predmsi1["sd1",i], GSE75315_predmsi1["sd1",i],
              GSE92921_predmsi1["sd1",i])
  mss_n <- c(GSE103340_predmss1["n2",i], GSE11543_predmss1["n2",i], GSE13067_predmss1["n2",i], GSE13294_predmss1["n2",i], GSE143985_predmss1["n2",i],
             GSE146889_predmss1["n2",i], GSE156915_predmss1["n2",i], GSE18088_predmss1["n2",i], GSE185055_predmss1["n2",i], GSE2138_predmss1["n2",i],
             GSE24550_predmss1["n2",i], GSE24551_predmss1["n2",i], GSE24795_predmss1["n2",i], GSE25071_predmss1["n2",i], GSE26682_predmss1["n2",i],
             GSE27544_predmss1["n2",i], GSE29638_predmss1["n2",i], GSE30378_predmss1["n2",i], GSE35566_predmss1["n2",i], GSE35896_predmss1["n2",i], 
             GSE39084_predmss1["n2",i], GSE41258_predmss1["n2",i], GSE4459_predmss1["n2",i], GSE75315_predmss1["n2",i],
             GSE92921_predmss1["n2",i])
  mss_mean <- c(GSE103340_predmss1["mean2",i], GSE11543_predmss1["mean2",i], GSE13067_predmss1["mean2",i], GSE13294_predmss1["mean2",i], GSE143985_predmss1["mean2",i],
                GSE146889_predmss1["mean2",i], GSE156915_predmss1["mean2",i], GSE18088_predmss1["mean2",i], GSE185055_predmss1["mean2",i], GSE2138_predmss1["mean2",i],
                GSE24550_predmss1["mean2",i], GSE24551_predmss1["mean2",i], GSE24795_predmss1["mean2",i], GSE25071_predmss1["mean2",i], GSE26682_predmss1["mean2",i],
                GSE27544_predmss1["mean2",i], GSE29638_predmss1["mean2",i], GSE30378_predmss1["mean2",i], GSE35566_predmss1["mean2",i], GSE35896_predmss1["mean2",i], 
                GSE39084_predmss1["mean2",i], GSE41258_predmss1["mean2",i], GSE4459_predmss1["mean2",i], GSE75315_predmss1["mean2",i],
                GSE92921_predmss1["mean2",i])
  mss_sd <- c(GSE103340_predmss1["sd2",i], GSE11543_predmss1["sd2",i], GSE13067_predmss1["sd2",i], GSE13294_predmss1["sd2",i], GSE143985_predmss1["sd2",i],
              GSE146889_predmss1["sd2",i], GSE156915_predmss1["sd2",i], GSE18088_predmss1["sd2",i], GSE185055_predmss1["sd2",i], GSE2138_predmss1["sd2",i],
              GSE24550_predmss1["sd2",i], GSE24551_predmss1["sd2",i], GSE24795_predmss1["sd2",i], GSE25071_predmss1["sd2",i], GSE26682_predmss1["sd2",i],
              GSE27544_predmss1["sd2",i], GSE29638_predmss1["sd2",i], GSE30378_predmss1["sd2",i], GSE35566_predmss1["sd2",i], GSE35896_predmss1["sd2",i], 
              GSE39084_predmss1["sd2",i], GSE41258_predmss1["sd2",i], GSE4459_predmss1["sd2",i], GSE75315_predmss1["sd2",i],
              GSE92921_predmss1["sd2",i])
  continuousdf <- data.frame(study = samname, n1 = msi_n, mean1 = msi_mean, sd1 = msi_sd,
                             n2 = mss_n, mean2 = mss_mean, sd2 = mss_sd)

  predict_metadata[[i]] <- metacont(msi_n, msi_mean, msi_sd, mss_n, mss_mean, mss_sd,
                                    data = continuousdf,
                                    # ROM refers to ratio of mean，as know as fold change
                                    sm = "ROM",
                                    # T not log，F log
                                    backtransf = F,
                                    # method.smd = "Hedges",
                                    studlab = samname,
  )
  
  predict_table[[i]] <- data.frame(
    "Study" = predict_metadata[[i]][["studlab"]],
    "Total(MSI-H)" = predict_metadata[[i]][["n.e"]],
    "Mean(MSI-H)" = predict_metadata[[i]][["mean.e"]],
    "SD(MSI-H)" = predict_metadata[[i]][["sd.e"]],
    "Total(MSI-H)" = predict_metadata[[i]][["n.c"]],
    "Mean(MSI-H)" = predict_metadata[[i]][["mean.c"]],
    "SD(MSI-H)" = predict_metadata[[i]][["sd.c"]],
    "logFC" = predict_metadata[[i]][["TE"]],
    "95%CI" = paste0("(",
                     predict_metadata[[i]][["lower"]], ";", predict_metadata[[i]][["upper"]], ")")
  )
  
  # Because the output of the RMA function is logified, if we want to de log the result, we need to use predict (, transfer=exp) to convert it to a normal value
  con_meta[[i]] <-rma.uni(n1i = continuousdf$n1,n2i = continuousdf$n2,
                          m1i = continuousdf$mean1, m2i = continuousdf$mean2,
                          sd1i = continuousdf$sd1, sd2i = continuousdf$sd2,
                          # FE refers to Fixed effect model
                          # MD refers to mean difference
                          # ROM refers to ratio of mean
                          data = continuousdf,measure = "ROM",method="FE",slab=samname, digits = 8)
  
}
names(con_meta) <- drug
names(predict_metadata) <- drug
names(predict_table) <- drug


# Sort the forestplot from top to bottom
order(c(con_meta[[1]][["b"]], con_meta[[2]][["b"]], con_meta[[3]][["b"]],
        con_meta[[4]][["b"]], con_meta[[5]][["b"]], con_meta[[6]][["b"]],
        con_meta[[7]][["b"]], con_meta[[8]][["b"]], con_meta[[9]][["b"]],
        con_meta[[10]][["b"]], con_meta[[11]][["b"]], con_meta[[12]][["b"]],
        con_meta[[13]][["b"]], con_meta[[14]][["b"]]), decreasing = T)

# Manually build a text table and use the R package forestplot to draw a forest map
tabletext_predict <- cbind(
  c("Drugs", NA, drug),
  c("logFC(95%CI)", NA, 
    paste0(sprintf("%0.4f", con_meta[[8]][["b"]]), " (", sprintf("%0.4f", con_meta[[8]][["ci.lb"]]), "; ", sprintf("%0.4f", con_meta[[8]][["ci.ub"]]), ")"),
    paste0(sprintf("%0.4f", con_meta[[5]][["b"]]), " (", sprintf("%0.4f", con_meta[[5]][["ci.lb"]]), "; ", sprintf("%0.4f", con_meta[[5]][["ci.ub"]]), ")"), 
    paste0(sprintf("%0.4f", con_meta[[9]][["b"]])," (", sprintf("%0.4f", con_meta[[9]][["ci.lb"]]), "; ", sprintf("%0.4f", con_meta[[9]][["ci.ub"]]), ")"),
    paste0(sprintf("%0.4f", con_meta[[1]][["b"]]), " (", sprintf("%0.4f", con_meta[[1]][["ci.lb"]]), "; ", sprintf("%0.4f", con_meta[[1]][["ci.ub"]]), ")"),
    paste0(sprintf("%0.4f", con_meta[[12]][["b"]]), " (", sprintf("%0.4f", con_meta[[12]][["ci.lb"]]), "; ", sprintf("%0.4f", con_meta[[12]][["ci.ub"]]), ")"),
    paste0(sprintf("%0.4f", con_meta[[6]][["b"]])," (", sprintf("%0.4f", con_meta[[6]][["ci.lb"]]), "; ", sprintf("%0.4f", con_meta[[6]][["ci.ub"]]), ")"),
    paste0(sprintf("%0.4f", con_meta[[4]][["b"]]), " (", sprintf("%0.4f", con_meta[[4]][["ci.lb"]]), "; ", sprintf("%0.4f", con_meta[[4]][["ci.ub"]]), ")"),
    paste0(sprintf("%0.4f", con_meta[[2]][["b"]]), " (", sprintf("%0.4f", con_meta[[2]][["ci.lb"]]), "; ", sprintf("%0.4f", con_meta[[2]][["ci.ub"]]), ")"),
    paste0(sprintf("%0.4f", con_meta[[13]][["b"]])," (", sprintf("%0.4f", con_meta[[13]][["ci.lb"]]), "; ", sprintf("%0.4f", con_meta[[13]][["ci.ub"]]), ")"),
    paste0(sprintf("%0.4f", con_meta[[10]][["b"]]), " (", sprintf("%0.4f", con_meta[[10]][["ci.lb"]]), "; ", sprintf("%0.4f", con_meta[[10]][["ci.ub"]]), ")"),
    paste0(sprintf("%0.4f", con_meta[[11]][["b"]]), " (", sprintf("%0.4f", con_meta[[11]][["ci.lb"]]), "; ", sprintf("%0.4f", con_meta[[11]][["ci.ub"]]), ")"),
    paste0(sprintf("%0.4f", con_meta[[7]][["b"]])," (", sprintf("%0.4f", con_meta[[7]][["ci.lb"]]), "; ", sprintf("%0.4f", con_meta[[7]][["ci.ub"]]), ")"),
    paste0(sprintf("%0.4f", con_meta[[3]][["b"]]), " (", sprintf("%0.4f", con_meta[[3]][["ci.lb"]]), "; ", sprintf("%0.4f", con_meta[[3]][["ci.ub"]]), ")"),
    paste0(sprintf("%0.4f", con_meta[[14]][["b"]]), " (", sprintf("%0.4f", con_meta[[14]][["ci.lb"]]), "; ", sprintf("%0.4f", con_meta[[14]][["ci.ub"]]), ")")),
  c("I2(%)", NA, sprintf("%0.4f", con_meta[[8]][["I2"]]), sprintf("%0.4f", con_meta[[5]][["I2"]]),
    sprintf("%0.4f", con_meta[[9]][["I2"]]), sprintf("%0.4f", con_meta[[1]][["I2"]]), 
    sprintf("%0.4f", con_meta[[12]][["I2"]]), sprintf("%0.4f", con_meta[[6]][["I2"]]),
    sprintf("%0.4f", con_meta[[4]][["I2"]]), sprintf("%0.4f", con_meta[[2]][["I2"]]),
    sprintf("%0.4f", con_meta[[13]][["I2"]]), sprintf("%0.4f", con_meta[[10]][["I2"]]),
    sprintf("%0.4f", con_meta[[11]][["I2"]]), sprintf("%0.4f", con_meta[[7]][["I2"]]),
    sprintf("%0.4f", con_meta[[3]][["I2"]]), sprintf("%0.4f", con_meta[[14]][["I2"]])))

cochrane_from_rmeta <- 
  structure(list(
    mean  = c(NA, NA, con_meta[[8]][["b"]], con_meta[[5]][["b"]], con_meta[[9]][["b"]],
              con_meta[[1]][["b"]], con_meta[[12]][["b"]], con_meta[[6]][["b"]],
              con_meta[[4]][["b"]], con_meta[[2]][["b"]], con_meta[[13]][["b"]],
              con_meta[[10]][["b"]], con_meta[[11]][["b"]], con_meta[[7]][["b"]],
              con_meta[[3]][["b"]], con_meta[[14]][["b"]]), 
    lower = c(NA, NA, con_meta[[8]][["ci.lb"]], con_meta[[5]][["ci.lb"]], con_meta[[9]][["ci.lb"]],
              con_meta[[1]][["ci.lb"]], con_meta[[12]][["ci.lb"]], con_meta[[6]][["ci.lb"]],
              con_meta[[4]][["ci.lb"]], con_meta[[2]][["ci.lb"]], con_meta[[13]][["ci.lb"]],
              con_meta[[10]][["ci.lb"]], con_meta[[11]][["ci.lb"]], con_meta[[7]][["ci.lb"]],
              con_meta[[3]][["ci.lb"]], con_meta[[14]][["ci.lb"]]),
    upper = c(NA, NA, con_meta[[8]][["ci.ub"]], con_meta[[5]][["ci.ub"]], con_meta[[9]][["ci.ub"]],
              con_meta[[1]][["ci.ub"]], con_meta[[12]][["ci.ub"]], con_meta[[6]][["ci.ub"]],
              con_meta[[4]][["ci.ub"]], con_meta[[2]][["ci.ub"]], con_meta[[13]][["ci.ub"]],
              con_meta[[10]][["ci.ub"]], con_meta[[11]][["ci.ub"]], con_meta[[7]][["ci.ub"]],
              con_meta[[3]][["ci.ub"]], con_meta[[14]][["ci.ub"]])),
    .Names = c("mean", "lower", "upper"), row.names = c(NA, -11L), class = "data.frame")

# calling Arial fonts
theme(text = element_text(family = "Arial"))


forestplot(tabletext_predict, cochrane_from_rmeta, 
           # forest plot position
           graph.pos = 3,
           # Add horizontal lines
           hrzl_lines = list("3" = gpar(lty=2, col = "#000044")), 
           # zero line
           zero = -2,
           # zero line width
           lwd.zero = 2.5,
           title = "Meta analysis of differences in chemotherapeutic response by MSI status",
           # x-axis
           xticks = c(-5, -4, -3, -2, -1, 0, 1), lwd.xaxis = 2,
           xlab = "logFC(MSI-H/MSS/MSI-L)", 
           new_page = TRUE, col = fpColors(box="royalblue",line="darkblue"),
           # line
           # Line type and width of confidence interval
           lty.ci = 1, lwd.ci = 1.5,          
           # Add small vertical lines at both ends of the confidence interval
           ci.vertices = TRUE,  
           # The height of the small vertical line
           ci.vertices.height = 0.2, 
           boxsize = 0.5,
           # Align Text Left
           align = "l",
           # The length of the vector is equal to the number of rows in the chart. If true, the rows are bold, and a line is added below the row, but it is not displayed when the color is not set
           is.summary= c(T,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F),  
           txt_gp = fpTxtGp(label = gpar(cex = 1.15, fontfamily="Arial"), # Text font size
                            ticks = gpar(cex = 1.0), # Axis scale size
                            xlab = gpar(cex = 1.2))) # Axis lable Size
dev.off()




