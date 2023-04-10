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
# using ddr gene sets from MSigDB creat ddr_final.gmt
gsea <- getGmt("/mnt/DATA/ytj/GSEA/ddr_final.gmt")


# 3. data preprocessing ---------------------------------------------------

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

tissue <- c("Primary Tumor" , "Metastatic")

phe1 <- phe[phe$sample_type%in%tissue,]


# match sample
inter_syl <- intersect(phe1[,1],sur[,1])
sur1 <- sur[sur[,1]%in%inter_syl,]
phe2 <- phe1[phe1[,1]%in%inter_syl,]
# check duplication
sur1[duplicated(sur1[,1]),]

# some patients have more than one sample data and TCGA barcode
b <- sur1[duplicated(sur1[,3]),]
# sorted the most suitable TCGA sample for
c <- separate(data = b , col = 1 ,
              into = c("1","2","3","4") , sep = "-")

d <- separate(data = sur1 , col = 1 ,
              into = c("1","2","3","4") , sep = "-")

e <- d[d[,6]%in%c[,6],]

patient <- unique(e[,6])
g <- data.frame()
for (i in 1:length(patient)) {
  f <- e[e$X_PATIENT%in%patient[i],]
  if("01A"%in%f$`4`){
    g[i,1:7] <- f[f$`4`%in%"01A",]  
  }
  else{if("06A"%in%f$`4`){
    g[i,1:7] <- f[f$`4`%in%"06A",]}
    else{if("01B"%in%f$`4`){
      g[i,1:7] <- f[f$`4`%in%"01B",]}
      else{if("01C"%in%f$`4`){
        g[i,1:7] <- f[f$`4`%in%"01C",]}
        else{if("01R"%in%f$`4`){
          g[i,1:7] <- f[f$`4`%in%"01R",]}else{}}}}}
}
h <- sur1[!sur1[,3]%in%b[,3],]
j <- unite(g, col = "sample" , sep = "-" , 1:4)
sur2 <- rbind(h,j)

inter_syl2 <- intersect(sur2[,1] , phe2[,1])
inter_syl3 <- intersect(colnames(log_tpm) , inter_syl2)
sur3 <- sur2[sur2[,1]%in%inter_syl3,]
phe3 <- phe2[phe2[,1]%in%inter_syl3,]
log_tpm1 <- log_tpm[,inter_syl3]

cancer <- unique(phe3$project_id)

cancer_syl <- list()
for (i in 1:length(cancer)) {
  cancer_syl[[i]] <- phe3[phe3$project_id%in%cancer[i],1]
}

exp_pan <- list()
for (i in 1:length(cancer)) {
  exp_pan[[i]] <- log_tpm1[,colnames(log_tpm1)%in%cancer_syl[[i]]]
}
names(exp_pan) <- cancer

# 3. GSVA -----------------------------------------------------------------


exp_pan <- exp_pan[-c(15,25,33,35)]
cancer <- cancer[-c(15,25,33,35)]
pan_gsva <- list()
for (i in 1:length(cancer)) {
  pan_gsva[[i]] <- gsva(exp_pan[[i]] , gsea , kcdf = "Gaussian", method = "ssgsea", parallel.sz = 45)
}
names(pan_gsva) <- cancer

# combine phe3 and sur3
phe_sur <- cbind(phe3 , sur3 )
rownames(phe_sur) <- phe_sur[,1]
phe_sur <- phe_sur[,-1]

# combine phe_sur and pan_gsva
pan_all <- list()
for (i in 1:length(cancer)) {
  
  exp_pre <- t(pan_gsva[[i]])
  pan_all[[i]] <- cbind(exp_pre , phe_sur[rownames(exp_pre),])
}
names(pan_all) <- cancer



# 4. cox regression ----------------------------------------------------------------

ddr <- rownames(pan_gsva[[1]])

# combine COAD and READ
pan_all1 <- pan_all
TCGA_COAD <- pan_all1$`TCGA-COAD`
TCGA_READ <- pan_all1$`TCGA-READ`
TCGA_CRC <- rbind(TCGA_COAD , TCGA_READ)
TCGA_CRC$project_id <- "TCGA-CRC"
pan_all1 <- pan_all1[-c(27 , 29)]
pan_all1[[31]] <- TCGA_CRC
names(pan_all1)[31] <- "TCGA-CRC"

cancer <- names(pan_all1)

# find best cutoff for the stratification of high and low ssGSEA score
bestcutoff <- function(datavector, clintable) {
  # datavector is a vector of ssGSEA scores of a ddr pathway in all samples 
  # clinitable refers to clinical data, OS and OS.time in specific
  # probe of "quantile" need to be filled with that what kind of percentile do you need to divide this vector into how many parts with "by" as the span
  breaks <- quantile(datavector, probs = seq(0.25, 0.75, by= 0.01))
  # cutoff() is created below
  cutoff.table <- t(sapply(breaks, 
                           function(z) cutoff(datavector = datavector, cutpoint = z, clintable = clintable)))
  colnames(cutoff.table) <- c("cutoff", "pvalue")
  # order, and output the most significant cutoff value to bestcutoff
  cutoff.table[order(cutoff.table[, 2]), "cutoff"][1]
}

# cutpoint is the breaks of bestcutoff()
cutoff <- function(datavector, cutpoint, clintable) {
  
  # Cut is to divide the datavector into several parts according to breaks, and the reason for adding an additional min and max is due to the inherent logic of the function
  # Only in this way can the best cutoff value that comes out last fall on any value in the cutpoint
  # Label refers to the naming of each copy that can be customized for you，such as large, medium and small
  term <- cut(x = datavector, breaks = c(min(datavector), cutpoint, max(datavector)), labels = F, 
              include.lowest = T)
  # test for definition of best cutoff
  cox <- summary(coxph(Surv(OS.time, OS) ~ term, data = clintable))
  # The second one is the p-value, and the output table is the p-value corresponding to a cutoff
  c(cutpoint, cox$sctest[3])
}


pan_all2 <- pan_all1
for (c in 1:length(cancer)) {
  for (i in 1:length(ddr)) {
    
    cutoff.point <-  as.numeric(bestcutoff(datavector = pan_all2[[c]][,i], clintable = pan_all2[[c]][,10:19]))
    # print(cutoff.point)
    for (j in 1:nrow(pan_all2[[c]])){
      
      if (pan_all2[[c]][j,i] >= cutoff.point){
        
        pan_all2[[c]][j,i] = "High"
        
      } else {
        
        pan_all2[[c]][j,i] = "Low"
        
      }
    }
    pan_all2[[c]][,i] = factor( pan_all2[[c]][,i] , levels = c("Low", "High"))
    
  }
}


cox <- list()
bcd <- data.frame()
for (j in 1:length(ddr)) {
  for (i in 1:length(cancer)) {
    
    abc <- summary(coxph(Surv(pan_all2[[i]]$OS.time , pan_all2[[i]]$OS) ~ pan_all2[[i]][,j]))
    hr <- t(as.data.frame(abc$conf.int[,c(1,3,4)]))
    p <- t(as.data.frame(abc$coefficients[,c(3:5)]))
    bcd[i,1:6] <- cbind(hr , p)
    bcd[i,7] <- cancer[i]
  }
  bcd[,8] <- ddr[j]
  colnames(bcd) <- c("exp(coef)" , "lower .95" , "upper .95" ,  "se(coef)" ,  "z" , "p-value" , 
                     "Cancer" , "Pathway")
  # result of cox regression of all tumor except the 33th tumor-CCSK, which is full of null and outlier, are normal 
  bcd <- bcd[-33,]
  cox[[j]] <- bcd
  
}
names(cox) <- ddr


# change name
names(cox)[which(names(cox)=="REACTOME_BASE_EXCISION_REPAIR")] = "BER"
names(cox)[which(names(cox)=="REACTOME_SUMOYLATION_OF_DNA_DAMAGE_RESPONSE_AND_REPAIR_PROTEINS")] = "DDR"
names(cox)[which(names(cox)=="REACTOME_NUCLEOTIDE_EXCISION_REPAIR")] = "NER"
names(cox)[which(names(cox)=="REACTOME_MISMATCH_REPAIR")] = "MMR"
names(cox)[which(names(cox)=="GOMF_SINGLE_STRANDED_DNA_BINDING")] = "SSB"
names(cox)[which(names(cox)=="REACTOME_NONHOMOLOGOUS_END_JOINING_NHEJ")] = "NHEJ"
names(cox)[which(names(cox)=="KEGG_HOMOLOGOUS_RECOMBINATION")] = "HR"
names(cox)[which(names(cox)=="REACTOME_FANCONI_ANEMIA_PATHWAY")] = "FA"
names(cox)[which(names(cox)=="REACTOME_DOUBLE_STRAND_BREAK_REPAIR")] = "DSB"

cox_log <- list()
for (i in 1:length(ddr)) {
  cox[[i]] <- separate(cox[[i]] , col = 7 , into = c("TCGA" , "Cancer") , sep = "-")
  cox[[i]][,9] <- names(cox)[i]
  # solving the problem of ggplot coordinate axis sorting only by letters
  cox_log[[i]] <- cox[[i]]
  cox_log[[i]][,1:3] <- log10(cox_log[[i]][,1:3])
  cox_log[[i]] <- cox_log[[i]][order(cox_log[[i]]$`exp(coef)` , decreasing = F),]
  cox_log[[i]]$Cancer <- factor(cox_log[[i]]$Cancer,levels=cox_log[[i]]$Cancer)
  
}
names(cox_log) <- names(cox)


p <- rbind(cox_log[[1]][,c(1,6,8,9)] , cox_log[[2]][,c(1,6,8,9)] , cox_log[[3]][,c(1,6,8,9)]  , cox_log[[4]][,c(1,6,8,9)] ,
           cox_log[[5]][,c(1,6,8,9)] , cox_log[[6]][,c(1,6,8,9)] , cox_log[[7]][,c(1,6,8,9)] , cox_log[[8]][,c(1,6,8,9)] , 
           cox_log[[9]][,c(1,6,8,9)])



# 4. forest plot -----------------------------------------------------------------


# palette
colors<-c("#015C92","#2D82B5","#53A7D8","#88CDF6",
          "#A5DEE5","#FCAF67","#FA8743","#D53E0F","#7E0306")
# forest plot
plot <- list()
for (i in 1:length(ddr)) {
  
  def  <- ggplot(cox_log[[i]], aes(`exp(coef)`, `Cancer`)) +
    geom_point(size=2.4, col= colors[i]) +
    geom_vline(aes(xintercept = 0)) +
    geom_errorbarh(aes(xmax =`upper .95`, xmin = `lower .95`), height = 0.1) +
    xlab('logHR') + ylab('')+
    coord_cartesian(xlim=c(-1,1))+ # The visible range of the image on the coordinate axis
    scale_x_continuous(breaks=seq(-1, 1, 1))+# The interval displayed on the coordinate axis is 200
    ggtitle(names(cox_log)[[i]]) +
    theme(plot.title = element_text(hjust = 0.4) )+# Title centered
    theme(plot.title = element_text(colour = colors[i]))+
    theme(panel.background = element_rect(fill = "#F4F5EF"))# Change the drawing background color
  plot[[i]] <- def
}


# splicing
forest <- plot_grid(plotlist=list(plot[[1]] , plot[[2]] , plot[[3]] , plot[[4]] , plot[[5]] , plot[[6]] , 
                                  plot[[7]] , plot[[8]] , plot[[9]]), 
                    nrow=1, align='v')

# Waffle Pie Chart
p[,5] <- ""
colnames(p)[5] <- "logHR"
p[,5][which(p[,2]>0.05)] <- "Insignificant"
p[,5][which(p[,1]<0&p[,5]=="")] <- "Protective"
p[,5][which(p[,1]>0&p[,5]=="")] <- "Risky"
p[,5] <- factor(p[,5] , levels = c("Insignificant" , "Protective" , "Risky"))
p_plot <- ggplot(p, aes(`Cancer`, `Pathway`)) +
  geom_tile(aes(fill = `logHR`), colour = "white")+
  scale_fill_manual(values=c("#F6F2E8","#08A0B5","#FAD6A8"))+
  theme(axis.text.x = element_text(angle=90))+
  theme(axis.text.y = element_text(size = 7))+
  theme(panel.background = element_rect(fill = "#F4F5EF"))# Change the drawing background color

plot_grid(p_plot , forest , nrow = 2 , rel_heights = c(1,4))




