setwd("/home/mrna/Documents")

if (!require("BioCManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(version = "3.16")

install.packages("RCurl")
install.packages("curl")
install.packages("XML")
install.packages("xml2")
install.packages("httr")

if (!require("BioCManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biostrings")

if (!require("BioCManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("BiocFileCache")

if (!require("BioCManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GenomicRanges")

if (!require("BioCManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("rvest")

if (!require("BioCManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("KEGGREST")

if (!require("BioCManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SummarizedExperiment")

if (!require("BioCManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("AnnotationDbi")

if (!require("BioCManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")

if (!require("BioCManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")

library(TCGAbiolinks)
library(SummarizedExperiment)

setwd("/home/mrna/Documents/TCGA-Assembler-2-master/TCGA-Assembler")
source("Module_A.R")

RNASeqRawData<-DownloadRNASeqData(cancerType = "LIHC", assayPlatform = "gene.normalized_RNAseq", saveFolderName = "./")

source("Module_B.R")

GeneExpData<-ProcessRNASeqData(inputFilePath = RNASeqRawData[1],
                                          outputFileName ="READ__illuminahiseq_rnaseqv2__GeneExp", outputFileFolder = "./", dataType = "geneExp", verType = "RNASeqV2")


# Query platform Illumina HiSeq 
query <- GDCquery(
  project = "TCGA-LIHC", 
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification"
  )

# Download with platform IlluminaHiSeq_RNASeqV2
GDCdownload(query)

# Prepare expression matrix with geneID in the rows and samples (barcode) in the columns
# rsem.genes.results as values
LIHC.Rnaseq.SE <- GDCprepare(query)

LIHCMatrix <- assay(LIHC.Rnaseq.SE,"unstranded") 
# For gene expression if you need to see a boxplot correlation and AAIC plot to define outliers you can run
LIHC.RNAseq_CorOutliers <- TCGAanalyze_Preprocessing(LIHC.Rnaseq.SE)

# normalization of genes
dataNorm <- TCGAanalyze_Normalization(
  tabDF = LIHC.RNAseq_CorOutliers, 
  geneInfo =  geneInfoHT
)

# quantile filter of genes
dataFilt <- TCGAanalyze_Filtering(
  tabDF = dataNorm,
  method = "quantile", 
  qnt.cut =  0.25
)

# selection of normal samples "NT"
samplesNT <- TCGAquery_SampleTypes(
  barcode = colnames(dataFilt),
  typesample = c("NT")
)

# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(
  barcode = colnames(dataFilt), 
  typesample = c("TP")
)

# Diff.expr.analysis (DEA)
dataDEGs <- TCGAanalyze_DEA(
  mat1 = dataFilt[,samplesNT],
  mat2 = dataFilt[,samplesTP],
  Cond1type = "Normal",
  Cond2type = "Tumor",
  fdr.cut = 0.01 ,
  logFC.cut = 1,
  method = "glmLRT"
)

# DEGs table with expression values in normal and tumor samples
dataDEGsFiltLevel <- TCGAanalyze_LevelTab(
  FC_FDR_table_mRNA = dataDEGs,
  typeCond1 = "Tumor",
  typeCond2 = "Normal",
  TableCond1 = dataFilt[,samplesTP],
  TableCond2 = dataFilt[,samplesNT]
)

BiocManager::install("EDASeq")

install.packages("ggplot")
install.packages("ggrepel")

dataDEGs$diffexpressed<-"No"
# if log2Foldchange &gt; 0.6 and pvalue &lt; 0.05, set as &quot;UP&quot;
dataDEGs$diffexpressed[dataDEGs$logFC>1 & dataDEGs$PValue<0.05] <-"UP"
# if log2Foldchange &lt; -0.6 and pvalue &lt; 0.05, set as &quot;DOWN&quot;
dataDEGs$diffexpressed[dataDEGs$logFC<1 & dataDEGs$PValue<0.05] <-"DOWN"
library(ggplot2)
library(ggrepel)

# Re-plot but this time color the points with &quot;diffexpressed&quot;
p<-ggplot(data = dataDEGs,aes(x=logFC,y=-
                                  log10(PValue),col=diffexpressed))+geom_point()+theme_minimal()
p

# Add lines as before...
p2<-p+geom_vline(xintercept = c(-0.25, 0.25),col="red")+geom_hline(yintercept = -log10(0.05),col="red")
p2
mycolors<-c("blue","red","black")
names(mycolors)<c("DOWN","UP","NO")
p3<-p2+scale_color_manual(values = mycolors)
p3

table(dataDEGs$diffexpressed)

dataDEGs$SYMBOL<-rownames(dataDEGs)
dataDEGs$delabel<-NA
dataDEGs$delabel[dataDEGs$diffexpressed!="NO"]<-dataDEGs$SYMBOL[dataDEGs$diffexpressed!="NO"]
write.csv(dataDEGs, "DE_genes.csv")

# plot adding up all layers we have seen so far
png("age_volcano.png",width = 20,height = 10,units = 'in',res = 300)

#p4<-ggplot(data = dataDEGs,aes(x=logFC,y=-log10(PValue),col=diffexpressed,label=delabel))+
  #geom_point()+
  #theme_minimal()+
  #geom_text_repel()+
  #scale_color_manual(values = c("blue","black","red"))+
  #geom_vline(xintercept = c(-0.25,0.25),col="red")+
 # geom_hline(yintercept = -log10(0.05),col="red")
#p4

library(writexl)
write_xlsx(dataDEGs,"/home/mrna/Documents/TCGA-Assembler-2-master/TCGA-Assembler/DEG.xlsx")
write_xlsx(dataDEGsFiltLevel,"/home/mrna/Documents/TCGA-Assembler-2-master/TCGA-Assembler/DEG_NT.xlsx")








