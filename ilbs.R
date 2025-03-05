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

setwd("/home/mrna/Documents/TCGA-Assembler-2-master/TCGA-Assembler")
source("Module_A.R")

RNASeqRawData<-DownloadRNASeqData(cancerType = "LIHC", assayPlatform = "gene.normalized_RNAseq", saveFolderName = "./")

source("Module_B.R")

GeneExpData<-ProcessRNASeqData(inputFilePath = RNASeqRawData[1],
                                          outputFileName ="READ__illuminahiseq_rnaseqv2__GeneExp", outputFileFolder = "./", dataType = "geneExp", verType = "RNASeqV2")


# Query platform Illumina HiSeq with a list of barcode 
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
library(writexl)
write_xlsx(dataDEGs,"/home/mrna/Documents/TCGA-Assembler-2-master/TCGA-Assembler/DEG.xlsx")
write_xlsx(dataDEGsFiltLevel,"/home/mrna/Documents/TCGA-Assembler-2-master/TCGA-Assembler/DEG_NT.xlsx")








