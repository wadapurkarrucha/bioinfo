setwd("/home/mrna/Documents")

if (!require("BioCManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(version = "3.17")

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

if (!require("BioCManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")

library(TCGAbiolinks)
library(SummarizedExperiment)
library(edgeR)

query.miRNA <- GDCquery(
  project = "TCGA-LIHC",
  experimental.strategy = "miRNA-Seq",
  data.category = "Transcriptome Profiling", 
  data.type = "miRNA Expression Quantification"
)

GDCdownload(query = query.miRNA)

dataAssy.miR <- GDCprepare(
  query = query.miRNA
)
rownames(dataAssy.miR) <- dataAssy.miR$miRNA_ID

# using read_count's data 
read_countData <-  colnames(dataAssy.miR)[grep("count", colnames(dataAssy.miR))]
dataAssy.miR <- dataAssy.miR[,read_countData]
colnames(dataAssy.miR) <- gsub("read_count_","", colnames(dataAssy.miR))

dataFilt <- TCGAanalyze_Filtering(
  tabDF = dataAssy.miR,
  method = "quantile", 
  qnt.cut =  0.25
)   

# selection of normal samples "NT"
dataSmNT_miR <- TCGAquery_SampleTypes(
  barcode = colnames(dataFilt),
  typesample = c("NT")
)

# selection of tumor samples "TP"
dataSmTP_miR <- TCGAquery_SampleTypes(
  barcode = colnames(dataFilt), 
  typesample = c("TP")
)

dataDEGs <- TCGAanalyze_DEA(
  mat1 = dataFilt[,dataSmNT_miR],
  mat2 = dataFilt[,dataSmTP_miR],
  Cond1type = "Normal",
  Cond2type = "Tumor",
  fdr.cut = 0.01 ,
  logFC.cut = 1,
  method = "glmLRT"
)  


BiocManager::install("EDASeq")

install.packages("ggplot")
install.packages("ggrepel")

dataDEGs$diffexpressed<-"No"
# if log2Foldchange &gt; 0.6 and pvalue &lt; 0.05, set as &quot;UP&quot;
dataDEGs$diffexpressed[dataDEGs$logFC>1 & dataDEGs$PValue<0.05] <-"UP"
# if log2Foldchange &lt; -0.6 and pvalue &lt; 0.05, set as &quot;DOWN&quot;
dataDEGs$diffexpressed[dataDEGs$logFC< -1 & dataDEGs$PValue<0.05] <-"DOWN"
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
write.csv(dataDEGs, "DE_miRNA_NEW.csv")

# plot adding up all layers we have seen so far
png("age_volcano_28.08.png",width = 20,height = 10,units = 'in',res = 300)

#p4<-ggplot(data = dataDEGs,aes(x=logFC,y=-log10(PValue),col=diffexpressed,label=delabel))+
#geom_point()+
#theme_minimal()+
#geom_text_repel()+
#scale_color_manual(values = c("blue","black","red"))+
#geom_vline(xintercept = c(-0.25,0.25),col="red")+
# geom_hline(yintercept = -log10(0.05),col="red")
#p4

library(writexl)
write_xlsx(dataDEGs,"/home/mrna/Documents/TCGA-Assembler-2-master/TCGA-Assembler/DE_miRNA_28.08.xlsx")
write_xlsx(dataDEGsFiltLevel,"/home/mrna/Documents/TCGA-Assembler-2-master/TCGA-Assembler/DEG_NT_28.08.xlsx")












