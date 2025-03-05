setwd("D:/MIT_PhD/DESeq")

install.packages("writexl")
install.packages("readxl")
install.packages("survminer")
#install.packages("dplyr")
#install.packages("CNTools")
#install.packages("TCGA2STAT_1.0.tar.gz", repos = NULL, type = "source")
#########install the package#############
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")

#if (!requireNamespace("BiocManager"))
#  install.packages("BiocManager")
#BiocManager::install("RTCGAToolbox")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("CNTools")

library(TCGAbiolinks)
library(survminer)
library(writexl)
library(readxl)
#library(dplyr)
#library(RTCGAToolbox)
#library(CNTools)
#library(TCGA2STAT)
#Sys.setenv(TAR="C:/cygwin64/bin/tar", R_GZIPCMD="C:/cygwin64/bin/gzip")
########Download the data#############
#rnaseq_os.ov <- getTCGA(disease="OV", data.type="RNASeq", type="RPKM", clinical=TRUE)

#getFirehoseDatasets()
#readData<-getFirehoseData(dataset="OV", runDate=NULL,forceDownload = FALSE,
#                          Clinic=TRUE, Mutation=FALSE, Methylation=FALSE, RNASeq2GeneNorm=TRUE)

#data(readData)
#getData(readData, "RNASeq2GeneNorm")
#OVdata <- getFirehoseData(dataset="OV",
#                           runDate="20140416",gistic2Date="20140115",
#                           RNASeqGene=TRUE,clinical=TRUE,mRNAArray=TRUE,Mutation=FALSE)

#View(OVdata)
CancerProject<-"TCGA-OV"
DataDir<-paste0("../GDC/",gsub("-","_",CancerProject))
FileNameData<-paste0(DataDir, "_","HTSeq_Counts",".rda")

query<-GDCquery(project=CancerProject,
                data.category = "Transcriptome Profiling",
                data.type="Gene Expression Quantification",
                workflow.type="HTSeq - Counts")

samplesDown<-getResults(query,cols = c("cases"))
typeof(samplesDown)

###################Selection of tumor samples###############
dataSmTp<-TCGAquery_SampleTypes(barcode=samplesDown,typesample = "TP")

###################Selection of normal samples###############
dataSmRC<-TCGAquery_SampleTypes(barcode=samplesDown,typesample = "TR")

dataSmTp_short<-dataSmTp[1:374]
dataSmRC_short<-dataSmRC[1:5]

############Download data for selected barcodes################
queryDown<-GDCquery(project=CancerProject,
                    data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    workflow.type = "HTSeq - Counts",
                    barcode=c(dataSmTp_short, dataSmRC_short))

#GDCdownload(query = queryDown)

#dataPrep=GDCprepare(query = queryDown,save=TRUE)
#dataPrep<-TCGAanalyze_Preprocessing(object = dataPrep,cor.cut = 0.6,datatype = "HTseq - Counts")
GDCdownload(queryDown,method = "api", files.per.chunk = 379)
data<-GDCprepare(queryDown)
data<-TCGAanalyze_Preprocessing(data,cor.cut = 0.6,datatype = "HTSeq - Counts")
#df1<-data.frame(data)
#typeof(df)
#write_xlsx(df1,"D:/MIT_PhD/DESeq/Expr_data.xlsx")
#clinical <- data.frame(data@colData)

############Normalization of genes#############
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("EDASeq")
library(EDASeq)
dataNorm<-TCGAanalyze_Normalization(tabDF = data,geneInfo = geneInfoHT,method = "gcContent")
typeof(dataNorm)

#df7<-c(dataNorm)
df3<-as.data.frame(dataNorm)
typeof(df3)
View(dataNorm)
View(df3)
is.data.frame(df3)
#write_xlsx(df3,"D:/MIT_PhD/DESeq/datanorm_3.xlsx")
write.csv(df3,"D:/MIT_PhD/DESeq/datanorm_4.csv")
#Merging expr and clinical data

df4=read.csv("D:/MIT_PhD/DESeq/Clinical info_4.csv")
View(df4)
merge_df<-merge(df3,df4)
#merge_df2<-full_join(df3,df4)
df5<-join

boxplot(data,outline=FALSE)
boxplot(dataNorm,outline = FALSE)

BiocManager::install("edgeR")
library('edgeR')

#####Quantile filter of genes##############
dataFilt<-TCGAanalyze_Filtering(tabDF = dataNorm,method = "quantile",qnt.cut = 0.25)
df2<-data.frame(dataFilt)
write_xlsx(df2,"D:/MIT_PhD/DESeq/Expr_data_2.xlsx")

#############Differential expression analysis###################
dataDEG<-TCGAanalyze_DEA(mat1 = dataFilt[,dataSmTp_short],
                         mat2 =  dataFilt[,dataSmRC_short],
                         Cond1type = "Tumor",
                         Cond2type = "Recurrent",
                         method = "glmLRT"
)

dataDEG<-TCGAanalyze_DEA(mat1 = dataFilt[,dataSmTp_short],
                         mat2 =  dataFilt[,dataSmRC_short],
                         Cond1type = "Tumor",
                         Cond2type = "Recurrent",
                         fdr.cut = 0.5,
                         method = "glmLRT"
)                         
typeof(dataDEG)
library(RColorBrewer)
hist(dataDEG$PValue,col=brewer.pal(3,name="Set2")[1],
     main = "Primary tumor vs Recurrent tumor",xlab = "p-values"
)

##################Volcano Plot######################
dataDEG$diffexpressed<-"No"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
dataDEG$diffexpressed[dataDEG$logFC>2 & dataDEG$PValue<0.05]<- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
dataDEG$diffexpressed[dataDEG$logFC<2 & dataDEG$PValue<0.05]<- "DOWN"
library(ggplot2)
library(ggrepel)
# Re-plot but this time color the points with "diffexpressed"
p<-ggplot(data = dataDEG,aes(x=logFC,y=-log10(PValue),col=diffexpressed))+geom_point()+theme_minimal()
p
# Add lines as before...
p2<-p+geom_vline(xintercept = c(-0.25, 0.25),col="red")+geom_hline(yintercept = -log10(0.05),col="red")
p2
mycolors<-c("blue","red","black")
names(mycolors)<c("DOWN","UP","NO")
p3<-p2+scale_color_manual(values = mycolors)
p3

table(dataDEG$diffexpressed)

dataDEG$SYMBOL<-rownames(dataDEG)
dataDEG$delabel<-NA
dataDEG$delabel[dataDEG$diffexpressed!="NO"]<-dataDEG$SYMBOL[dataDEG$diffexpressed!="NO"]
write.csv(dataDEG, "DE_genes.csv")
BiocManager::install("ggrepel")
#library(ggrepel)
# plot adding up all layers we have seen so far
png("age_volcano.png",width = 20,height = 10,units = 'in',res = 300)

p4<-ggplot(data = dataDEG,aes(x=logFC,y=-log10(PValue),col=diffexpressed,label=delabel))+
  geom_point()+
  theme_minimal()+
  geom_text_repel()+
  scale_color_manual(values = c("blue","black","red"))+
  geom_vline(xintercept = c(-0.25,0.25),col="red")+
  geom_hline(yintercept = -log10(0.05),col="red")
p4

#######################SURVIVAl ANALYSIS###################
clin.ov <- GDCquery_clinic("TCGA-OV", "clinical")
TCGAanalyze_survival(clin.ov,
                     "figo_stage",
                     main = "TCGA Set\n OV",height = 10, width=10, legend = "Clinical stage")

TCGAanalyze_survival(clin.ov,
                     "age_at_index",
                     main = "TCGA Set\n OV",height = 10, width=10, legend = "Age at index")


TCGAanalyze_survival(clin.ov,
                     "vital_status",
                     main = "TCGA Set\n OV",height = 10, width=10, legend = "Vital status")

TCGAanalyze_survival(clin.ov,
                     "age_at_diagnosis",
                     main = "TCGA Set\n OV",height = 10, width=10, legend = "Age at diagnosis")



######### Correlating gene expression and Survival Analysis  ##########
dataOVcomplete <- log2(OV_rnaseqv2)

tokenStop<- 1

tabSurvKMcomplete <- NULL

for( i in 1: round(nrow(dataOVcomplete)/100)){
  message( paste( i, "of ", round(nrow(dataOVcomplete)/100)))
  tokenStart <- tokenStop
  tokenStop <-100*i
  tabSurvKM<-TCGAanalyze_SurvivalKM(clin.ov,
                                    dataOVcomplete,
                                    Genelist = rownames(dataOVcomplete)[tokenStart:tokenStop],
                                    Survresult = F,
                                    ThreshTop=0.67,
                                    ThreshDown=0.33)
  
  tabSurvKMcomplete <- rbind(tabSurvKMcomplete,tabSurvKM)
}

tabSurvKMcomplete <- tabSurvKMcomplete[tabSurvKMcomplete$pvalue < 0.01,]
typeof(tabSurvKMcomplete)
tabSurvKMcomplete <- tabSurvKMcomplete[order(tabSurvKMcomplete$pvalue, decreasing=F),]

tabSurvKMcompleteDEGs <- tabSurvKMcomplete[
  rownames(tabSurvKMcomplete) %in% dataDEGsFiltLevel$mRNA,
  ]





dev.off()










