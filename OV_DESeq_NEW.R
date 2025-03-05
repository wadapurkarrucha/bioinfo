setwd("D:/MIT_PhD/DESeq")

#########install the package#############
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)

########Download the data#############

CancerProject<-"TCGA-OV"
DataDir<-paste0("../GDC/",gsub("-","_",CancerProject))
FileNameData<-paste0(DataDir, "_","HTSeq_Counts",".rda")

query<-GDCquery(project=CancerProject,
                data.category = "Transcriptome Profiling",
                data.type="Gene Expression Quantification",
                workflow.type="HTSeq - Counts")

samplesDown<-getResults(query,cols = c("cases"))

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

############Normalization of genes#############
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("EDASeq")
library(EDASeq)
dataNorm<-TCGAanalyze_Normalization(tabDF = data,geneInfo = geneInfoHT,method = "gcContent")

boxplot(data,outline=FALSE)
boxplot(dataNorm,outline = FALSE)

BiocManager::install("edgeR")
library('edgeR')

#####Quantile filter of genes##############
dataFilt<-TCGAanalyze_Filtering(tabDF = dataNorm,method = "quantile",qnt.cut = 0.25)

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

dev.off()











