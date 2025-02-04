version
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

install.packages(devtools)
BiocManager::install(c("RColorBrewer", "pheatmap","tidyverse","minfi","minfiData","remotes","shinyMethyl","shinyMethylData","limma"))
BiocManager::install("jokergoo/IlluminaHumanMethylationEPICv2manifest")
BiocManager::install("liftOver")
BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
BiocManager::install("liftOver")
BiocManager::install("DMRcate")
BiocManager::install("EnhancedVolcano")
BiocManager::install("missMethyl")
install.packages("gt")
BiocManager::install("ggrepel")
devtools::install_github('Mikata-Project/ggthemr')
install.packages("~/Desktop/PhD/2024 First Term/DNA Methylation/jokergoo-IlluminaHumanMethylationEPICv2anno.20a1.hg38-702e52a.tar.gz", repos = NULL, type = "source")

#loading those libraries
library(pheatmap)
library(gt)
library(tidyverse)
library(minfi)
library(limma)
library(missMethyl)
library(minfiData)
library(DMRcate)
library(shinyMethyl)
library(shinyMethylData)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(IlluminaHumanMethylationEPICv2manifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(RColorBrewer)
library(dplyr)
library(readr)
library(ggrepel)
library(ggthemr)
library(EnhancedVolcano)

library(devtools)
install_github("achilleasNP/IlluminaHumanMethylationEPICmanifest") 
install_github("achilleasNP/IlluminaHumanMethylationEPICanno.ilm10b5.hg38")


# --------------------------- SHORT EXAMPLE------------------------------#
#you can skip to the detailed analysis

#Setting the working directory
setwd("/Users/egeezen/Desktop/PhD/2024 First Term/DNA Methylation/EPICv2_samples")

#load example data
baseDir <- system.file("extdata", package = "minfiData")

#2 folders for each slide and a csv file samplesheet
list.files(baseDir)

#the folder contains green and red idat files, each slide here has 3 arrays so 6 files
list.files(file.path(baseDir, "5723646052"))

#you can read the directory that returns a dataframe
patients <- read.metharray.sheet(baseDir)

#now the all info is in the RG set
RGset <- read.metharray.exp(targets = patients)
#just to look at it as dataframe (extracts phenotype data)
pd <- pData(RGset)
pd[,1:4]

#you can also just look at the sample sheet, skip=7 skips first 7 rows
#sentrix_id is slide number and sentrix_position is the array number
patients2 <- read.csv(file.path(baseDir, "SampleSheet.csv"), 
                     stringsAsFactors = FALSE, skip=7)
#but there is no basename column, that is why we can manually add it
patients2$Basename <- file.path(baseDir, patients2$Sentrix_ID, 
                               paste0(patients2$Sentrix_ID, 
                                      patients2$Sentrix_Position))
#you can check the annotations used 
annotation(RGset)
                     
#quality control!
qcReport(RGset, sampNames= pd$Sample_Name,
         sampGroups = pd$Sample_Group,pdf="QC_report_sample.pdf")

#you can also use shinyMethly package for further QC
#first create shinymethlyset from your RGset
ShinyMethylSet <- shinySummarize(RGset)

#then you can start the application
runShinyMethyl(ShinyMethylSet)

#additionally MDS plot can be drawn (uses Euclidian distances)
names <- pd$Sample_Name
groups <- pd$Sample_Group
mdsPlot(RGset, sampNames=names, sampGroups=groups,pch=1)

#density plot can be drawn independently
densityPlot(RGset, sampGroups=groups)

#density bean plot can be drawn independently
par(mar=c(5,6,4,2)) #for the alignment of the next figure
densityBeanPlot(RGset, sampNames=names, sampGroups=groups)

#now we can do filtering to get rid of non-specifics
detP <- detectionP(RGset) #detects the p values 
dim(detP)
failed <- detP>0.01
colMeans(failed) # Fraction of failed positions per sample
sum(rowMeans(failed)>0.5) # How many positions failed in >50% of samples?
#filtering out the ones with p-value smaller than 0.01     
absent<- apply(detP, 1, function(x) sum(x>0.01)) #if none of the values are bigger than 0.01 it will return 0  
RGset_filtered <- RGset[absent==0,]
dim(RGset)
dim(RGset_filtered)

#we can d preprocessing here, after filtering!
ssNoob_filtered <- preprocessNoob(RGset_filtered,dyeCorr = TRUE,dyeMethod="single")

#getting methylated and unmethylated sites
methylated_sites <- getMeth(ssNoob_filtered)
unmethylated_sites <- getUnmeth(ssNoob_filtered)

#beta and M (logit of beta) values of the sites
beta_values <- getBeta(ssNoob_filtered)
m_values <- getM(ssNoob_filtered)

#we can plot betas by type(infiniium), only one column at a time can be plotted
plotBetasByType(ssNoob_filtered[,1], main = "R02C02")

#a quality control on methylation object
qc <- getQC(ssNoob_filtered)
ssNoob_filtered <- addQC(ssNoob_filtered, qc = qc)
plotQC(qc)

#if you want to look at specific genomic positions you can plot singleCpGs
cpgs <- c("cg00050873", "cg00212031", "cg26684946", "cg00128718") #examples
par(mfrow=c(2,2)) #to align the plot
par(mar=c(4,5,1,2)) #for the alignment of the next figure
plotCpg(ssNoob_filtered, cpg=cpgs, pheno=groups, type="categorical")

#quality control on methylation object, it does a series of controls and can remove outliers
#if you are going to do this, there is no need to do getQC, getSex or fixMethOutliers
qc_from_minfiQC <- minfiQC(ssNoob_filtered, fixOutliers = TRUE, verbose =TRUE)
plotQC(ssNoob_filtered)

#to find differentially methylated positions from methyl object (F-test for categorical)
dmp <- dmpFinder(ssNoob_filtered, pheno=groups, type="categorical")
sum(dmp$qval < 0.05, na.rm=TRUE) #how many sites are differentially methylated

#you can get the genomic annotation
annot <- as.data.frame(getAnnotation(ssNoob_filtered, dropNonMapping = TRUE))
#only get the annotations for the ones that are significant
significant_dmp <-as.data.frame(annot[rownames(annot) %in% rownames(dmp)[dmp$qval<0.05],])
#to look at a specific gene
CNIH3_gene <- significant_dmp[grep("CNIH3", significant_dmp$UCSC_RefGene_Name), ]#to plot the associated loci
cpgs <- rownames(CNIH3_gene)
par(mfrow=c(3,3)) #to align the plot
par(mar=c(2,2,2,2)) #for the alignment of the next figure
plotCpg(ssNoob_filtered, cpg=cpgs, pheno=groups, type="categorical")

#now we can map to genome to to obtain GenomicMethyl object
genomic_methly_set <- mapToGenome(ssNoob_filtered)

#we can obtain GenomicRatio Object from GenomicMetyl object (what argument for keeping beta and M values)
genomic_ratio_set <- ratioConvert(genomic_methly_set, what = "both")

##you can remove the the probes with SNPs 
genomic_ratio_set_filtered <- dropLociWithSnps(genomic_ratio_set)

#you can remove X and Y chromosomes from the data
#first get the annotation of the genome
#then pick the X and Y chromosomes
#keep <- !(featureNames(genomic_ratio_set_filtered) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
#table(keep)
#genomic_ratio_set_filtered <- genomic_ratio_set_filtered[keep,]

#you can calculate the beta values or the m values after this as well! 
#you can do the annotations here as well 
#annot <- getAnnotation(genomic_ratio_set,dropNonMapping = TRUE)

# --------------------------- SHORT EXAMPLE IS OVER -----------------------------#

#detailed analysis stars here but you can skip to reading final_m_values If you already have the m-values#

setwd("/Users/egeezen/Desktop/PhD/2024 First Term/DNA Methylation/Merged_samples/three_of_them_separately")
baseDir <- "/Users/egeezen/Desktop/PhD/2024 First Term/DNA Methylation/Merged_samples/three_of_them_separately"

#skip to reading final_m_values if you have it!

#read the files from the working directory and create a "patients" reference dataframe
list.files(baseDir)
patients_1<- read_csv("SampleSheet.csv", col_types=cols(.default = "c")) #it reads as a tibble, change the datastructure to df
patients_1$Basename <- file.path(baseDir, patients_1$slide, 
                               paste0(patients_1$slide, "_",
                                      patients_1$array))

#read the RG set
RGset_EPICv2 <- read.metharray.exp(targets = patients_1)
pd <- pData(RGset_EPICv2)
names <- pd$Patient_Name
groups <- pd$Patient_Group
pd

#put name to the RG set
patients_1$ID <- paste(patients_1$Final_Group,patients_1$Patient_Name,sep=".")
sampleNames(RGset_EPICv2) <- patients_1$ID
RGset_450k

#if we want to merge
RGset <- combineArrays(RGset_EPICv1, RGset_450k, outType = "IlluminaHumanMethylation450k")
#we have to add annotations to the RGset object, both manifest(the chip information) and annotation(the genome version)
annotation(RGset_EPICv2)["array"] = "IlluminaHumanMethylationEPICv2"
annotation(RGset_EPICv2)["annotation"] = "20a1.hg38"
annEPICv2 <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
annEPICv2_df <- as.data.frame(annEPICv2)

#For EPICv1 and hg38(this is not correct probably, maybe do a lift over from 19 to 28)
annotation(RGset)["array"] = "IlluminaHumanMethylationEPIC"
#RGset@annotation = c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b5.hg38")

#for 450k methylation this is done automatically?

#the loci names are different than the epicv1 names but there is also a column in the annotation
#file that is called EPICv1_Loci, maybe that can be used to merge EPICv2&v1 samples?
manifestEPICv2 = getManifest(IlluminaHumanMethylationEPICv2manifest)

densityPlot(RGset, sampGroups=groups,pal = brewer.pal(8, "Dark2"),legend=FALSE)
legend("top", legend = levels(factor(groups)),
       text.col=brewer.pal(8,"Dark2"))
par(mar=c(3,3,3,3)) #for the alignment of the next figure

#I don't think this is working well ? The colors are always wrong
densityBeanPlot(RGset, sampGroups=patients$Patient_Group,sampNames=rev(names),pal = brewer.pal(3, "Dark2"))

#In this qc report, the bean plot is also wrong
qcReport(RGset, sampNames= patients$ID, sampGroups = patients$Patient_Group,
         pdf="QC_report_sample_EPICv1.pdf")

ShinyMethylSet <- shinySummarize(RGset)
runShinyMethyl(ShinyMethylSet)

detP <- detectionP(RGset)
failed <- detP>0.01
colMeans(failed) # Fraction of failed positions per sample
sum(rowMeans(failed)>0.5) # How many positions failed in >50% of samples?
failed <- NULL

#you can make a bar plot to look at the quality of the samples
pal <- brewer.pal(8,"Dark2")
pal <- brewer.pal(8,"Pastel2")
barplot(colMeans(detP), col=pal[factor(groups)], las=2,
        cex.names=1, cex.axis = 0.8,
        ylab="Mean detection p-values",names.arg=names,ylim=c(0,0.01)) #can add ylim here
abline(h=0.01,col="red") #this line is for 0.01 value
legend("top", legend=levels(factor(groups)), fill=pal,
       bty="n")
dev.off()

#getting rid of the low quality ones, if they exist in any of the samples
absent <- apply(detP, 1, function(x) sum(x>0.01))
RGset_filtered <- RGset[absent==0,]

#getting the methylation data analysis using single sample Noob method
#this method is specifically very important if you want to combine different sets of data!!
#single sample method should not do across samples normalization, so you can combine EPICv2 m_values with the EPICv1 + 450k m_values for example 
ssNoob_filtered <- preprocessNoob(RGset_filtered, dyeCorr = TRUE,dyeMethod="single")

#after normalization we can look at the density plot
densityPlot(getBeta(ssNoob_filtered), sampGroups=groups,main="Density Plot of Beta Values", legend=FALSE)
legend("top", legend = levels(factor(groups)),
       text.col=brewer.pal(8,"Dark2"), bty="n",text.font=2)
dev.off()

#quick quality control if you want
qc <- getQC(ssNoob_filtered)
ssNoob_filtered <- addQC(ssNoob_filtered, qc = qc)
plotQC(qc)

#getting the methylated and unmethylated sites
methylated_sites <- getMeth(ssNoob_filtered)
unmethylated_sites <- getUnmeth(ssNoob_filtered)

#this part is important, you can do nearly everything with the beta values and m_values from now on.
#I will only use m_values or b_values from now, will not use any of the objects
beta_values_Epicv2 <- getBeta(ssNoob_filtered) #I dont know the offset here
beta_values_illumina_Epicv2 <- getBeta(ssNoob_filtered,type="Illumina") #this is different than the default, a certain offset (0.100?)
beta_values_illumina_Epicv2 <- as.data.frame(beta_values_illumina_Epicv2)

#two ways to calculate m_values one is from the function, one is taking the logit of beta values
#I still have questions about this part...
m_values <- getM(ssNoob_filtered) # if put type = "", computes M-values as the logarithm of Meth/Unmeth
#or calculate the logit of beta_values_illumina
beta_to_m <- function(beta_values) {
  m_values <- log2(beta_values / (1 - beta_values))
  return(m_values)
}

m_values_illumina <- as.data.frame(lapply(beta_values_illumina,beta_to_m))
rownames(m_values_illumina) <- rownames(beta_values_illumina)

#maybe we can get back to here by saving the csv file.
write.table(m_values_illumina_e2, file="three_of_them_separately/m_values_illumina_epicv2.csv", sep=",", row.names=TRUE)
m_values <- read.csv("m_values_EPICv1_and_450k.csv") #here use the read.csv function not read_csv..

#you dont have to save beta values you can use this function to obtain beta values
m_to_beta <- function(m_values) {
  2^m_values / (2^m_values + 1)
}

beta_values_x <- as.data.frame(lapply(m_values, m_to_beta))
rownames(beta_values) <- rownames(m_values)

#if you want to combine EPICv2 data with the 450k data
#merge the annotation of annEPICv2 with the m_values gathered from EPICv2 analysis
merged_EPICv2 <- merge(annEPICv2_df, m_values, by="row.names") #somehow they don't 100% match, this is problem from Illumina's side ?
#only get the ones that have a corresponding methyl450loci 
merged_EPICv2 <- merged_EPICv2[(merged_EPICv2$Methyl450_Loci != ""),] #now reduced to 390k from 900k
#we can look at the duplicated rows, because there are duplicates in EPICv2 for the same loci
duplicated_rows <- merged_EPICv2[duplicated(merged_EPICv2$Methyl450_Loci) | duplicated(merged_EPICv2$Methyl450_Loci, fromLast = TRUE), ]
#aggregating and taking the averages of m_values of the replicates
merged_EPICv2_avg <- aggregate(. ~ Methyl450_Loci, 
                               data = merged_EPICv2[, c(40,44:53)],
                               FUN = mean)
rownames(merged_EPICv2_avg) <- merged_EPICv2_avg$Methyl450_Loci
merged_EPICv2_avg$Methyl450_Loci <- NULL
final_m_values <- merge(beta_values_illumina, merged_EPICv2_avg,by="row.names")
write.table(m_values_illumina, file="m_values_illumina_all_merged.csv", sep=",", row.names=TRUE)
write.table(final_b_values, file="beta_values_illumina_all_merged.csv", sep=",", row.names=TRUE)

#------------------------------------#
# You can just start here if you have the final m values! (or beta values, you can calculate the m values fro it easily)
beta_values_illumina <- read.csv("~/Desktop/PhD/2024 First Term/DNA Methylation/Merged_samples/three_of_them_separately/beta_values_illumina_all_merged.csv")
rownames(beta_values_illumina) <- beta_values_illumina$Row.names
beta_values_illumina$Row.names <- NULL
m_values_illumina <- as.data.frame(beta_to_m(as.matrix(beta_values_illumina)))
rownames(m_values_illumina) <- rownames(beta_values_illumina)

#also read the patients reference data
patients <- read.csv("~/Desktop/PhD/2024 First Term/DNA Methylation/Merged_samples/three_of_them_separately/patients_new.csv")

#you can remove X and Y chromosomes from the data
#hg38(the lifted over version from hg18) for 450k :) #because our data is now reduced to 450k form, we can use this 
ann450k_hg38 <- read.delim("~/Desktop/PhD/2024 First Term/DNA Methylation/450k_samples/HM450.hg38.manifest.gencode.v36.tsv", header = TRUE)
row.names(ann450k_hg38) <- ann450k_hg38$probeID
#we can pick the X and Y chromosomes from the annotation file
keep <- !(rownames(m_values_illumina) %in% ann450k_hg38$probeID[ann450k_hg38$CpG_chrm %in% c("chrX","chrY")])
#now we can filter our m_values by removing sex chromosome loci, this will reduce variabilty significantly
m_values_filtered  <- m_values_illumina[keep,]
#if you have an outlier, get rid of it!:
m_values_filtered_wo_outlier <- m_values_filtered[,-43] #ra.shk482 is an outlier in this case
patients_wo_outlier <- patients[-43,]

#lets look at PCA, plotMDS is from limma package, mdsPlot is from minfi, the labels is optional
#pal
pal <- c("#B5C8A4","#7FB3D5","#E2A37E",  "#E6C0E2", "#A3D6C1","#E8B2A3", "#D6D1A1")
#like previously said, we only need m_values_from now on
plotMDS(m_values_illumina, top=10000, gene.selection="common",cex= 1.0, pch = 16,
        col=pal[factor(patients$sex)]) #pch is for dots, it is ignored if you give labels= patients$Patient_Name
legend("top", legend=levels(factor(patients$sex)), text.col=pal,
       bg="white", text.font=2,cex=0.8,bty="n")
dev.off()


#examining the higher dimensions of variance
par(mfrow=c(1,3))

plotMDS(m_values_filtered, top=10000, gene.selection="common",labels=patients$Patient_Name,
        col=pal[factor(patients$location)],dim=c(1,3))
legend("bottom",legend = levels(factor(patients$location)), text.col = pal,
       bg = "white", text.font = 2, bty = "n")
plotMDS(m_values_filtered, top=10000, gene.selection="common",labels=patients$Patient_Name,
        col=pal[factor(patients$location)],dim=c(2,3))
legend("bottom", legend=levels(factor(patients$location)), text.col=pal,
       bg="white", text.font=2,bty="n")
plotMDS(m_values_filtered, top=10000, gene.selection="common",labels=patients$Patient_Name,
        col=pal[factor(patients$location)],dim=c(3,4))
legend("bottomleft", legend=levels(factor(patients$location)), text.col=pal,
       bg="white", text.font=2,bty="n")

dev.off()

#hierarchical clustering heat map, in default uses euclidian distance
m_values_fixed <- m_values_filtered_wo_outlier
colnames(m_values_fixed) <- sub(".*\\.", "", colnames(m_values_fixed)) #getting rid of OA,RA, etc.
vsd_cor <- cor(m_values_fixed)
vsd_cor <- as.data.frame(vsd_cor)
patients_for_pheatmap <- as.data.frame(patients_wo_outlier)[,c(6,8),drop=FALSE] #somehow "patients" is a tibble not a dataframe
rownames(patients_for_pheatmap) <- patients_wo_outlier$Patient_Name
pheatmap(vsd_cor,annotation_col= patients_for_pheatmap,
         main = "Hierarchical Heatmap of Patient Samples", filename ="after_x_y_removal/pheatmap_corr_wo_outlier.jpeg",width= 10, height = 8)
graphics.off()

cpgs <- graph_vector[9:11]
par(mfrow=c(2,2)) #how many plots in one figure
par(mar=c(4,5,1,2)) #margins
plotCpg(final_beta_values, cpg=cpgs, pheno=patients_wo_outlier$Final_Group, type="categorical",ylab="Beta Values",ylim=c(0,1))

#to find differentially methylated positions from methyl object (F-test for categorical)
#we cant specify the specific comparisons, it compares all three of them, that is why limma is better
dmp <- dmpFinder(ssNoob_filtered, pheno=groups, type="categorical") 
sum(dmp$qval < 0.05, na.rm=TRUE) #how many sites are differentially methylated (0?)

# finding differentially methylated positions in limma 
# Healthy vs RA and the others ------
# this is the factor of interest
cellType <- factor(patients_wo_outlier$Final_Group)
# this is the joint location effect (we don't have biological replicates.. is this needed?)
jointLocation <- factor(patients_wo_outlier$location)
# creating a design matrix based on factor of interest
design <- model.matrix(~0 + cellType, data=patients_wo_outlier)
colnames(design) <- c(levels(cellType))

# or use it like this (doesn't give the same result, look further into this)
design <- model.matrix(~0 + cellType + jointLocation, data=patients_wo_outlier)
colnames(design) <- c("Healthy","OA","RA","VeRA","finger","knee","PIP","shoulder","thumb","wrist")

# if you used the first design matrix use this: remove the batch effect
adjusted_m_values <- removeBatchEffect(m_values_filtered_wo_outlier, batch=jointLocation,design=design)
# save the adjusted beta values while we are at it 
adjusted_beta_values <- 2^adjusted_m_values / (2^adjusted_m_values + 1)
write.table(adjusted_beta_values, file ="batch_effect_removed_beta_values.csv")

# fit the model
fit <- lmFit(adjusted_m_values, design) #if you used the first design matrix use the m_values_filtered_wo_outlier
# be careful about the x vs y, the controls are x most of the time
# but the logFC will be x to y difference not y to x, maybe we can do it the opposite way? 
contMatrix <- makeContrasts(VeRA-Healthy,RA-Healthy,OA-Healthy,RA-VeRA,RA-OA, levels=design)

# fit the contrasts
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)
# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit2))

# Joint location differences!: 
# getting only the RA patients
patients_RA <- patients_wo_outlier[patients_wo_outlier$Final_Group == "RA",]
patients_RA$final_location <- c("knee","shoulder","hand","shoulder","shoulder","hand","shoulder","hand","knee","knee","hand")
cellType <- factor(patients_RA$final_location)
#this is the individual effect (we don't have biological replicates.. is this needed?)
individual <- factor(patients_RA$Patient_Name)
#creating a design matrix
design <- model.matrix(~0+cellType, data=patients_RA)
colnames(design) <- c(levels(cellType))
# getting the specific cols from m_values
m_values_RA <- m_values_filtered_wo_outlier[colnames(m_values_filtered_wo_outlier) %in% patients_RA$ID ]
fit <- lmFit(m_values_RA, design)
#be careful about the x vs y, the controls are x most of the time
#but the logFC will be x to y difference not y to x, maybe we can do it the opposite way? 
contMatrix <- makeContrasts(hand-knee,shoulder-knee,hand-shoulder,levels=design)

# fit the contrasts
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)
# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit2))

#-----------------------------------------#

#then getting the annotations
#if you want to do for 450k and EPICv1 maybe we can just get the ones from EPICv2 annotation because it is hg38
cg22779972 <- as.data.frame(annEPICv2)
non_unique_names <- annEPICv2_df$Methyl450_Loci[duplicated(annEPICv2_df$Methyl450_Loci)]
annEPICv2_df_uni <- annEPICv2_df[!annEPICv2_df$Methyl450_Loci %in% non_unique_names, ]
rownames(annEPICv2_df_uni) <- annEPICv2_df_uni$Methyl450_Loci
#half of it is not even unique - there are duplicates and also some portion does not even have a Methyl450 counterpart

get_annots <- function(df,annot) {
  # Set the remaining unique values as row names
  df <- as.data.frame(merge(annot, df, by = "row.names",all.y=TRUE))
  df <- df[order(df$adj.P.Val), ]
  #writing the table as a csv file
  rownames(df) <- df$Row.names
  df <- df[,-1]
  return(df)
}

#getting the top genes (there is a problem with genelist, manifest file doesn't match one hundred percent)
DMPs_resolv_vs_RA <- topTable(fit2,  num=10, coef=1 ,sort.by="P")
DMPs_resolv_vs_RA <- get_annots(DMPs_resolv_vs_RA,annEPICv2_df_uni)
write.table(DMPs_resolv_vs_RA, file="DMPs_resolv_vs_RA.csv",sep=",",row.names = FALSE)

DMPs_RA_vs_estRA <- topTable(fit2, num=10,coef=6, sort.by="P")
DMPs_RA_vs_estRA <- get_annots(DMPs_RA_vs_estRA,annEPICv2_df_uni)

#---------------------------------------#

#or maybe just get the hg38(the lifted over version from hg18) for 450k :)
ann450k_hg38 <- read.delim("~/Desktop/PhD/2024 First Term/DNA Methylation/450k_samples/HM450.hg38.manifest.gencode.v36.tsv", header = TRUE)
row.names(ann450k_hg38) <- ann450k_hg38$probeID
#or if you use only EPIC samples you can use this as well:
annEPIC_hg38 <- read.delim("EPIC.hg38.manifest.gencode.v36.tsv", header = TRUE)
row.names(annEPIC_hg38) <- annEPIC_hg38$probeID

DMPs_Healthy_vs_VeRA <- topTable(fit2, num=Inf,coef=1, sort.by="P")
DMPs_Healthy_vs_VeRA <- get_annots(DMPs_Healthy_vs_VeRA,ann450k_hg38)
write.table(DMPs_Healthy_vs_VeRA, file="DMPs_Healthy_vs_VeRA.csv",sep=",",row.names = FALSE)

DMPs_Healthy_vs_RA <- topTable(fit2, num=Inf,coef=2, sort.by="P")
DMPs_Healthy_vs_RA <- get_annots(DMPs_Healthy_vs_RA,ann450k_hg38)
write.table(DMPs_Healthy_vs_RA, file="DMPs_Healthy_vs_RA.csv",sep=",",row.names = FALSE)

DMPs_Healthy_vs_OA <- topTable(fit2, num = Inf, coef=3, sort.by="P")
DMPs_Healthy_vs_OA <- get_annots(DMPs_Healthy_vs_OA,ann450k_hg38)
write.table(DMPs_Healthy_vs_OA, file="DMPs_Healthy_vs_OA.csv",sep=",",row.names = FALSE)

DMPs_VeRA_vs_RA <- topTable(fit2, num = Inf, coef=4, sort.by="P")
DMPs_VeRA_vs_RA <- get_annots(DMPs_VeRA_vs_RA,ann450k_hg38)
write.table(DMPs_VeRA_vs_RA, file="DMPs_VeRA_vs_RA.csv",sep=",",row.names = FALSE)

# for joint specific:
DMPs_shoulder_vs_hand <- topTable(fit2, num=Inf,coef=3, sort.by="P")
DMPs_shoulder_vs_hand <- get_annots(DMPs_shoulder_vs_hand,ann450k_hg38)
write.table(DMPs_shoulder_vs_hand, file="after_x_y_removal/joint_location/DMPs_shoulder_vs_hand.csv",sep=",",row.names = FALSE)
# etc
#---------------------------------------#

#look at the significant loci with graphs (we can find loci associated with specific genes)
par(mfrow=c(2,2))
sapply(DMPs_knee_vs_shoulder$probeID[1:4], function(cpg){
  plotCpg(m_values_RA, cpg=cpg, pheno=patients_RA$final_location, ylab = "Beta values")
})

dev.off()

genomic_methly_set <- mapToGenome(ssNoob_filtered)
genomic_ratio_set <- ratioConvert(genomic_methly_set, what = "both")

##you can remove the the probes with SNPs 
#genomic_ratio_set_filtered <- dropLociWithSnps(genomic_ratio_set)

# Create Venn diagram
set_1 <- rownames(DMPs_Healthy_vs_RA_sig)
set_2 <- rownames(DMPs_Healthy_vs_RA_sig_2nd)

library(VennDiagram)
venn.plot <- venn.diagram(
  x = list(set1 = set_1, set2 = set_2),
  category.names = c("Healthy vs RA from \n 1st analysis", "Healthy vs RA from \n 2nd analysis"), filename=NULL,
    
  # Circles
  lwd = 1,
  lty = 'blank',
  fill = c("#FFCBCB","#3C5B6F"),        cex = 1.0,
  fontface = "bold",
  fontfamily = "sans",        cat.cex = 0.9,
  cat.fontface = "bold",
  cat.default.pos = "outer",
)

# Plot Venn diagram
grid.draw(venn.plot)

#heatmap for differentially methylated positions -> unique to RA
#the first requirement for significance is padj < 0.05
DMPs_Healthy_vs_RA_sig <- subset(DMPs_Healthy_vs_RA, DMPs_Healthy_vs_RA$adj.P.Val <= 0.05 )
DMPs_Healthy_vs_OA_sig <- subset(DMPs_Healthy_vs_OA, DMPs_Healthy_vs_OA$adj.P.Val <= 0.05 )

#filtering so that the delta beta is higher than 0.1, which is 10 percent increase in methylation level
#you can either pick the RA
check_beta <- final_beta_values[rownames(final_beta_values) %in% DMPs_Healthy_vs_RA_sig$probeID,]
#getting healthy and ra from the beta values dataframe
healthy_data_ra <- check_beta[, grepl("^Healthy", colnames(check_beta))]
ra_data <- check_beta[, grepl("^RA", colnames(check_beta))]
# Calculate the averages of healthy and ra
avg_ra <- rowMeans(as.data.frame(ra_data))
avg_healthy_ra <- rowMeans(as.data.frame(healthy_data_ra))
#abs_difference should be higher than 0.2
abs_diff <- avg_ra - avg_healthy_ra
sig_loci <- names(abs_diff[abs_diff > 0.1 | abs_diff < -0.1])

#also do the calculations for OA as well!
check_beta <- final_beta_values[rownames(final_beta_values) %in% DMPs_Healthy_vs_OA_sig$probeID,]
oa_data <- check_beta[, grepl("^OA", colnames(check_beta))]
healthy_data_oa <- check_beta[, grepl("^Healthy", colnames(check_beta))]
avg_oa <- rowMeans(as.data.frame(oa_data))
avg_healthy_oa <- rowMeans(as.data.frame(healthy_data_oa))
abs_diff_1 <- avg_oa - avg_healthy_oa
sig_loci_oa <- names(abs_diff_1[abs_diff_1 > 0.1 | abs_diff_1 < -0.1])

#unique to RA sig loci
sig_loci_unique_to_ra <- setdiff(sig_loci, sig_loci_oa)

#save the sig_loci as dataframe
write.table(DMPs_Healthy_vs_RA_sig[DMPs_Healthy_vs_RA_sig$probeID %in% sig_loci,], "DMPs_Healthy_vs_all_RA_with_deltaB.csv",sep=",")
write.table(DMPs_Healthy_vs_RA_sig[DMPs_Healthy_vs_RA_sig$probeID %in% sig_loci_unique_to_ra,], "DMPs_Healthy_vs_unique_to_RA_with_deltaB.csv",sep=",")

#Heatmap drawing starts here
b_values_fixed <- final_beta_values[row.names(final_beta_values) %in% sig_loci_unique_to_ra,] #or unique genes
b_values_fixed <- b_values_fixed[,c(1:9,11,12,35:38,20:25,10,26:34,43,13:19,39:42)] #order according to patient type
#b_values_fixed <- b_values_fixed[,c(1:13,16,20:23,33,34,36,42,24,25,14,15,18,19,26,28,29,31,27,32,35,37:39,17,40,41,30,43)] #order according to joint location
colnames(b_values_fixed) <- sub(".*\\.", "", colnames(b_values_fixed)) #getting rid of OA,RA, etc.
patients_for_pheatmap <- as.data.frame(patients_wo_outlier)[,c(8,6,11),drop=FALSE] #somehow "patients" is a tibble not a dataframe
rownames(patients_for_pheatmap) <- patients_wo_outlier$Patient_Name
patients_for_pheatmap <- patients_for_pheatmap[c(1:9,11,12,35:38,20:25,10,26:34,43,13:19,39:42),,drop=FALSE]
#patients_for_pheatmap <- patients_for_pheatmap[c(1:13,16,20:23,33,34,36,42,24,25,14,15,18,19,26,28,29,31,27,32,35,37:39,17,40,41,30,43),,drop=FALSE]

#better color palette if we want 
custom_colors <- colorRampPalette(c("blue", "white", "red"))(20)

pheatmap(b_values_fixed, cluster_rows=TRUE, show_rownames=FALSE,scale = "row", color=custom_colors,
         cluster_cols=FALSE, annotation_col=patients_for_pheatmap,file="dmrs_heatmap_unique_to_RA_sig_patient_type.pdf",width=10,height=10,
         main = "(Unique to) RA vs Healthy Differentially Methylated Locations", treeheight_row = 0)


# we will pick the ones that are around the TSS regions by -500 and +100 bps 
# function to genes based on the distance to TSS
filter_genes_around_TSS <- function(distToTSS, genesUniq) {
  distances <- as.numeric(unlist(strsplit(distToTSS, ";")))
  genes <- unlist(strsplit(genesUniq, ";"))
  #filter distances
  genes_in_range <- genes[distances >= -500 & distances <= 100]
  #return the unique genes
  paste(unique(genes_in_range), collapse = ";")
}

# applying the function to each row to create the new column
DMPs_Healthy_vs_unique_to_RA_with_deltaB$genes_around_TSS <- mapply(filter_genes_around_TSS, 
                                                                    DMPs_Healthy_vs_unique_to_RA_with_deltaB$distToTSS, 
                                                                    DMPs_Healthy_vs_unique_to_RA_with_deltaB$genesUniq)
# then we can filter out the rows that have no data (NA or "")
DMPs_Healthy_vs_unique_to_RA_TSS <- DMPs_Healthy_vs_unique_to_RA_with_deltaB[DMPs_Healthy_vs_unique_to_RA_with_deltaB$genes_around_TSS != "",]
DMPs_Healthy_vs_unique_to_RA_TSS <- DMPs_Healthy_vs_unique_to_RA_TSS[DMPs_Healthy_vs_unique_to_RA_TSS$genes_around_TSS != "NA",]
# removing "NA" strings from the genes_around_TSS column
DMPs_Healthy_vs_unique_to_RA_TSS$genes_around_TSS <- gsub(";NA", "", DMPs_Healthy_vs_unique_to_RA_TSS$genes_around_TSS)

# we can split the singular genes from the dataframe:
DMPs_split <- DMPs_Healthy_vs_unique_to_RA_TSS %>%
  separate_rows(genes_around_TSS, sep = ";")
# replacing the genesUniq with the individual genes from genes_around_TSS
DMPs_split$genesUniq <- DMPs_split$genes_around_TSS
DMPs_split$genes_around_TSS <- NULL
# after putting delta B values, we will take average of the ones that have the same gene, multiple loci
averaged_data <- DMPs_split %>%
  group_by(genesUniq) %>%
  summarise(mean_deltaB = mean(deltaB, na.rm = TRUE))  %>%
  arrange(desc(mean_deltaB))

averaged_data <- as.data.frame(averaged_data)
colnames(averaged_data) <- c("Genes","Delta B")
rownames(averaged_data) <- averaged_data$Genes
averaged_data$Genes <- NULL
averaged_data$`Delta B` <- averaged_data$`Delta B`*100

write.csv(averaged_data, "DMPs_Healthy_vs_unique_to_RA_percentage.csv")

# now we can draw a heatmap
b_values_fixed <- final_beta_values[row.names(final_beta_values) %in% DMPs_split$probeID,] #or unique genes
b_values_fixed <- b_values_fixed[,c(1:9,11,12,35:38,20:25,10,26:34,43,13:19,39:42)] #order according to patient type
#b_values_fixed <- b_values_fixed[,c(1:13,16,20:23,33,34,36,42,24,25,14,15,18,19,26,28,29,31,27,32,35,37:39,17,40,41,30,43)] #order according to joint location
colnames(b_values_fixed) <- sub(".*\\.", "", colnames(b_values_fixed)) #getting rid of OA,RA, etc.
patients_for_pheatmap <- as.data.frame(patients_wo_outlier)[,c(8,6,11),drop=FALSE] #somehow "patients" is a tibble not a dataframe
rownames(patients_for_pheatmap) <- patients_wo_outlier$Patient_Name
patients_for_pheatmap <- patients_for_pheatmap[c(1:9,11,12,35:38,20:25,10,26:34,43,13:19,39:42),,drop=FALSE]
#patients_for_pheatmap <- patients_for_pheatmap[c(1:13,16,20:23,33,34,36,42,24,25,14,15,18,19,26,28,29,31,27,32,35,37:39,17,40,41,30,43),,drop=FALSE]

#better color palette if we want 
custom_colors <- colorRampPalette(c("blue", "white", "red"))(20)

pheatmap(b_values_fixed, cluster_rows=TRUE, show_rownames=FALSE,scale = "row", color=custom_colors,
         cluster_cols=FALSE, annotation_col=patients_for_pheatmap,file="dmrs_promoters_healthy_vs_RA_unique_to_RA.pdf",width=10,height=10,
         main = "(Unique to) RA vs Healthy Differentially Methylated Promoters", treeheight_row = 0)

# or we can only visualize the percentage difference (average) on the methylation levels in promoters

pheatmap(averaged_data, cluster_rows=FALSE, show_rownames=FALSE,scale = "none", color=custom_colors,
         cluster_cols=FALSE,file="dmrs_promoters_healthy_vs_RA_unique_to_RA_percentage.pdf",width=1,height=10,
         main = "Per", treeheight_row = 0)

#volcano plot for but instead of logFC in y axis, delta average beta ? 
#remove x and y chromosome loci from the final beta values df
final_beta_values_filtered <- final_beta_values[rownames(final_beta_values) %in% DMPs_Healthy_vs_RA$probeID,]
#getting healthy and ra from the beta values dataframe
healthy_volcano <- final_beta_values_filtered[, grepl("^Healthy", colnames(final_beta_values_filtered))]
ra_volcano <- final_beta_values_filtered[, grepl("^RA", colnames(final_beta_values_filtered))]
oa_volcano <- final_beta_values_filtered[, grepl("^OA", colnames(final_beta_values_filtered))]
# Calculate the averages of healthy and ra, make sure that it is a dataframe
avg_all_ra <- rowMeans(as.data.frame(ra_volcano)) #rownames is a dataframe function, I don't know but sometimes it doesn't work(probably about tibbles)
avg_all_healthy <- rowMeans(as.data.frame(healthy_volcano))
avg_all_oa <- rowMeans(as.data.frame(oa_volcano))
#differences are calculated
difference_volcano <- avg_all_ra - avg_all_healthy
#now add delta b to the dmps_healthy_vs_ra
matching <- match(DMPs_Healthy_vs_RA$probeID, names(difference_volcano))
DMPs_Healthy_vs_RA$deltaB <- difference_volcano[matching]

#volcano plot drawing starts here:
#add a column to classify points as significant or not
DMPs_Healthy_vs_RA$Classification <- with(DMPs_Healthy_vs_RA, ifelse(adj.P.Val < 0.05 & deltaB > 0.1, "Hypermethylated", 
                                                                     ifelse(adj.P.Val < 0.05 & deltaB < -0.1, "Hypomethylated", 
                                                                            "Not Significant")))
#or for unique to RA, we should get sig_loci exclusive to the unique RA
DMPs_Healthy_vs_RA$Classification <- with(DMPs_Healthy_vs_RA, ifelse(adj.P.Val < 0.05 & probeID %in% sig_loci_unique_to_ra & deltaB > 0.1, "Hypermethylated", 
                                                                     ifelse(adj.P.Val < 0.05 & probeID %in% sig_loci_unique_to_ra & deltaB < 0.1, "Hypomethylated", 
                                                                            "Not Significant")))


#calculating the -log10(padj)
DMPs_Healthy_vs_RA$log_adjP <- -log10(DMPs_Healthy_vs_RA$adj.P.Val)

#selecting points for geom text label
#label_data <- DMPs_Healthy_vs_RA[c(1,13,25,42,15,3,178,713,1210,12291,11233,37706,33675,23938,36695),,drop=FALSE]
label_data <- subset(DMPs_Healthy_vs_RA, (deltaB > 0.1 | deltaB < -0.1) & Classification != "Not Significant" )
rownames(label_data) <- NULL
label_data$rownames <- rownames(label_data)
label_data <- label_data[,c(5,6,19,21,22)]

#all RA label data
label_data <- label_data[c(5579,5174,12859,12607,9228,9371,12607,4342,4532,1078,5409,12773,4667,5356,6383,1423,9069,1,25,13,41,15,58,3,876,553,153,50,9281,7698,155,533),]
#unique to RA label data
label_data <- label_data[c(5483,4960,3332,929,989,117,2581,1284,6240,3051,2156,4983,3744,1669,3458,2965,5038,3287,2033,4055,2709,4326,5104,3022,5866,1366,1624,2711,3139,587,3522),]

ggthemr('fresh')
ggplot(DMPs_Healthy_vs_RA, aes(x = deltaB, y = log_adjP)) +
  geom_point(aes(color = Classification), size = 1,alpha=0.8) + #alpha 0.5-0.8 looks "cool"
  #scale_alpha_manual(values=c(c("Hypermethylated" = 1, "Hypomethylated" = 1, "Not Significant" = 0.25)),labels = c("Hypermethylated", "Hypomethylated", "Not Significant and/or\nNot Unique to RA")) + # if you want to add transparency
  scale_color_manual(values = c("Hypermethylated" = "#EE4E4E", "Hypomethylated" = "#2A629A", "Not Significant" = "lightgray"),labels = c("Hypermethylated", "Hypomethylated", "Not Significant")) +
  geom_vline(xintercept = c(-0.1,0.1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(
  title = "Volcano Plot of Differential Methylation \nRA compared to Healthy",
  x = "Delta Beta",
  y = "-log10(Adjusted P-Value)",
  fill = "Classification"
  ) +   # Increase plot margins if necessary
  ylim(0,NA) +
  geom_label_repel(
    data = label_data,
    aes(label = genesUniq),
    box.padding = unit(0.1, "lines"),
    segment.color = "#2B2A4C",
    max.overlaps= 20
  ) +
  guides(color = guide_legend(title = "Methylation Level", override.aes = list(size = 3)),alpha= FALSE)

#-----------------------------Creating volcano plot with two different comparisons ----- START

#create a column describing uniqueness of significant loci
DMPs_Healthy_vs_RA$unique_col <- "N/A"
DMPs_Healthy_vs_RA$unique_col[DMPs_Healthy_vs_RA$probeID %in% unique_to_RA$probeID] <- "Unique to RA"
DMPs_Healthy_vs_RA$unique_col[DMPs_Healthy_vs_RA$probeID %in% common$probeID] <- "Common in RA & OA"

#create a column for gene labels
subset_df <- DMPs_Healthy_vs_RA[(DMPs_Healthy_vs_RA$unique_col == "Unique to RA" & grepl("protein_coding", DMPs_Healthy_vs_RA$transcriptTypes)),]
ordered_subset <- subset_df[order(subset_df$logFC), "probeID"]
ordered_subset_reverse <- subset_df[order(subset_df$logFC,decreasing=TRUE), "probeID"]
condition <- c(head(ordered_subset, 15),head(ordered_subset_reverse,15))
DMPs_Healthy_vs_RA$delabel <- ifelse(DMPs_Healthy_vs_RA$probeID %in% condition, DMPs_Healthy_vs_RA$genesUniq, NA)

# Plot the volcano plot
theme_set(theme_classic(base_size = 11) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
              plot.title = element_text(hjust = 0.5)
            )
)
ggplot(data = DMPs_Healthy_vs_RA, aes(x = logFC, y = -log10(adj.P.Val), col= unique_col,label=delabel)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size =1) +
  scale_color_manual(values = c("#E6A4B4", "gray", "#51829B"), 
                     labels = c("Common in RA & OA", "Not Significant", "Unique to RA")) +
  labs(color = 'Unique or Not') +
  geom_text_repel(max.overlaps = Inf) # To show all labels 
  
#-----------------------Creating volcano plot with two different comparisons ----- END

#Finding Differentially Methylated Regions!!DMRcate
#We will do lift-over from hg19 to hg38. I wont use lifted-over external hg38 because the CpGs are *possibly* not updated
library(liftOver)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
#we extract the genome locations from the hg19
Human450klocs <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Locations 
#Granges object is created so that we can do DMRcate 
Human450kRhg19 <- GRanges(paste(Human450klocs$chr, Human450klocs$pos, sep=":"))
names(Human450kRhg19) <- rownames(Human450klocs)

genome(Human450kRhg19) <- "hg19"
seqlevelsStyle(Human450kRhg19) <- "UCSC"

#getting the liftOver chain
ch <- import.chain("~/Desktop/PhD/2024 First Term/DNA Methylation/Merged_samples/hg19ToHg38.over.chain")
Human450kRhg38 <- unlist(liftOver(Human450kRhg19, ch)) #this performs lifting over

#removing x and y from again because somehow some probes are retained
data.noXY <- rmSNPandCH(as.matrix(m_values_filtered_wo_outlier), rmXY=TRUE)

#creating myAnnotation
myAnnotation <- cpg.annotate(object = data.noXY, datatype = "array", what = "M",
                             analysis.type = "differential", design = design, 
                             contrasts = TRUE, cont.matrix = contMatrix,
                             coef = "RA - Healthy", arraytype = "450K",fdr=0.05)

#lifting over the myAnnotation object
myAnnotation.hg19 <- myAnnotation@ranges
retain <- names(myAnnotation.hg19) %in% names(Human450kRhg38)
myAnnotation.hg19 <- myAnnotation.hg19[retain]
myAnnotation.hg38 <- Human450kRhg38[names(myAnnotation.hg19)]
values(myAnnotation.hg38) <- values(myAnnotation.hg19)
#now we have our final annotation file!
myAnnotation <- new("CpGannotated", ranges=myAnnotation.hg38)

#DMRcate analysis
DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
results.ranges <- extractRanges(DMRs, genome = "hg38")
write.table(as.data.frame(results.ranges), file="after_x_y_removal/results.ranges_healthy_vs_OA.csv", sep=",", row.names=FALSE)
library(rtracklayer)
export(results.ranges,"after_x_y_removal/results_ranges_healthy_vs_OA.bed",format="bed",index=FALSE)

#b_values_are extracted as a bed file with 4 columns
#chr, chr start, chr end, value
#they are sorted, overlapping regions are averaged and returned
bed_file <- function(b_values_filtered,string){
df <- b_values_filtered[,grepl(string,colnames(b_values_filtered))]
df$avg <- rowMeans(df)*100
df <- merge(df, ann450k_hg38, by = "row.names", all.x = TRUE)
df <- df[,c("CpG_chrm","CpG_beg","CpG_end","avg")]
#b_values_healthy_vs_ra$name <- "."
#b_values_healthy_vs_ra$strand <- "."
#b_values_healthy_vs_ra <- b_values_healthy_vs_ra[, c(1,2,3,5,4,6)]
df <- na.omit(df)
df <- df[order(df$CpG_chrm, df$CpG_beg), ]
df <- aggregate(avg ~ CpG_chrm + CpG_beg + CpG_end, data = df, FUN = mean)
df <- df[order(df$CpG_chrm, df$CpG_beg), ] #should be sorted again after aggregation
return(df)
}
b_values_OA <- bed_file(b_values_filtered, "\\bOA\\.") #"\\bRA\\."
write.table(b_values_OA, file="after_x_y_removal/b_values_from_OA.bedGraph",sep="\t", row.names = FALSE, col.names=FALSE,quote = FALSE)

#GO Analysis
sigCpGs <- DMPs_Healthy_vs_RA$probeID[DMPs_Healthy_vs_RA$adj.P.Val<0.05]
all_loci <- DMPs_Healthy_vs_RA$probeID
#you can do this with a DMRs as well (results.ranges)
gst <- gometh(sig.cpg=sigCpGs, all.cpg=all_loci, plot.bias=TRUE)
top_pathways <- gst[order(gst$FDR), ]
top_pathways <- top_pathways[top_pathways$FDR < 0.05 & top_pathways$ONTOLOGY == "BP" ,]
top_pathways_second_col <- top_pathways[1:10,2]
top_pathways_rownames <- row.names(top_pathways)
new_df <- data.frame(rownames = top_pathways_rownames, first_column = top_pathways_second_col)
colnames(new_df) <- c("GO_Number","TERM")

new_df %>%gt() %>%
tab_header(title = "TOP 10 GO Biological Processes Pathways") %>%
tab_style(style = list(cell_fill(color = "#b2f7ef"),
cell_text(weight = "bold")),
locations = cells_body(columns = GO_Number))%>%
tab_style(style = list(cell_fill(color = "#ffefb5"),
cell_text(weight = "bold")), 
locations = cells_body(columns = TERM))

#extract the data -> the raw methylation data (methylation status with the b-values)
bedfile <- function (df, beta_values){
  matching <- match(df$probeID, names(beta_values))
  df$deltaB <- beta_values[matching]
  df <- df[,c("CpG_chrm","CpG_beg","CpG_end","deltaB")]
  df <- na.omit(df)
  df <- df[order(df$CpG_chrm, df$CpG_beg), ]
  df <- aggregate(deltaB ~ CpG_chrm + CpG_beg + CpG_end, data = df, FUN = mean)
  df <- df[order(df$CpG_chrm, df$CpG_beg), ] #should be sorted again after aggregation
  return(df)
}

extract_healthy <- bedfile(DMPs_Healthy_vs_RA, avg_all_healthy)
write.table(extract_healthy, file="b_values_from_healthy.bedGraph",sep="\t", row.names = FALSE, col.names=FALSE,quote = FALSE)
extract_ra <- bedfile(DMPs_Healthy_vs_RA, avg_all_ra)
write.table(extract_ra, file="b_values_from_ra.bedGraph",sep="\t", row.names = FALSE, col.names=FALSE,quote = FALSE)
extract_oa <- bedfile(DMPs_Healthy_vs_OA, avg_all_oa)
write.table(extract_oa, file="b_values_from_oa.bedGraph",sep="\t", row.names = FALSE, col.names=FALSE,quote = FALSE)

#extract the data -> the significant ones (methylation status with the b-values)
extract_sig_RA <- bedfile(DMPs_Healthy_vs_RA_sig[DMPs_Healthy_vs_RA_sig$probeID %in% sig_loci,], abs_diff[names(abs_diff) %in% sig_loci])
write.table(extract_sig_RA, file="RA_vs_Healthy_all_deltaB.bedGraph",sep="\t", row.names = FALSE, col.names=FALSE,quote = FALSE)

extract_sig_RA_unique <- bedfile(DMPs_Healthy_vs_RA_sig[DMPs_Healthy_vs_RA_sig$probeID %in% sig_loci_unique_to_ra,], abs_diff[names(abs_diff) %in% sig_loci_unique_to_ra])
write.table(extract_sig_RA_unique, file="RA_vs_Healthy_unique_deltaB.bedGraph",sep="\t", row.names = FALSE, col.names=FALSE,quote = FALSE)
