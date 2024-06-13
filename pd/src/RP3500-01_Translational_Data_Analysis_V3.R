################################################################################
###  Historical version of the pd analysis for reference
################################################################################

library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(RColorBrewer)
library(openxlsx)
library(readxl)
library(ggrepel)

#####
theme_Publication <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#984ea3","#ef3b2c","#915d44","#a6cee3","#fb9a99","#7fc97f","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#984ea3","#ef3b2c","#915d44","#a6cee3","#fb9a99","#7fc97f","#ffff33")), ...)
  
}


GeneCol = c(ATM = "#386cb0",
            BRCA1 = "#fdb462",
            BRCA2 = "#984ea3",
            CDK12 = "#ef3b2c", 
            `Other STEP2` = "#915d44"
            
)

### read in clinical data file
load("C://Users/IanSilverman/Repare Therapeutics/Repare-R&D - Clinical-Team/RP-3500-01_Patient Status Update/Genomics Tracker/Patient Tracker/Tracker_Extract/2021-09-14/Data.Rda")
Clinical <- Data
load("C://Users/IanSilverman/Repare Therapeutics/Repare-R&D - Clinical-Team/RP-3500-01_Patient Status Update/Genomics Tracker/Patient Tracker/Tracker_Extract/2021-09-14/Genomic.Rda")

### Set Dirs
inDir <- "C://Users/IanSilverman/Repare Therapeutics/Repare-R&D - Clinical-Team/RP-3500-01_Patient Status Update/Genomics Tracker/Translational Data/"
setwd(inDir)

#### Read in Data
### Change
PairedBiopsy <- read_excel("C://Users/IanSilverman/Repare Therapeutics/Repare-R&D - Clinical-Team/Biomarkers and Diagnostics Team/RP-3500/RP-3500-01 TRESR/TRESR Clinical Analyses/Histowiz PD IHC Data/Cumulative_RP-3500-01_paired biopsy scoresheet_01Sep21.xlsx",sheet = 1)
BaselineBiopsy <- read_excel("Data/Histowiz15468_ScoreSheet.xlsx",sheet = 2)
IHC <- read_excel("Data/RP3500-01_Translational_Data.xlsx",sheet = 2)

### Fix Patient ID
IHC$`SUBJECT ID` <- gsub("^(.{4})(.*)$","\\1-\\2",gsub("-","",IHC$`SUBJECT ID`))

#### Remove Duplicates
PairedBiopsy <- PairedBiopsy[!duplicated(PairedBiopsy[,c(1:14)]),]
PairedBiopsy[PairedBiopsy$Stain == "p-histone (??H2AX)","Stain"] <- "pH2AX"

### Remove 18 and 23
# Remove <- c("1001-0018","1005-0023")
# PairedBiopsy <- PairedBiopsy[!(PairedBiopsy$SUBJID %in% Remove),]
### Remove unplanned
# PairedBiopsy <- PairedBiopsy[PairedBiopsy$VISIT != "Unplanned",]

### Remove QC Fail
PairedBiopsy <- PairedBiopsy[is.na(PairedBiopsy$`Exclude from dataset for QC reasons`),]
### Change to Baseline and C1D10
PairedBiopsy$VISIT <- gsub("Screening","Pre-Tx",PairedBiopsy$VISIT)
PairedBiopsy$VISIT <- gsub("Cycle 2 Day 10|Unplanned","On-Tx",PairedBiopsy$VISIT)


###Remove tumors with no scores
BaselineBiopsy <- BaselineBiopsy[!is.na(BaselineBiopsy$`Score 1`),]
Remove <- c("1004-0041_s")
BaselineBiopsy <- BaselineBiopsy[!(BaselineBiopsy$`Histowiz Sample ID` %in% Remove),]
BaselineBiopsy$Source <- "Baseline"
### Pull Baseline From Paired
PairedBiopsyPreTx <- PairedBiopsy[PairedBiopsy$VISIT == "Pre-Tx" & PairedBiopsy$Stain != "pKAP1",]
colnames(PairedBiopsyPreTx) <- gsub("H-Score","H Score",colnames(PairedBiopsyPreTx))
colnames(PairedBiopsyPreTx) <- gsub("SUBJID","Subject ID",colnames(PairedBiopsyPreTx))
PairedBiopsyPreTx$Source <- "Paired"
BaselineBiopsyCombined <- rbind(BaselineBiopsy[,c("Subject ID","H Score","Source")],PairedBiopsyPreTx[,c("Subject ID","H Score","Source")])
write.xlsx(BaselineBiopsyCombined,"Analysis/RP3500-01_BaselineBiopsy.xlsx")


BaselineBiopsyCombined<- merge(Genomic,BaselineBiopsyCombined,by.x="Patient ID",by.y="Subject ID")
write.xlsx(BaselineBiopsyCombined,"Analysis/RP3500-01_BaselineBiopsy.xlsx")

BaselineBiopsyCombined[is.na(BaselineBiopsyCombined$`Enrollment Gene`),"Enrollment Gene"] <- "Missing"
BaselineBiopsyCombined[is.na(BaselineBiopsyCombined$`Type of loss`),"Type of loss"] <- "missing"

BaselineBiopsyCombined$GeneGroup <- BaselineBiopsyCombined$`Enrollment Gene`
BaselineBiopsyCombined[grep("BRCA1|BRCA2|ATM|CDK12",BaselineBiopsyCombined$GeneGroup,invert=T),"GeneGroup"] <- "Other STEP2"
BaselineBiopsyCombined$GeneGroup2 <- BaselineBiopsyCombined$GeneGroup
BaselineBiopsyCombined[grep("BRCA1|BRCA2",BaselineBiopsyCombined$GeneGroup,),"GeneGroup2"] <- "BRCA1/2"
BaselineBiopsyCombined[grep("CDK12|Other",BaselineBiopsyCombined$GeneGroup,),"GeneGroup2"] <- "Other"


p1 <- ggplot(BaselineBiopsyCombined, aes(x=`Enrollment Gene`,y=`H Score`)) + 
  geom_boxplot(aes(fill=GeneGroup),outlier.shape = NA) +
  geom_jitter(aes(shape=`Type of loss`),size = 3, width = 0.2,height = 0) +
  ylab("Baseline ??H2AX H-Score") +
  theme_Publication() +   theme(axis.text.x=element_text(angle=45,hjust=1)) +
  scale_fill_Publication() + scale_colour_Publication() +
  facet_grid(~GeneGroup,scales="free",space="free",switch = "y")

png(paste("Plots/BaselineBiopsy_H-Score.png",sep=""),res = 300,width=10,height=5,units = "in")
print(p1)
dev.off()

p1 <- ggplot(BaselineBiopsyCombined, aes(x=reorder(`Patient ID`,`H Score`),y=`H Score`)) + 
  geom_bar(aes(fill=Source),stat="identity",position="dodge") +
  ylab("Baseline ??H2AX H-Score") +
  theme_Publication() +   theme(axis.text.x=element_text(angle=90,hjust=1)) +
  scale_fill_Publication() + scale_colour_Publication() +
  facet_grid(~GeneGroup,scales="free",space="free",switch = "y")

png(paste("Plots/BaselineBiopsy_H-Score_All.png",sep=""),res = 300,width=10,height=5,units = "in")
print(p1)
dev.off()


### Paired Analysis
### If Score is 0 set to 1
PairedBiopsy[PairedBiopsy$`H-Score`== 0,"H-Score"] <- 1
PairedBiopsy[PairedBiopsy$`%3+`== 0,"%3+"] <- 1


###
PairedBiopsy <- PairedBiopsy %>% group_by(Stain,SUBJID) %>% mutate(`H-Score % Change` = (`H-Score`[VISIT == "On-Tx"] /`H-Score`[VISIT == "Pre-Tx"] ) -1) %>% ungroup()
PairedBiopsy[PairedBiopsy$VISIT == "Pre-Tx","H-Score % Change"] <- 0

PairedBiopsy<- merge(Clinical,PairedBiopsy,by.x="Patient ID",by.y="SUBJID")
write.xlsx(PairedBiopsy,"Analysis/RP3500-01_PairedBiopsy.xlsx")



### 
PairedBiopsyH2AX <- PairedBiopsy[PairedBiopsy$Stain != "pKAP1",]
PairedBiopsyKAP1 <- PairedBiopsy[PairedBiopsy$Stain == "pKAP1",]

PairedBiopsyH2AXTest <- dcast(PairedBiopsyH2AX,`Patient ID`~VISIT+Stain,value.var="H-Score") 
PairedBiopsyKAP1Test <- dcast(PairedBiopsyKAP1,`Patient ID`~VISIT+Stain,value.var="H-Score") 

wilcox.test(PairedBiopsyH2AXTest$`Pre-Tx`,PairedBiopsyH2AXTest$`On-Tx`,paired = T)
wilcox.test(PairedBiopsyKAP1Test$`Pre-Tx`,PairedBiopsyH2AXTest$`On-Tx`,paired = T)

t.test(PairedBiopsyH2AXTest$`Pre-Tx`,PairedBiopsyH2AXTest$`On-Tx`,paired = T)
t.test(PairedBiopsyKAP1Test$`Pre-Tx`,PairedBiopsyH2AXTest$`On-Tx`,paired = T)

### Paired LIne Plots
PairedBiopsyH2AX$VISIT <- factor(PairedBiopsyH2AX$VISIT, levels = c('Pre-Tx', 'On-Tx'))
PairedBiopsyKAP1$VISIT <- factor(PairedBiopsyKAP1$VISIT, levels = c('Pre-Tx', 'On-Tx'))


p1 <- ggplot(PairedBiopsyH2AX,  aes(x=VISIT,y=`H-Score`)) +
  geom_line(aes(col=`Enrollment Gene`,group = `Patient No`)) +
  geom_point(aes(col=`Enrollment Gene`)) +
  ylab("??H2AX H-Score") +
  xlab("") +
  theme_Publication() +
  scale_colour_Publication() +   
  scale_x_discrete(expand = c(.1, .1)) 

png(paste("Plots/H2AX_H-Score_Line.png",sep=""),res = 300,width=5,height=4,units = "in")
print(p1)
dev.off()



p1 <- ggplot(PairedBiopsyKAP1,  aes(x=VISIT,y=`H-Score`)) +
  geom_line(aes(col=`Enrollment Gene`,group = `Patient No`)) +
  geom_point(aes(col=`Enrollment Gene`)) +
  ylab("p-KAP1 H-Score") +
  xlab("") +
  theme_Publication() +
  scale_colour_Publication() + 
  scale_x_discrete(expand = c(.1, .1)) 

png(paste("Plots/KAP1_H-Score_Line.png",sep=""),res = 300,width=5,height=4,units = "in")
print(p1)
dev.off()


p1 <- ggplot(PairedBiopsyH2AX,  aes(x=VISIT,y=as.numeric(`H-Score % Change`))) +
  geom_line(aes(col=`Enrollment Gene`,group = `Patient No`))+
  geom_point(data=PairedBiopsyH2AX[PairedBiopsyH2AX$VISIT == "On-Tx",],aes(col=`Enrollment Gene`)) +
  ylab("??H2AX H-Score % Change") +
  xlab("") +
  theme_Publication() +
  scale_colour_Publication() + 
  scale_y_continuous(labels = scales::percent) +
  scale_x_discrete(expand = c(.1, .1))

png(paste("Plots/H2AX_H-Score_PC_Line.png",sep=""),res = 300,width=5,height=4,units = "in")
print(p1)
dev.off()




p1 <- ggplot(PairedBiopsyKAP1,  aes(x=VISIT,y=as.numeric(`H-Score % Change`))) +
  geom_line(aes(col=`Enrollment Gene`,group = `Patient No`))+
  geom_point(data=PairedBiopsyKAP1[PairedBiopsyKAP1$VISIT == "On-Tx",],aes(col=`Enrollment Gene`)) +
  ylab("p-KAP1 H-Score % Change") +
  xlab("") +
  theme_Publication() +
  scale_colour_Publication() + 
  scale_y_continuous(labels = scales::percent) +
  scale_x_discrete(expand = c(.1, .1))


png(paste("Plots/KAP1_H-Score_PC_Line.png",sep=""),res = 300,width=5,height=4,units = "in")
print(p1)
dev.off()


### Baseline
PairedBiopsyH2AXBaseline <- PairedBiopsyH2AX[PairedBiopsyH2AX$VISIT == "Pre-Tx",]
PairedBiopsyKAP1Baseline <- PairedBiopsyKAP1[PairedBiopsyKAP1$VISIT == "Pre-Tx",]

p1 <- ggplot(PairedBiopsyH2AXBaseline, aes(y=`Patient No`,x=`H-Score`)) + 
  geom_bar(aes(fill=`Enrollment Gene`),stat="identity",position="dodge") +
  # geom_text(aes(x=- 20, label = paste(`Dose at biopsy`,"mg"))) +
  xlab("Baseline ??H2AX H-Score") +
  ylab("") +
  theme_Publication() +
  scale_fill_Publication() + 
  scale_x_continuous(limits = c(-25,max(PairedBiopsyH2AXBaseline$`H-Score`))) +
  facet_grid(GeneGroup2~.,scales="free",space="free",switch = "y")

png(paste("Plots/H2AX_H-Score_Baseline.png",sep=""),res = 300,width=8,height=5,units = "in")
print(p1)
dev.off()


p1 <- ggplot(PairedBiopsyKAP1Baseline, aes(y=`Patient No`,x=`H-Score`)) + 
  geom_bar(aes(fill=`Enrollment Gene`),stat="identity",position="dodge") +
  geom_text(aes(x=-20, label = paste(`Dose at biopsy`,"mg"))) +
  xlab("Baseline p-KAP1 H-Score") +
  ylab("") +
  theme_Publication() +
  scale_fill_Publication() + 
  scale_x_continuous(limits = c(-25,max(PairedBiopsyKAP1Baseline$`H-Score`))) +
  facet_grid(GeneGroup2~.,scales="free",space="free",switch = "y")

png(paste("Plots/KAP1_H-Score_Baseline.png",sep=""),res = 300,width=8,height=5,units = "in")
print(p1)
dev.off()




### Paired Only
PairedBiopsyH2AX$VISIT <- factor(PairedBiopsyH2AX$VISIT, levels = c('On-Tx','Pre-Tx'))
PairedBiopsyKAP1$VISIT <- factor(PairedBiopsyKAP1$VISIT, levels = c('On-Tx','Pre-Tx'))


p1 <- ggplot(PairedBiopsyH2AX,  aes(y=`Patient No`,x=`H-Score`)) +
  geom_bar(aes(fill=VISIT,group = VISIT),stat="identity",position="dodge") +
  xlab("??H2AX H-Score") +
  geom_text(aes(x=-20, label = paste(`Dose at biopsy`,"mg"))) +
  ylab("") +
  theme_Publication() +
  scale_fill_Publication()  +
  scale_x_continuous(limits = c(-25,max(PairedBiopsyH2AX$`H-Score`))) +
  facet_grid(GeneGroup2~.,scales="free",space="free",switch = "y")

png(paste("Plots/H2AX_H-Score.png",sep=""),res = 300,width=8,height=5,units = "in")
print(p1)
dev.off()

p1 <- ggplot(PairedBiopsyKAP1,  aes(y=`Patient No`,x=`H-Score`)) +
  geom_bar(aes(fill=VISIT,group = VISIT),stat="identity",position="dodge") +
  xlab("p-KAP1 H-Score") +
  geom_text(aes(x=-20, label = paste(`Dose at biopsy`,"mg"))) +
  ylab("") +
  theme_Publication() +
  scale_fill_Publication()  +
  scale_x_continuous(limits = c(-25,max(PairedBiopsyKAP1$`H-Score`))) +
  facet_grid(GeneGroup2~.,scales="free",space="free",switch = "y")

png(paste("Plots/KAP1_H-Score.png",sep=""),res = 300,width=8,height=5,units = "in")
print(p1)
dev.off()






################### IHC analysis



ATM <- IHC %>% filter(`TEST NAME` == "ATM" & `SUBJECT ID` != "1001-0018")
ATM <- dcast(ATM,`SUBJECT ID`+ `VISIT ID` + `BODY SITE` + `COLLECTION DATE` + `TEST NAME`~`ANALYTE ID`,value.var = "ANALYTE VALUE")
ATM <- ATM[ATM$`Sample Evaluable` == "Yes",]
ATM <- merge(Clinical,ATM,by.x="Patient ID",by.y="SUBJECT ID",all.y=T)
ATM <- ATM[!is.na(ATM$`Enrollment Gene`),]
# ATM[is.na(ATM$`germline status`),"germline status"]<- "TBD"

p1 <- ggplot(ATM, aes(y=reorder(`Patient No`,as.numeric(`H-Score`) ),x=as.numeric(`H-Score`))) + 
  geom_bar(stat="identity",position="dodge") +
  geom_point(aes(x=-10,shape = `Type of loss`),size=3) +
  geom_point(aes(x=-25,shape=`germline status`),size=3) +
  geom_vline(xintercept = 0,color="red",linetype = 1,size=2) +
  geom_text(aes(label=`Alteration (S)`,x=as.numeric(`H-Score`)+ 10), hjust = 0) +
  xlab("ATM H-Score") +   ylab("") +
  theme_Publication() +
  scale_fill_Publication() + 
  scale_x_continuous(limits = c(-25,max(as.numeric(ATM$`H-Score`)) * 1.3)) +
  scale_shape_manual(breaks=c("biallelic","CHIP","monoallelic","germline","somatic"),
                     values=c(15,14,0,17,2))

png(paste("Plots/ATM_H-Score.png",sep=""),res = 300,width=7,height=7,units = "in")
print(p1)
dev.off()
