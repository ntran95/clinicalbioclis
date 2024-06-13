
library(dplyr)
library(readxl)

source("functions_ctDNA.R")

tracker <- "~/Repare Therapeutics/Repare-R&D - ClinicalBioinformatics/shared/RP6306/RP6306-01/clinical/20220715_RP6306-01 Patient Tracker.xlsx"
clinical_list <- readxlAllSheets(xlsxFile = tracker)

### Keep only relevant tabs
tabs_to_keep <- c("Genomics", "SNIPDx")
clinical_list <- clinical_list[ names(clinical_list) %in% tabs_to_keep ]

### 

genomics <- select(clinical_list$Genomics, `Patient ID`, `Enrollment Gene`, `Alteration Classification`, `TP53 Status (curated)`, `KRAS Status (curated)`)
snipdx <- select(clinical_list$SNIPDx, `Patient ID`, `CCNE1 CN Call`, `Type of loss`, `TP53 mutant`, `KRAS mutant`)

comb <- merge(genomics, snipdx, all = T)

### 

comb <- filter(comb, !is.na(`Enrollment Gene`) )
comb$`CCNE1 CN Call` <- ifelse(comb$`Enrollment Gene` == "CCNE1",comb$`CCNE1 CN Call`, NA )
comb$`KRAS mutated` <- ifelse(is.na(comb$`KRAS mutant`) & is.na(comb$`KRAS Status (curated)`), "No", "Yes")
comb$`TP53 mutated` <- ifelse(is.na(comb$`TP53 mutant`) & is.na(comb$`TP53 Status (curated)`), "No", "Yes")
comb <- select(comb, `Patient ID`, `Enrollment Gene`, `Alteration Classification`, `CCNE1 CN Call`, `Type of loss`, `KRAS mutated`, `TP53 mutated`)
comb <- comb %>% rename(`FBXW7/PPP2R1A type of loss` = `Type of loss`)

write.table(comb, "2022-07-15_Genomics_StatsUpdate.txt", row.names = F, quote = F, sep = "\t", na = "")


tracker <- "~/Repare Therapeutics/Repare-R&D - ClinicalBioinformatics/shared/RP6306/RP6306-02/clinical/20220715_RP6306-02 Patient Tracker.xlsx"
clinical_list <- readxlAllSheets(xlsxFile = tracker)

### Keep only relevant tabs
tabs_to_keep <- c("Genomics", "SNIPDx")
clinical_list <- clinical_list[ names(clinical_list) %in% tabs_to_keep ]

### 

genomics <- select(clinical_list$Genomics, `Patient ID`, `Enrollment Gene`, `Alteration Classification`, `TP53 Status (curated)`, `KRAS Status (curated)`)
snipdx <- select(clinical_list$SNIPDx, `Patient ID`, `CCNE1 CN Call`, `Type of loss`, `TP53 mutant`, `KRAS mutant`)

comb <- merge(genomics, snipdx, all = T)

### 

comb <- filter(comb, !is.na(`Enrollment Gene`) )
comb$`CCNE1 CN Call` <- ifelse(comb$`Enrollment Gene` == "CCNE1",comb$`CCNE1 CN Call`, NA )
comb$`KRAS mutated` <- ifelse(is.na(comb$`KRAS mutant`) & is.na(comb$`KRAS Status (curated)`), "No", "Yes")
comb$`TP53 mutated` <- ifelse(is.na(comb$`TP53 mutant`) & is.na(comb$`TP53 Status (curated)`), "No", "Yes")
comb <- select(comb, `Patient ID`, `Enrollment Gene`, `Alteration Classification`, `CCNE1 CN Call`, `Type of loss`, `KRAS mutated`, `TP53 mutated`)
comb <- comb %>% rename(`FBXW7/PPP2R1A type of loss` = `Type of loss`)

write.table(comb, "2022-07-15_RP6306-02_Genomics_StatsUpdate.txt", row.names = F, quote = F, sep = "\t", na = "")
