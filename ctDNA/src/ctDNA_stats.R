
library(dplyr)
ctdna <- readRDS("~/Repare Therapeutics/Repare-R&D - ClinicalBioinformatics/shared/RP6306/RP6306-01/ctDNA/processing/Tempus/data/2022-06-17_RP6306-01_filteredCtDNA.rds")

temp <- ctdna@vafsummary
temp <- select(temp, `Patient ID`, TimeFromC1D1, mVAF, VisitID) %>% unique()

bl <- filter(temp, VisitID == "Baseline") %>% rename( `Baseline mVAF` = mVAF) %>% select( - VisitID, - TimeFromC1D1)
temp <- filter(temp, VisitID != "Baseline") %>% select(-VisitID)

temp <- merge(temp, bl)
temp <- temp %>% arrange( `Patient ID`, TimeFromC1D1)
write.csv(temp, "2022-06-17_RP6306-01_StatsUpdateFile.csv", row.names = F)
