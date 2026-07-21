

#------------------- Set Directory------------------------------------------------
setwd("D:/4.Caz-Micro-v4/0.Meta_Data")
getwd()
#---------------------------------------------------------------------------------
#-----------------load packages----------------------------------------------------
#-----------------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(scales)
library(grid)
library(stringr)
library(tibble)
library(tidyr)
library(phyloseq)
#-----------------------------------------------------------------------------
list.files()
rm(list = ls());gc()
#dev.off()
#-------------------------Read MetaDataFiles-----------------------------------------
#--------------- -----Part 01 -------------------------------------------------------
Treatment_Data_A <-read.csv("D:/4.Caz-Micro-v4/0.Meta_Data/Treatmentcode_A.csv",check.names = F)
Treatment_Data_B <-read.csv("D:/4.Caz-Micro-v4/0.Meta_Data/Treatmentcode_B.csv",check.names = F)
Treatment_Data_C <- read.csv("D:/4.Caz-Micro-v4/0.Meta_Data/CAZ_Data_Unblinded -IndexCase.csv",
  check.names = FALSE) %>%
  select(IndexID, TotMem) %>%
  rename(CAZ_HH = IndexID) %>%
  filter(CAZ_HH != "")


Treatment_Data <- Treatment_Data_A %>%inner_join(Treatment_Data_B, by='PID')
#Treatment_Data_P <- Treatment_Data %>%group_by(TreatID) %>%mutate(Type = ifelse(n() == 1, "Unique", "Multiple"),HH_Num = n()) %>%ungroup()%>%distinct(TreatID, .keep_all = TRUE) %>%mutate(CAZ_HH = sub("^(CAZ-[^-]+)-.*$", "\\1", PID))
Treatment_Data_P <- Treatment_Data %>%
  group_by(TreatID) %>%
  mutate(HH_Count = n(),Type = ifelse(n() == 1, "Unique", "Multiple")) %>%ungroup() %>%
  distinct(TreatID, .keep_all = TRUE) %>%
  mutate(CAZ_HH = sub("^(CAZ-[^-]+)-.*$", "\\1", PID))%>%inner_join(Treatment_Data_C,by="CAZ_HH")%>%
  mutate(Total_Family_Mem = case_when(
        TotMem == 1  ~ "One",
        TotMem == 2  ~ "Two",
        TotMem == 3  ~ "Three",
        TotMem == 4  ~ "Four",
        TotMem == 5  ~ "Five",
        TotMem == 6  ~ "Six",
        TotMem == 7  ~ "Seven",
        TotMem == 8  ~ "Eight",
        TotMem == 9  ~ "Nine",
        TotMem == 10 ~ "Ten",
        TotMem == 11 ~ "Eleven",
        TotMem == 12 ~ "Twelve",
        TotMem == 13 ~ "Thirteen",
        TotMem == 14 ~ "Fourteen",
        TotMem == 15 ~ "Fifteen",
        TRUE ~ as.character(TotMem)))

write.csv(Treatment_Data_P,"D:/4.Caz-Micro-v4/0.Meta_Data/Treatment_Data_P.csv",row.names = F)
#-----------------------------------------------------------------------------------------------------
#------------------------------- Part 02 ( Merge all the MetaData) ----------------------------------------------
Caz_Meta_Data <- read.csv("D:/4.Caz-Micro-v4/0.Meta_Data/CAZ-metadata-simple.csv",check.names = F) %>%filter(!TP %in% c("NC", "EC"),SampleID != "12204-DS-0009")#%>%filter(SampleID != "12204-DS-0009")# As this sample has problem with run
dim(Caz_Meta_Data)

Merged_Meta<-  Caz_Meta_Data %>%inner_join(Treatment_Data_P,by="CAZ_HH")%>%mutate(Gender= ifelse(Sex==1,"Male","Female"))%>%
                mutate(Age_Group = case_when(
                AgeYr  <=4~ "1-4 years",
                AgeYr  > 4 ~ ">=5 years", 
                TRUE~ "Other/Unknown"))

#================ Keep in mind that while analysing age-specific data ,==================
#================ filter the unique samples as these are the true representative of actual age----

write.csv(Merged_Meta,"D:/4.Caz-Micro-v4/0.Meta_Data/Merged_CAZ_Meta.csv",row.names = F)

Merged_Meta <- read.csv("D:/4.Caz-Micro-v4/0.Meta_Data/Merged_CAZ_Meta.csv",check.names = F)
dim(Merged_Meta)
#----------------------------------Metadata Processing Done---------------------------------------



