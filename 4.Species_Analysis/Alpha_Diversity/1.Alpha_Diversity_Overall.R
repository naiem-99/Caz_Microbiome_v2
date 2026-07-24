#------------------- Set Directory------------------------------------------------
setwd("D:/4.Caz-Micro-v4/Species")
getwd()
#-----------------load packages----------------------------------------------------
library(dplyr)
library(ggplot2)
library(scales)
library(grid)
library(stringr)
library(tibble)
library(phyloseq)
library(vegan)
library(ggpubr)
#===========================================================================
list.files()
rm(list = ls());gc()
Merged_Meta <- read.csv("D:/4.Caz-Micro-v4/0.Meta_Data/Merged_CAZ_Meta.csv",check.names = F)
dim(Merged_Meta)

Caz_Species <-read.csv("D:/4.Caz-Micro-v4/Species/species_Bacteria_matrix.csv", check.names = F) %>%
             select(-taxonomy_id, -taxonomy_lvl, -contains(".bracken_frac")) %>%
             rename_with(~ gsub("\\.bracken_num", "", .))%>%column_to_rownames("name") %>%   
             t() %>%                          # transpose
             as.data.frame() %>%              # matrix → data.frame
             rownames_to_column("SampleID") %>%
             mutate(SampleID = str_replace(SampleID,"(DS-)(\\d+)$",function(x) {
             prefix <- str_match(x, "(DS-)")[,2]
             num    <- str_match(x, "(\\d+)$")[,2]
             paste0(prefix, str_pad(num, width = 4, pad = "0"))}))

dim(Caz_Species)

#=====================================================================================
dir.create("D:/4.Caz-Micro-v4/Species/Processed_File_Species")
write.csv(Caz_Species,"D:/4.Caz-Micro-v4/Species/Processed_File_Species/Caz_Species_processed.csv",row.names = F)
#=====================================================================================
Merged_Species <- Merged_Meta%>%left_join(Caz_Species , by = "SampleID")
dim(Merged_Species)
write.csv(Merged_Species,"D:/4.Caz-Micro-v4/Species/Processed_File_Species/Caz_Species_merged.csv",row.names = F)
#------------------------------------------------------------------------------------------------------
#--------------------------------Merge the data--------------------------------------------------------
Merged_Meta <- read.csv("D:/4.Caz-Micro-v4/0.Meta_Data/Merged_CAZ_Meta.csv",check.names = F)
dim(Merged_Meta)
#-- the phyloseq require the matrix file (abundence)
Abundence_Matrix_Species <- Caz_Species %>%dplyr::filter(SampleID %in% Merged_Species$SampleID) %>%column_to_rownames("SampleID") %>%data.matrix()
dim(Abundence_Matrix_Species)

Merged_Meta_req_Species <- Merged_Meta %>% filter(SampleID!="12204-DS-0009")%>%column_to_rownames("SampleID")

# ===============================------------------------------------------------------
OTU_Species <- otu_table(Abundence_Matrix_Species, taxa_are_rows = FALSE)
META_Species <- sample_data(Merged_Meta_req_Species)
# ===============================
# MATCH SAMPLE ORDER (CRITICAL)--------------------------------------------------------
#----- ------------------option 3 ------------------------------------------------------
Abundence_Matrix_Species <- Abundence_Matrix_Species[rownames(META_Species), , drop = FALSE]

#setdiff(rownames(META_Species), rownames(Abundence_Matrix_Species))
#[1] "12204-DS-0009"

# --------------rebuild after matching----------------------------------------------------
OTU_Species <- otu_table(Abundence_Matrix_Species, taxa_are_rows = FALSE)



# -----------------------------BUILD PHYLOSEQ----------------------------------------------
# =============================== =============================== =========================
tax_Species <- data.frame(Species = colnames(Abundence_Matrix_Species))
rownames(tax_Species) <- colnames(Abundence_Matrix_Species)
tax_Species <- tax_table(as.matrix(tax_Species))
Physeq_Species<- phyloseq(OTU_Species, META_Species, tax_Species)


#---------------------------------------Filter and Save PHYLOSEQ----------------------------------------------
Physeq_Species <- Physeq_Species %>% prune_samples(sample_sums(.) > 0, .) %>% prune_taxa(taxa_sums(.) > 0, .)
saveRDS(Physeq_Species,"D:/4.Caz-Micro-v4/Species/Processed_File_Species/Physeq_Species.rds")
#No need to use this object while calculating alpha and beta div-----------------------------
#==========================================================================================
Physeq_Species <- readRDS("D:/4.Caz-Micro-v4/Species/Processed_File_Species/Physeq_Species.rds")
Physeq_Species_TC <- transform_sample_counts(Physeq_Species,function(x) x / sum(x))
saveRDS(Physeq_Species_TC,"D:/4.Caz-Micro-v4/Species/Processed_File_Species/Physeq_Species_transform_sample_counts.rds")
#--------------------------------------------------------------------------
Physeq_Species_TC <- readRDS("D:/4.Caz-Micro-v4/Species/Processed_File_Species/Physeq_Species_transform_sample_counts.rds")
#------------------------------------------------------------------------------------------
#------------------------Abundence Calculation--------------------------------------------------------
Physeq_Species_LF <- Physeq_Species_TC %>%psmelt 
saveRDS(Physeq_Species_LF,"D:/4.Caz-Micro-v4/Species/Processed_File_Species/Physeq_Species_LF.rds")
Physeq_Species_LF <- readRDS("D:/4.Caz-Micro-v4/Species/Processed_File_Species/Physeq_Species_LF.rds")
#-----------------------------------------------------------------------------------------
#------------------------ Alpha diversity from Phyloseq species level-----------------------------------------------------------------------------------------------------------------
Alpha_div_Overall <- phyloseq::estimate_richness(Physeq_Species,measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"))%>%rownames_to_column("SampleID")
#-------------------------------Process The Data-----------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
Alpha_div_Overall$SampleID <- gsub("^X", "", Alpha_div_Overall$SampleID)   # remove leading X
Alpha_div_Overall$SampleID <- gsub("\\.", "-", Alpha_div_Overall$SampleID) # replace . with -
#----------------------------Merge Metadata with Alpha Diversity------------------------------------------------
Meta_data_Alpha <- as(sample_data(Physeq_Species), "data.frame")
Alpha_merged_Meta <- merge(Alpha_div_Overall, Meta_data_Alpha, by.x = "SampleID", by.y = "row.names")
#write.csv(Alpha_merged_Meta,"D:/2.Caz_Microbiome_02/0.Raw_Data/Alpha_merged_Meta.csv",row.names = F)
write.csv(Alpha_merged_Meta,"D:/4.Caz-Micro-v4/Species/Processed_File_Species/Alpha_merged_Meta_estimate_richness.csv",row.names = F)
#-----------------------------Linear mixed model--------------------
Alpha_merged_Meta <-read.csv("D:/4.Caz-Micro-v4/Species/Processed_File_Species/Alpha_merged_Meta_estimate_richness.csv", check.names =F)

#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
summary(Alpha_merged_Meta$Observed)
summary(Alpha_merged_Meta$Chao1)
#=========================================================================
#---------------------- Plot the Data-------------------------------------
#------------------------ Richness----------------------------------------
####################################################################------
library(ggplot2)
library(ggpubr)
library(rstatix)
library(scales)
library(dplyr)ZZA
#https://rpubs.com/cly/compare2condition_independent
# Ensure correct factor order
Alpha_merged_Meta$Treatmentcode <- factor(Alpha_merged_Meta$Treatmentcode,levels = c("Placebo", "Azithromycin"))
Alpha_merged_Meta$TP <- factor(Alpha_merged_Meta$TP,levels = c("D-01", "D-04", "D-180"))
# -----------------------------
# Statistical test
# -----------------------------
#---- Normality check 
library(rstatix)
stat.test_Shannon<- Alpha_merged_Meta %>%group_by(TP)%>%
                    wilcox_test(Shannon ~ Treatmentcode) %>%
                    adjust_pvalue(method = "BH") %>%
                    add_significance()%>%
                    add_xy_position(x = "Treatmentcode")

#  -------------------------------------------------------------------------------
# Plot
# -------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------

#plot_Shannon <- ggplot(Alpha_merged_Meta,aes(x = Treatmentcode,y = Shannon ,fill = Treatmentcode)) +geom_boxplot(width = 0.55,outlier.shape = NA,color = "black",alpha = 0.85) +facet_wrap(~TP, nrow = 1) +scale_fill_manual(values = c( "Placebo" = "#00BFC4",       "Azithromycin" = "#F8766D" ))+ labs(title = "Shannon",x = NULL,y = "Shannon Index") +stat_pvalue_manual(stat.test_Shannon,label = "p.adj.signif",tip.length = 0.01,hide.ns = FALSE,size = 5) +theme_bw()+theme(legend.position = "none",strip.background = element_blank(),strip.text = element_text(face = "bold", size = 16),  axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black"), axis.title.x = element_text(size = 16, face = "bold"),axis.title.y = element_text(size = 16, face = "bold"),  plot.title = element_text(face = "bold", size = 20, hjust = 0.5),panel.spacing = unit(1.2, "lines") )

plot_Shannon <- ggplot(
                Alpha_merged_Meta,
                aes(x = Treatmentcode, 
                    y = Shannon, 
                    fill = Treatmentcode)) +
                geom_jitter(aes(fill = Treatmentcode),
                color = "black",width = 0.15,size = 1.8,
                alpha = 0.85,shape = 21,stroke = 0.4) +
                geom_boxplot(width = 0.55,outlier.shape = NA,
                color = "black",alpha = 0.5) +                         
                facet_wrap(~TP, nrow = 1) +
                scale_fill_manual(values = c("Placebo" = "#648FFF",
                "Azithromycin" = "#EE6AA7"))+stat_pvalue_manual(
                stat.test_Shannon, label = "p.adj.signif",
                tip.length = 0.01,hide.ns = FALSE,size = 5) +
                labs(title = "Shannon",x= NULL,y= "Shannon Index") +
                theme_bw(base_size = 12) +                                       
                theme(panel.spacing = unit(0.6, "lines"),
                axis.text.x  = element_text(angle = 45, hjust = 1, size = 10),  
                axis.text.y  = element_text(size = 10),
                axis.title   = element_text(size = 11),
                strip.background = element_rect(fill = "grey90", color = "black"),
                strip.text       = element_text(face = "bold", size = 11),
                plot.title       = element_text(hjust = 0.5, face = "bold", size = 12),
                legend.position  = "none",
                panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
                panel.grid.minor = element_line(color = "grey92", linewidth = 0.2))



plot_Shannon

#------------------------------------------------------------------------------
dir.create("D:/4.Caz-Micro-v4/Species/1.Species_Analysis")
library(ggpubr)
ggsave("D:/4.Caz-Micro-v4/Species/1.Species_Analysis/1.Alpha_Diversity_Overall/alpha_diversity_plot_Shannon.png", plot = plot_Shannon,width = 7.8,height = 7.8,dpi = 800)
summary(Alpha_merged_Meta$Shannon)
#_---------------------------------Simpson---------------------------------------
# Observed Statistical test
# -------------------------------------------------------------------------------
stat.test_Simpson <- Alpha_merged_Meta %>%group_by(TP)%>%
                     wilcox_test(Simpson ~ Treatmentcode) %>%
                     adjust_pvalue(method = "BH") %>%
                     add_significance()%>%
                     add_xy_position(x = "Treatmentcode")

# -------------------Plot

#---------------------------------------------------------------------------------------------
plot_Simpson <- ggplot(
                Alpha_merged_Meta,
                aes(x = Treatmentcode, 
                    y = Simpson, 
                    fill = Treatmentcode)) +
                geom_jitter(
                aes(fill = Treatmentcode),
                color = "black",width = 0.15,
                size = 1.8,alpha = 0.85,
                shape = 21,stroke = 0.4) +
                geom_boxplot(width = 0.55,
                outlier.shape = NA,
                color = "black",alpha = 0.5) +  # ← transparent so jitter shows through
                facet_wrap(~TP, nrow = 1) +
                scale_fill_manual(values = c(
                  "Placebo"      = "#648FFF",
                  "Azithromycin" = "#EE6AA7"))+
                stat_pvalue_manual(
                stat.test_Simpson,
                label      = "p.adj.signif",
                tip.length = 0.01,
                hide.ns    = FALSE,
                size       = 5) +
                labs(title = "Simpson",x = NULL,y = "Simpson Index") +
                theme_bw(base_size = 12) +                                        # ← reduced from 16
                theme(panel.spacing       = unit(0.6, "lines"),
                axis.text.x         = element_text(angle = 45, hjust = 1, size = 10),  # ← explicit size
                axis.text.y         = element_text(size = 10),
                axis.title          = element_text(size = 11),
                strip.background    = element_rect(fill = "grey90", color = "black"),
                strip.text          = element_text(face = "bold", size = 11),
                plot.title          = element_text(hjust = 0.5, face = "bold", size = 12), # ← smaller title
                legend.position     = "none",
                panel.grid.major    = element_line(color = "grey85", linewidth = 0.4),
                panel.grid.minor    = element_line(color = "grey92", linewidth = 0.2))
 

plot_Simpson
#------------------------------------------------------------------------------
ggsave("D:/4.Caz-Micro-v4/Species/1.Species_Analysis/1.Alpha_Diversity_Overall/alpha_diversity_plot_Simpson.png", plot = plot_Simpson,width = 7.8,height = 7.8,dpi = 800)
#-----------------------------------------------------------------------------------
#_---------------------------------InvSimpson---------------------------------------

# Observed Statistical test
# -----------------------------
stat.test_InvSimpson <- Alpha_merged_Meta %>%group_by(TP)%>%wilcox_test(InvSimpson ~ Treatmentcode) %>%
  adjust_pvalue(method = "BH") %>%add_significance()%>%add_xy_position(x = "Treatmentcode")

#  -------------------------------------------------------------------------------
# Plot
# --------------------------------------------------------------------------------
#---------------------------------------------------------------------------------

plot_InvSimpson <- ggplot(Alpha_merged_Meta,
                  aes(x = Treatmentcode, y = InvSimpson, fill = Treatmentcode)) +
                  geom_jitter(aes(fill = Treatmentcode),
                  color = "black",
                  width = 0.15,
                  size = 1.8,
                  alpha = 0.85,
                  shape = 21,
                  stroke = 0.4) +
                  geom_boxplot(width = 0.55,
                  outlier.shape = NA,
                  color = "black",
                  alpha = 0.5) +                         # ← transparent so jitter shows through
                  facet_wrap(~TP, nrow = 1) +
                  scale_fill_manual(values = c("Placebo"= "#648FFF","Azithromycin" = "#EE6AA7"))+
                  stat_pvalue_manual(
                  stat.test_InvSimpson,
                  label      = "p.adj.signif",
                  tip.length = 0.01,
                  hide.ns    = FALSE,
                  size       = 5) +
                 labs(title = "InvSimpson",x = NULL,y = "InvSimpson Index")  +
                 theme_bw(base_size = 12) +                                        # ← reduced from 16
                 theme(panel.spacing       = unit(0.6, "lines"),
                 axis.text.x         = element_text(angle = 45, hjust = 1, size = 10),  # ← explicit size
                 axis.text.y         = element_text(size = 10),
                 axis.title          = element_text(size = 11),
                 strip.background    = element_rect(fill = "grey90", color = "black"),
                 strip.text          = element_text(face = "bold", size = 11),
                 plot.title          = element_text(hjust = 0.5, face = "bold", size = 12), # ← smaller title
                 legend.position     = "none",
                 panel.grid.major    = element_line(color = "grey85", linewidth = 0.4),
                 panel.grid.minor    = element_line(color = "grey92", linewidth = 0.2))


plot_InvSimpson

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
ggsave("D:/4.Caz-Micro-v4/Species/1.Species_Analysis/1.Alpha_Diversity_Overall/alpha_diversity_plot_InvSimpson.png", plot = plot_InvSimpson,width = 7.8,height = 7.8,dpi = 800)

All_Alpha_D_Plot<-plot_Shannon+plot_Simpson+plot_InvSimpson 
ggsave("D:/4.Caz-Micro-v4/Species/1.Species_Analysis/1.Alpha_Diversity_Overall/alpha_diversity_plot_All_Metrics.png", plot = All_Alpha_D_Plot,width = 16,height = 8,dpi = 800)
#-------------------------------------------------The End----------------------------------------------------------------------------------------------------------------------

