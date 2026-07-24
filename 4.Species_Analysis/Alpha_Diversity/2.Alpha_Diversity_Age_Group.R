#------------------- Set Directory------------------------------------------------
setwd("D:/4.Caz-Micro-v4/Species")
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
library(phyloseq)
library(ggplot2)
library(dplyr)
library(rstatix)
library(ggpubr)
#-----------------------------------------------------------------------------
list.files()
rm(list = ls());gc()
#------------------------------------------------------------------------------
Alpha_merged_Meta <-read.csv("D:/4.Caz-Micro-v4/Species/Processed_File_Species/Alpha_merged_Meta_estimate_richness.csv", check.names =F)

#----------------------- For unique Sample  Analysis ----------------------------------------
#=========================================================
# Unique Sample Subset
#=========================================================
Unique_sample_analysis_Alpha_Div <- Alpha_merged_Meta %>%filter(Type == "Unique")

# Age group ordering
Unique_sample_analysis_Alpha_Div$Age_Group <- factor(Unique_sample_analysis_Alpha_Div$Age_Group,levels = c("1-4 years", ">=5 years"))

# TP ordering
Unique_sample_analysis_Alpha_Div$TP <- factor( Unique_sample_analysis_Alpha_Div$TP, levels = c("D-01", "D-04", "D-180"))

# Treatment ordering (optional but recommended)
Unique_sample_analysis_Alpha_Div$Treatmentcode <- factor(Unique_sample_analysis_Alpha_Div$Treatmentcode,levels = c("Placebo", "Azithromycin"))


#=========================================================
#---------------- SHANNON -------------------------------
#=========================================================
stat.test_Shannon_uni <- Unique_sample_analysis_Alpha_Div %>%group_by(TP, Age_Group) %>%wilcox_test(Shannon ~ Treatmentcode) %>%
  adjust_pvalue(method = "BH") %>%add_significance() %>%add_xy_position(x = "Age_Group")


plot_Shannon_uni <- ggplot(
                    Unique_sample_analysis_Alpha_Div,
                    aes(x = Age_Group, y = Shannon, fill = Treatmentcode)) +
                    geom_boxplot(
                    width = 0.5,                              
                    outlier.shape = NA,
                    color = "black",
                    alpha = 0.5,
                    position = position_dodge(0.6)) +         
                    geom_jitter(
                    aes(fill = Treatmentcode),
                    color = "black",
                    size = 1.8,
                    alpha = 0.85,
                    shape = 21,
                    stroke = 0.4,
                    position = position_jitterdodge(
                    jitter.width = 0.1,
                    dodge.width  = 0.6)) +                 
                    facet_wrap(~TP, nrow = 1) +
                    scale_fill_manual(values = c(
                   "Placebo"      = "#836FFF",
                    "Azithromycin" = "#66CDAA")) +
                    stat_pvalue_manual(
                    stat.test_Shannon_uni,
                   label         = "p.adj.signif",
                   tip.length    = 0.01,
                   hide.ns       = FALSE,
                   size          = 4,
                   step.increase = 0.08) +
                   labs(title = "Shannon", x     = NULL,
                   y     = "Shannon Index") +
                   theme_bw(base_size = 16) +
                   theme(legend.position = "top" ,
                   panel.spacing    = unit(0.6, "lines"),
                   axis.text.x      = element_text(angle = 45, hjust = 1),
                   strip.background = element_rect(fill = "grey90", color = "black"),
                   strip.text       = element_text(face = "bold"),
                   plot.title       = element_text(hjust = 0.5, face = "bold"),
                   panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
                   panel.grid.minor = element_line(color = "grey92", linewidth = 0.2))

plot_Shannon_uni

ggsave("D:/4.Caz-Micro-v4/Species/1.Species_Analysis/2.Alpha_Diversity_Age/alpha_diversity_plot_Shannon_Age.png", plot = plot_Shannon_uni,width = 7.8,height = 7.8,dpi = 800)

#=========================================================
#---------------- SIMPSON -------------------------------
#=========================================================
stat.test_Simpson_uni <- Unique_sample_analysis_Alpha_Div %>%
                         group_by(TP, Age_Group) %>%
                         wilcox_test(Simpson ~ Treatmentcode) %>%
                         adjust_pvalue(method = "BH") %>%
                         add_significance() %>%
                         add_xy_position(x = "Age_Group")



plot_Simpson_uni <- ggplot(
                    Unique_sample_analysis_Alpha_Div,
                    aes(x = Age_Group, y = Simpson, fill = Treatmentcode)) +
                    geom_boxplot(
                    width = 0.5,                              
                    outlier.shape = NA,
                    color = "black",
                    alpha = 0.5,
                    position = position_dodge(0.6)) +         
                    geom_jitter(
                    aes(fill = Treatmentcode),
                    color = "black",
                    size = 1.8,
                    alpha = 0.85,
                    shape = 21,
                    stroke = 0.4,
                    position = position_jitterdodge(
                    jitter.width = 0.1,
                    dodge.width  = 0.6)) +                 
                    facet_wrap(~TP, nrow = 1) +
                    scale_fill_manual(values = c(
                    "Placebo"      = "#836FFF",
                     "Azithromycin" = "#66CDAA")) +
                     stat_pvalue_manual(stat.test_Simpson_uni,
                     label         = "p.adj.signif",
                     tip.length    = 0.01,
                     hide.ns       = FALSE,
                     size          = 4)+#step.increase = 0.08) +
                     labs(title = "Simpson", x = NULL, y = "Simpson Index") +
                     theme_bw(base_size = 16) +
                     theme(legend.position = "top" ,
                     panel.spacing    = unit(0.6, "lines"),
                     axis.text.x      = element_text(angle = 45, hjust = 1),
                     strip.background = element_rect(fill = "grey90", color = "black"),
                     strip.text       = element_text(face = "bold"),
                     plot.title       = element_text(hjust = 0.5, face = "bold"),
                     panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
                     panel.grid.minor = element_line(color = "grey92", linewidth = 0.2))


plot_Simpson_uni

ggsave("D:/4.Caz-Micro-v4/Species/1.Species_Analysis/2.Alpha_Diversity_Age/alpha_diversity_plot_Simpson_Age.png", plot = plot_Simpson_uni,width = 7.8,height = 7.8,dpi = 800)
#=========================================================
#---------------- INVSIMPSON -----------------------------
#=========================================================
stat.test_InvSimpson_uni <- Unique_sample_analysis_Alpha_Div %>%
                            group_by(TP, Age_Group) %>%
                            wilcox_test(InvSimpson ~ Treatmentcode) %>%
                            adjust_pvalue(method = "BH") %>%
                            add_significance() %>%
                            add_xy_position(x = "Age_Group")

plot_InvSimpson_uni <-  ggplot(
                         Unique_sample_analysis_Alpha_Div,
                         aes(x = Age_Group, y = InvSimpson, fill = Treatmentcode)) +
                         geom_boxplot(
                         width = 0.5,                              
                         outlier.shape = NA,
                         color = "black",
                         alpha = 0.5,
                         position = position_dodge(0.6)) +         
                         geom_jitter(
                         aes(fill = Treatmentcode),
                         color = "black",
                         size = 1.8,
                         alpha = 0.85,
                         shape = 21,
                         stroke = 0.4,
                         position = position_jitterdodge(
                         jitter.width = 0.1,
                         dodge.width  = 0.6)) +                  
                         facet_wrap(~TP, nrow = 1) +
                         scale_fill_manual(values = c(
                        "Placebo"      = "#836FFF",
                        "Azithromycin" = "#66CDAA")) +
                         stat_pvalue_manual(
                         stat.test_InvSimpson_uni,
                         label         = "p.adj.signif",
                         tip.length    = 0.01,
                         hide.ns       = FALSE,
                         size          = 4)+#step.increase = 0.08) +
                         labs(title = "InvSimpson", x = NULL, y = "InvSimpson Index") +
                         theme_bw(base_size = 16) +
                         theme(legend.position = "top" ,
                         panel.spacing    = unit(0.6, "lines"),
                         axis.text.x      = element_text(angle = 45, hjust = 1),
                         strip.background = element_rect(fill = "grey90", color = "black"),
                         strip.text       = element_text(face = "bold"),
                         plot.title       = element_text(hjust = 0.5, face = "bold"),
                         panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
                         panel.grid.minor = element_line(color = "grey92", linewidth = 0.2))

plot_InvSimpson_uni

ggsave("D:/4.Caz-Micro-v4/Species/1.Species_Analysis/2.Alpha_Diversity_Age/alpha_diversity_plot_InvSimpson_Age.png", plot = plot_InvSimpson_uni,width = 7.8,height = 7.8,dpi = 800)


# ─── Combined: side by side, one legend at top ────────────────────────────────
combined <- ggarrange(
  plot_Shannon_uni, plot_Simpson_uni, plot_InvSimpson_uni,
  ncol          = 3,
  nrow          = 1,
  common.legend = TRUE,
  legend        = "top"
)

combined <- annotate_figure(
  combined,
  top = text_grob(
    "",
    face = "bold", size = 17
  )
)

combined

# ─── Save individual plots ────────────────────────────────────────────────────
out <- "D:/4.Caz-Micro-v4/Species/1.Species_Analysis/2.Alpha_Diversity_Age"


# ─── Save combined ────────────────────────────────────────────────────────────
ggsave(file.path(out, "alpha_diversity_Age_combined.png"),
       combined, width = 18, height = 7, dpi = 800)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------
