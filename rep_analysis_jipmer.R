### JIPMER repeat analysis

library(tidyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(FSA)

source("./function_script.R")

commonVFs <- c("cbpD", "cps4A",	"cps4B",	"hysA",	"lytA",	"lytC",	"nanA",	"pavA",
               "pce",	"ply",	"psaA")

highGPSC <- c("1", "10", "6", "13", "23", "67", "904;9")

#SPRITE
spriteReps <- read.csv("final_cleaned_SPRITE.csv", header = T)

sprite_clean <- spriteReps %>%
  group_by(Gene) %>%
  filter(length(unique(Upstream)) > 3)
perform_analysis(resDF = sprite_clean, commonVFs = commonVFs, 
                 highGPSC = highGPSC, repType = "SPRITE")

gen_plot(repDF = sprite_clean, repType = "SPRITE", commonVFs = commonVFs, 
         highGPSC = highGPSC)

#RUP
rupReps <- read.csv("final_cleaned_RUP.csv", header = T)

rup_clean <- rupReps %>%
  group_by(Gene) %>%
  filter(length(unique(Upstream)) > 3)
perform_analysis(resDF = rup_clean, commonVFs = commonVFs, 
                 highGPSC = highGPSC, repType = "RUP")

gen_plot(repDF = rup_clean, repType = "RUP", commonVFs = commonVFs, 
         highGPSC = highGPSC)

#boxA
boxA_Reps <- read.csv("final_cleaned_boxA.csv", header = T)

boxA_clean <- boxA_Reps %>%
  group_by(Gene) %>%
  filter(length(unique(Upstream)) > 3) # atleast 3 unique values to be valid
perform_analysis(resDF = boxA_clean, commonVFs = commonVFs, 
                 highGPSC = highGPSC, repType = "boxA")

gen_plot(repDF = boxA_clean, repType = "boxA", commonVFs = commonVFs, 
         highGPSC = highGPSC)

#boxB
boxB_Reps <- read.csv("final_cleaned_boxB.csv", header = T)

boxB_clean <- boxB_Reps %>%
  group_by(Gene) %>%
  filter(length(unique(Upstream)) > 3)
perform_analysis(resDF = boxB_clean, commonVFs = commonVFs, 
                 highGPSC = highGPSC, repType = "boxB")

gen_plot(repDF = boxB_clean, repType = "boxB", commonVFs = commonVFs, 
         highGPSC = highGPSC)
#boxC
boxC_Reps <- read.csv("final_cleaned_boxC.csv", header = T)

boxC_clean <- boxC_Reps %>%
  group_by(Gene) %>%
  filter(length(unique(Upstream)) > 3)
perform_analysis(resDF = boxC_clean, commonVFs = commonVFs, 
                 highGPSC = highGPSC, repType = "boxC")

gen_plot(repDF = boxC_clean, repType = "boxC", commonVFs = commonVFs, 
         highGPSC = highGPSC)


## large df
commonVF_boxA <- boxA_clean %>%
  filter(Gene %in% commonVFs)
commonVF_boxB <- boxB_clean %>%
  filter(Gene %in% commonVFs)
commonVF_boxC <- boxC_clean %>%
  filter(Gene %in% commonVFs)
commonVF_rup <- rup_clean %>%
  filter(Gene %in% commonVFs)
commonVF_sprite <- sprite_clean %>%
  filter(Gene %in% commonVFs)

jip_reps_df <- data.frame(matrix(data=NA, ncol = 5, 
                                 nrow = sum(nrow(commonVF_boxA), nrow(commonVF_boxB), 
                                            nrow(commonVF_boxC), nrow(commonVF_rup),
                                            nrow(commonVF_sprite))))
colnames(jip_reps_df) <- c("sample", "Gene", "repeats", "GPSC" , "value")

jip_reps_df$sample <- c(commonVF_boxA$Sample,commonVF_boxB$Sample,commonVF_boxC$Sample,
                        commonVF_rup$Sample,commonVF_sprite$Sample)



jip_reps_df$Gene <- c(commonVF_boxA$Gene,commonVF_boxB$Gene,commonVF_boxC$Gene, 
                      commonVF_rup$Gene,commonVF_sprite$Gene)

jip_reps_df$repeats <- c(rep("boxA", times=nrow(commonVF_boxA)), rep("boxB", times=nrow(commonVF_boxB)),
                         rep("boxC", times=nrow(commonVF_boxC)), rep("RUP", times=nrow(commonVF_rup)),
                         rep("SPRITE", times=nrow(commonVF_sprite)))
jip_reps_df$value <- c(commonVF_boxA$Upstream,commonVF_boxB$Upstream,commonVF_boxC$Upstream,
                       commonVF_rup$Upstream,commonVF_sprite$Upstream)

jip_reps_df$GPSC <- c(commonVF_boxA$GPSC,commonVF_boxB$GPSC,commonVF_boxC$GPSC,
                       commonVF_rup$GPSC,commonVF_sprite$GPSC) 
ggplot(data = jip_reps_df, aes(x = repeats, y= value, fill = Gene)) + geom_boxplot() +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16),
        axis.title = element_text(size = 18), legend.text = element_text(size = 18),
        legend.title = element_text(size = 20))

## gpsc SPECIFIC DIFFERENCES
commonVF_boxA <- boxA_clean %>%
  filter(Gene %in% commonVFs & GPSC %in% highGPSC)
commonVF_boxB <- boxB_clean %>%
  filter(Gene %in% commonVFs & GPSC %in% highGPSC)
commonVF_boxC <- boxC_clean %>%
  filter(Gene %in% commonVFs & GPSC %in% highGPSC)
commonVF_rup <- rup_clean %>%
  filter(Gene %in% commonVFs & GPSC %in% highGPSC)
commonVF_sprite <- sprite_clean %>%
  filter(Gene %in% commonVFs & GPSC %in% highGPSC)

jip_reps_df <- data.frame(matrix(data=NA, ncol = 5, 
                                 nrow = sum(nrow(commonVF_boxA), nrow(commonVF_boxB), 
                                            nrow(commonVF_boxC), nrow(commonVF_rup),
                                            nrow(commonVF_sprite))))
colnames(jip_reps_df) <- c("sample", "Gene", "repeats", "GPSC" , "value")

jip_reps_df$sample <- c(commonVF_boxA$Sample,commonVF_boxB$Sample,commonVF_boxC$Sample,
                        commonVF_rup$Sample,commonVF_sprite$Sample)



jip_reps_df$Gene <- c(commonVF_boxA$Gene,commonVF_boxB$Gene,commonVF_boxC$Gene, 
                      commonVF_rup$Gene,commonVF_sprite$Gene)

jip_reps_df$repeats <- c(rep("boxA", times=nrow(commonVF_boxA)), rep("boxB", times=nrow(commonVF_boxB)),
                         rep("boxC", times=nrow(commonVF_boxC)), rep("RUP", times=nrow(commonVF_rup)),
                         rep("SPRITE", times=nrow(commonVF_sprite)))
jip_reps_df$value <- c(commonVF_boxA$Upstream,commonVF_boxB$Upstream,commonVF_boxC$Upstream,
                       commonVF_rup$Upstream,commonVF_sprite$Upstream)

jip_reps_df$GPSC <- c(commonVF_boxA$GPSC,commonVF_boxB$GPSC,commonVF_boxC$GPSC,
                      commonVF_rup$GPSC,commonVF_sprite$GPSC) 

jip_reps_df$GPSC <- as.factor(jip_reps_df$GPSC)
kruskall_stats_tests <- jip_reps_df %>%
  group_by(Gene) %>%
  # print(summary(value))
  summarise(pval = kruskal.test(value ~ GPSC)$p.value)
kruskall_stats_tests["adjP"] <- p.adjust(kruskall_stats_tests$pval)

dunn_test_res <- jip_reps_df %>%
  group_by(Gene) %>%
  summarise(adjpvals = dunn.test(x=value, g = GPSC, method = "bh")$P.adjust,
            comp = dunn.test(x=value, g = GPSC, method = "bh")$comparisons,
            zvals = dunn.test(x=value, g = GPSC, method = "bh")$Z)
plot(x=dunn_test_res$zvals, y = -log10(dunn_test_res$adjpvals))
abline(h=-log10(0.05), v = c(-1.96, 1.96))

write.csv(jip_reps_df, file = "JIPMER_2024_vfs_reps.csv", row.names = F)



