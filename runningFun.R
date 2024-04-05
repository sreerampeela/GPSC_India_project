## testing new script
source("function_script.R")
library(ggplot2)

# filters to use
commonGenes <- c("cps4A", "cps4B", "hysA", "lytA", "lytC", "nanA", "pavA",
               "pce", "ply", "psaA")
commonGPSC <- c("1", "10", "6", "13", "23", "67", "904;9")

## REPEAT ANALYSIS
#boxA
boxA_res <- read.csv("jip_2024_new_boxA.csv", header = T)
perform_analysis(resDF = boxA_res, commonVFs = commonGenes, 
                 highGPSC = commonGPSC, repType = "boxA")

#boxB
boxB_res <- read.csv("jip_2024_new_boxB.csv", header = T)
perform_analysis(resDF = boxB_res, commonVFs = commonGenes, 
                 highGPSC = commonGPSC, repType = "boxB")

#boxC
boxC_res <- read.csv("jip_2024_new_boxC.csv", header = T)
perform_analysis(resDF = boxC_res, commonVFs = commonGenes, 
                 highGPSC = commonGPSC, repType = "boxC")

#RUP
rup_res <- read.csv("jip_2024_new_RUP.csv", header = T)
perform_analysis(resDF = rup_res, commonVFs = commonGenes, 
                 highGPSC = commonGPSC, repType = "RUP")

#SPRITE
sprite_res <- read.csv("jip_2024_new_SPRITE.csv", header = T)
perform_analysis(resDF = sprite_res, 
                 commonVFs = commonGenes, 
                 highGPSC = commonGPSC, repType = "SPRITE")


# plotting function
gen_plot(repDF = boxA_res, repType = "boxA", 
          commonVFs = commonGenes, highGPSC = commonGPSC)

gen_plot(repDF = boxB_res, repType = "boxB", 
         commonVFs = commonGenes, highGPSC = commonGPSC)

gen_plot(repDF = boxC_res, repType = "boxC", 
         commonVFs = commonGenes, highGPSC = commonGPSC)

gen_plot(repDF = rup_res, repType = "rup", 
         commonVFs = commonGenes, highGPSC = commonGPSC)

gen_plot(repDF = sprite_res, repType = "sprite", 
         commonVFs = commonGenes, highGPSC = commonGPSC)
