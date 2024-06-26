## ICMR GPSC project data analysis - main file
library(xlsx)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(ggiraphExtra)
library(vegan)
library(picante)
library(knitr)


## reading DATA
df <- read.xlsx("allData_analysis_FINAL.xlsx", sheetName = "allData_analysis")
## look for empty rows
df <- df[1:197, ]

## demographic details
ages <- summary(df$PATIENT.AGE)
hist(df$PATIENT.AGE, ylim = c(0,60), xlab = "Patient Age (years)", 
     ylab = "Number of isolates", main = "Patient ages (N=197)")

df$SOURCE.OF.ISOLATE <- as.factor(df$SOURCE.OF.ISOLATE)
df$Year <- as.factor(df$Year)
plot(df$Year, ylim = c(0,100), xlab = "year", ylab = "Number of isolates", 
     main = "Isolates collected among study years", cex.axis = 1.25, 
     cex.names = 1.25, cex.lab = 1.5, cex.main = 2)

df$ageGrp <- cut(df$PATIENT.AGE, breaks = c(0, 5,14,46,60, 100), 
                 labels = c("<5yrs", "5-13yrs", "14-45yrs", "46-60yrs", ">60yrs"))
df_age <- data.frame(table(df$ageGrp))
barplot(df_age$Freq, names.arg = df_age$Var1, xlab = "Age group", 
        ylab = "Numbe of isolates", ylim = c(0,100), cex.axis = 1.25, 
        cex.names = 1.25, cex.lab = 1.5, 
        main = "Number of isolates among different age groups", cex.main = 2)

############ GPSC types #############
df$Strain <- as.factor(df$Strain)
df$Invasive.Non.invasive <- as.factor(df$Invasive.Non.invasive)
gpscs <- as.data.frame(table(df$Strain, df$Invasive.Non.invasive))
colnames(gpscs) <- c("GPSC", "Type_infection", "freq")
ggplot(data = gpscs, aes(x=freq, y=reorder(GPSC, freq), fill = Type_infection)) + 
  geom_bar(stat = "identity", position = "stack") + 
  labs(x = "Number of isolates", y = "GPSC") +
  theme(text = element_text(size = 20))

# GPSC diversity measures
diversity(table(df$Invasive.Non.invasive, df$Strain), index = "simpson")
diversity(table(df$Year, df$Strain), index = "simpson")
dfSmall <- df[ , c(6, 36, 80)]
df2 <- dfSmall[dfSmall$Year %in% c("2016", "2018", "2021", "2022", "2023"), ]
x <- vegdist(table(df2$Year, df2$Strain), method = "jaccard")


############## AST #################
## Analysed for only prospective isolates
##penicillin MEAN, median and IQR MIC
dfAST <- df[df$Year %in% c("2020", "2021", "2022", "2023"), ]
dfAST$Invasive.Non.invasive <- as.factor(dfAST$Invasive.Non.invasive)
dfAST$Year <- as.factor(dfAST$Year)

## penicillin MIC change
boxplot(dfAST$PnG ~ dfAST$Year, xlab = "Year", ylab = "Penicillin MIC (�g/mL)")
boxplot(dfAST$PnG ~ dfAST$Invasive.Non.invasive, xlab = "Type of infection", 
        ylab = "Penicillin MIC (�g/mL)")
meanPEN <- summarySE(data = dfAST, measurevar = "PnG", groupvars = "Year", na.rm = T)
kruskal.test(PnG ~ Year, data = dfAST)
wilcox.test(PnG ~ Invasive.Non.invasive, data = dfAST)
summary(dfAST$PnG)

##CEFTRIAXONE MIC MEAN, MEDIAN AND IQR
boxplot(dfAST$Ceft ~ dfAST$Year, xlab = "Year", ylab = "Ceftriaxone MIC (�g/mL)")
boxplot(dfAST$Ceft ~ dfAST$Invasive.Non.invasive, xlab = "Type of infection", 
        ylab = "Ceftriaxone MIC (�g/mL)")
meanCTX <- summarySE(data = dfAST, measurevar = "Ceft", groupvars = "Year", na.rm = T)
kruskal.test(Ceft ~ Year, data = dfAST)
wilcox.test(Ceft ~ Invasive.Non.invasive, data = dfAST)

## MDR
mdr <- data.frame(table(dfAST$Year, dfAST$MDR_YES))
colnames(mdr) <- c("Year", "MDR", "Freq")
mdr_pros <- mdr %>%
  filter(Year %in% c("2020", "2021", "2022", "2023"))
ggplot(data=mdr_pros, aes(x=Year, y=Freq, fill = MDR)) + geom_bar(stat = "identity")+
  labs(x="Year", y= "Number of isolates") + ylim(0, 100) +
  theme(text = element_text(size = 20))
diversity(table(dfAST$MDR_YES, dfAST$Strain), index = "simpson")

#AMR mechanisms
table(dfAST$macro_res, dfAST$ery_int)
table(dfAST$cot_res, dfAST$cot_int)
table(dfAST$lev_wgs, dfAST$lev_int)
table(df$PBP1a, df$Png_int)
table(df$PBP2b, df$Png_int)
table(df$PBP2x, df$Png_int)

#GPSC and AMR
gpsc_amr <- select(df, c(6, 17, 19, 21, 24, 26, 30, 32, 34, 36, "Strain"))
gpsc_amr <- gpsc_amr[gpsc_amr$Year %in% c("2020", "2021", "2022", "2023"), ]

# for multiple antibiotics
for (i in c(2:10)) {
  gpsc_amr2 <- table(gpsc_amr[ , "Strain"], gpsc_amr[ , i])
  print(colnames(gpsc_amr)[i])
  print(gpsc_amr2)
}

## MDR and GPSC
gpsc_amr2 <- data.frame(table(gpsc_amr$Strain, gpsc_amr$MDR_YES))
colnames(gpsc_amr2) <- c("GPSC", "MDR", "Frequency")
gpsc_amr2 <- gpsc_amr2 %>%
  group_by(GPSC) %>%
  filter(sum(Frequency) >= 5)
ggplot(data=gpsc_amr2, aes(x=GPSC, y=Frequency, fill=MDR))+
  geom_bar(stat = "identity")
chisq.test(gpsc_amr$MDR_YES, gpsc_amr$Strain)

########### serotypes ############
df$Serotype_wgs <- as.factor(df$Serotype_wgs)
#serotypes among age groups
sero_age <- data.frame(table(df$ageGrp, df$Serotype_wgs))
colnames(sero_age) <- c("Age_group", "Serotype", "Frequency")
sero_age$Age_group <- as.factor(sero_age$Age_group)

ggplot(data=sero_age, aes(x=Frequency, y=reorder(Serotype, Frequency), fill=Age_group))+
  geom_bar(stat = "identity",  position = "stack") + 
  labs(x = "Number of isolates", y = "Serotype") + xlim(0, 60) +
  theme(text = element_text(size = 20))

chisq.test(df$ageGrp, df$Serotype_wgs)

# serotypes and prevaccine
sero_prevacc <- data.frame(table(df$vaccine_era, df$Serotype_wgs))
colnames(sero_prevacc) <- c("vaccine_era", "Serotype", "Frequency")
sero_prevacc <- sero_prevacc %>%
  group_by(Serotype) %>%
  filter(sum(Frequency) >= 5)
ggplot(data=sero_prevacc, aes(x=Serotype, y=Frequency, fill=vaccine_era))+
  geom_bar(stat = "identity") + theme(axis.text=element_text(size=14),
                                      axis.title=element_text(size=18))
fisher.test(df$vaccine_era, df$PCV13)

# Invasive and PCV13 coverage
table(df$Invasive.Non.invasive, df$PCV13)

#MDR and prevacc
table(df$MDR_YES, df$vaccine_era)
fisher.test(df$MDR_YES, df$vaccine_era)

#AMR and serotypes ; from 2020 only
amr_sero <- select(df, c(6, 17, 19, 21, 24, 26, 30, 32, 34, 36, "Serotype_wgs"))

amr_sero <- amr_sero[amr_sero$Year %in% c("2020", "2021", "2022", "2023"), ]
table(amr_sero$Serotype_wgs, amr_sero$Png_int)
table(amr_sero$Serotype_wgs, amr_sero$ceft_int)
table(amr_sero$Serotype_wgs, amr_sero$ery_int)
table(amr_sero$Serotype_wgs, amr_sero$tet_int)

# pcv13 AND GPSC
gpsc_pcv <- data.frame(table(df$Strain, df$PCV13))
colnames(gpsc_pcv) <- c("GPSC", "PCV13_type", "Frequency")
gpsc_pcv <- gpsc_pcv %>%
  group_by(GPSC) %>%
  filter(sum(Frequency) >= 5)
ggplot(data=gpsc_pcv, aes(x=GPSC, y=Frequency, fill=PCV13_type))+
  geom_bar(stat = "identity")

