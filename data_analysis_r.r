# serotype paper data analysis
library(xlsx)
library(dplyr)
library(reshape2)
library(ggplot2)
library(tidyr)
### Reading a data file
df <- read.xlsx("merged_data_twostudies.xlsx", header= T,
                sheetName = "merged_data_2020_2023")
# convert all blank cells to NA
#df[df == " "] <- NA
#df

# year-wise plots
year_inv <- data.frame(table(df$Date, df$Invasive_Non.invasive))
colnames(year_inv) <- c("Year", "Disease", "Freq")
ggplot(data=year_inv, aes(x = Year, y = Freq, fill = Disease)) + 
  geom_bar(stat = "identity")+
  ylim(c(0, 75))

# pcv type and ast
ast_pcv_period <- select(df, c(4,5,7,9,11,13,17,19,21, 31))
# creating a new dataframe
ast_era <- data.frame(matrix(ncol = 2, nrow = 9))
row.names(ast_era) <- colnames(df)[c(4,5,7,9,11,13,17,19,21)]
colnames(ast_era) <- unique(df$PCV_PERIOD)

for (i in 1:9) {
  ast_var <- as.character(row.names(ast_era)[i])
  tab_data <- table(df$PCV_PERIOD, df[ , ast_var])
  #print(tab_data)
  if (ast_var == "MDR_INT") {
    pre_susceptible <- tab_data[2,1]*100/rowSums(tab_data)[2]
    post_susceptible <- tab_data[1,1]*100/rowSums(tab_data)[1]
  } else {
    pre_susceptible <- tab_data[2,ncol(tab_data)]*100/rowSums(tab_data)[2]
    post_susceptible <- tab_data[1,ncol(tab_data)]*100/rowSums(tab_data)[1]
  }
  ast_era[i,1] <- pre_susceptible
  ast_era[i,2] <- post_susceptible
}

ast_era[, "Antibiotic"] <- row.names(ast_era)
df_long <- melt(ast_era)
colnames(df_long) <- c("Antibiotic", "Vaccine_era", "Percentage_isolates")
ggplot(data=df_long, aes(x = Antibiotic, y = Percentage_isolates, fill =  Vaccine_era)) + 
  geom_bar(stat = "identity", position = "dodge")

# ast pcv13-nvt
ast_pcv13 <- select(df, c(4,5,7,9,11,13,17,19,21, 33))
ast_pcv13_df <- data.frame(matrix(ncol = 2, nrow = 10))
row.names(ast_pcv13_df) <- colnames(df)[c(4,5,7,9,11,13,17,19,21, 33)]
colnames(ast_pcv13_df) <- unique(df$PCV13_TYPES)
for (i in 1:9) {
  ast_var <- as.character(row.names(ast_pcv13_df)[i])
  tab_data <- table(ast_pcv13$PCV13_TYPES, df[ , ast_var])
  print(ast_var)
  print(tab_data)
  if (ast_var == "MDR_INT") {
    pcv_susceptible <- tab_data[2,1]*100/rowSums(tab_data)[2]
    nvt_susceptible <- tab_data[1,1]*100/rowSums(tab_data)[1]
  } else {
    pcv_susceptible <- tab_data[2,ncol(tab_data)]*100/rowSums(tab_data)[2]
    nvt_susceptible <- tab_data[1,ncol(tab_data)]*100/rowSums(tab_data)[1]
  }
  ast_pcv13_df[i,1] <- pcv_susceptible
  ast_pcv13_df[i,2] <- nvt_susceptible
  
}

ast_pcv13_df[, "Antibiotic"] <- row.names(ast_pcv13_df)


# Plotting
df_long <- melt(ast_pcv13_df)
colnames(df_long) <- c("Antibiotic", "PCV13_type", "Percentage_isolates")
ggplot(data=df_long, aes(x = Antibiotic, y = Percentage_isolates, fill =  PCV13_type)) + 
  geom_bar(stat = "identity", position = "dodge")

# AST vs pneumosil

ast_pneumosil <- select(df, c(4,5,7,9,11,13,17,19,21, 34))
ast_pneumosil_df <- data.frame(matrix(ncol = 2, nrow = 9))
row.names(ast_pneumosil_df) <- colnames(df)[c(4,5,7,9,11,13,17,19,21)]
colnames(ast_pneumosil_df) <- unique(ast_pneumosil$PNEUMOSIL)
for (i in 1:9) {
  ast_var <- as.character(row.names(ast_pneumosil_df)[i])
  tab_data <- table(ast_pneumosil$PNEUMOSIL, df[ , ast_var])
  print(ast_var)
  print(tab_data)
  if (ast_var == "MDR_INT") {
    pcv_susceptible <- tab_data[2,1]*100/rowSums(tab_data)[2]
    nvt_susceptible <- tab_data[1,1]*100/rowSums(tab_data)[1]
  } else {
    pcv_susceptible <- tab_data[2,ncol(tab_data)]*100/rowSums(tab_data)[2]
    nvt_susceptible <- tab_data[1,ncol(tab_data)]*100/rowSums(tab_data)[1]
  }
  ast_pneumosil_df[i,1] <- pcv_susceptible
  ast_pneumosil_df[i,2] <- nvt_susceptible
  
}

ast_pneumosil_df[, "Antibiotic"] <- row.names(ast_pneumosil_df)
df_long <- melt(ast_pneumosil_df)
colnames(df_long) <- c("Antibiotic", "Pneumosil_type", "Percentage_isolates")
ggplot(data=df_long, aes(x = Antibiotic, y = Percentage_isolates, fill =  Pneumosil_type)) + 
  geom_bar(stat = "identity", position = "dodge")

# PEN MIC
png_mic <- select(df, c(3,23))
boxplot(as.numeric(png_mic$PEN) ~ as.factor(png_mic$Date), 
        xlab = "YEAR", ylab = "Penicillin MIC (µg/ml)", ylim = c(0,10))


# Serotypes pre and post vacc
sero_prevacc <- select(df, c(22,23,25, 31))
sero_table <- sero_prevacc %>%
  group_by(PCV_PERIOD, Invasive_Non.invasive, Date) %>%
  count(SEROTYPE)
colnames(sero_table)[5] <- "Numb"

ggplot(data = sero_table, aes(x = SEROTYPE, y= Numb, 
                              fill = PCV_PERIOD)) + 
  geom_bar(stat = "identity", position="fill") + facet_wrap(~Invasive_Non.invasive) +
  ylab("Proportion of Isolates") + xlab("Serotype/group") + coord_flip()

sero_table[sero_table$Numb >= 5, ]
#CREATE TABLE for each variable and run fisher exact test
ceftnm <- table(df$ceft_int_nm, df$PCV_PERIOD)
ceftnm
fisher.test(ceftnm)

pngnm <- table(df$Png_int_nm, df$PCV_PERIOD)
pngnm
fisher.test(pngnm)

ery <- table(df$ery_int, df$PCV_PERIOD)
ery
fisher.test(ery)

ceftmen <- table(df$ceft_int_men, df$PCV_PERIOD)
ceftmen
fisher.test(ceftmen)

clin <- table(df$cln_int, df$PCV_PERIOD)
clin
fisher.test(clin)

tet <- table(df$tet_int, df$PCV_PERIOD)
tet
fisher.test(tet)

lev <- table(df$lev_int, df$PCV_PERIOD)
lev
fisher.test(lev)

cot <- table(df$cot_int, df$PCV_PERIOD)
cot
fisher.test(cot)

chlor <- table(df$chlor_int, df$PCV_PERIOD)
chlor
fisher.test(chlor)

## create a loop to run fisher test
#antibiotics <- c("Png_int_men", "Png_int_nm", "ery_int", "ceft_int_men", "ceft_int_nm",
                 "cln_int", "tet_int", "lev_int", "cot_int", "chlor_int")
for (i in 2:11) {
  varval <- table(df[,i], df$PCV13_TYPES)
  print(colnames(df)[i])
  print(varval)
  print(fisher.test(varval))
}

##PCV13 SEROGROUPS
for (i in 2:11) {
  varval <- table(df[,i], df$PCV13_GRPS)
  print(colnames(df)[i])
  print(varval)
  print(fisher.test(varval))
}

for (i in 2:11) {
  varval <- table(df[,i], df$PNEUMOSIL_TYPES)
  print(colnames(df)[i])
  print(varval)
  print(fisher.test(varval))
}

for (i in 2:11) {
  varval <- table(df[,i], df$PNEUMOSIL_GRPS)
  print(colnames(df)[i])
  print(varval)
  print(fisher.test(varval))
}



#library(ggplot2)

#png_pcv13 <- ggplot(df, aes(x=df$PCV13_TYPES, y = df$Png_int_men[df$Png_int_men == "Yes"]))
#png_pcv13 <- png_pcv13 + geom_bar(stat = "identity", position = "dodge")
#png_pcv13

data_long <- reshape(df, idvar = "PCV13_TYPES", timevar = "Png_int_men", direction = "wide")
data_long
row.names(data_long) <- data_long$PCV13_TYPES
data_long <- data_long[, 2:ncol(data_long)]
data_matrix <- as.matrix(data_long)
barplot(height = tet,                       
        beside = TRUE, legend.text  = TRUE)


### box plot of png
df2 <- read.csv("png_data.csv", header = T)
print(df2$X2015)
png_box <- boxplot(df2, ylab = "PEN_MIC")
PNG_ANALYSIS <- boxplot.stats(df2$X2015)

df3 <- read.csv("png_data_LONG.csv", header = T)
print(df3$PEN_MIC)
library(ggplot2)
years <- as.factor(df3$YEAR)
ggplot(df3, aes(x = years , y = PEN_MIC)) +
  geom_boxplot() + labs(x = "YEAR", y = "PEN MIC") + theme_classic()

## pcv13 coverage among children <5 across pre and post-vaccine era
pcv13_ages <- df %>%
  group_by(Invasive_Non.invasive, age_group, PCV_PERIOD) %>%
  count(PCV13_TYPES)

colnames(pcv13_ages)[5] <- "Number_isolates"

pneumosil_ages <- df %>%
  group_by(Invasive_Non.invasive, age_group, PCV_PERIOD) %>%
  count(PNEUMOSIL)

# MDR and serotypes
predominant_serotypes <- c("14", "15B", "19A", "19F", "23F", "3", "6A", "6B")

common_sero_mdr <- df %>%
  filter(SEROTYPE %in% predominant_serotypes) %>%
  select(c(4,5,7,9,11,13,17,19,21, 25))

df_mdr_sero <- table(common_sero_mdr$SEROTYPE, common_sero_mdr$MDR_INT)
total_mdr = 139
total_non_mdr = 161
mdr_sero_p <- c()
odds_ratios <- c()
for (i in 1:nrow(df_mdr_sero)) {
  print(row.names(df_mdr_sero)[i])
  tabular_data <- matrix(c(df_mdr_sero[i, 2], df_mdr_sero[i, 1], 
                    total_mdr - df_mdr_sero[i, 2],  
                    total_non_mdr - df_mdr_sero[i, 1]), nrow = 2, ncol = 2, byrow = T)
  print(tabular_data)
  fisher_test <- fisher.test(tabular_data)
  mdr_sero_p <- c(mdr_sero_p, fisher_test$p.value)
  odds_ratios <- c(odds_ratios, fisher_test$estimate)
  
  
}

# AST PCV13 stats
# using ast_pcv13 dataset
ast_pcv13_df <- data.frame(matrix(ncol = 2, nrow = 9))
row.names(ast_pcv13_df) <- colnames(df)[c(4,5,7,9,11,13,17,19,21)]
colnames(ast_pcv13_df) <- unique(df$PCV13_TYPES)
pvals_pcv13_ast <- c()
odds_pcv13_ast <- c()
for (i in 1:9) {
  ast_var <- as.character(row.names(ast_pcv13_df)[i])
  tab_data <- table(ast_pcv13$PCV13_TYPES, df[ , ast_var])
  print(ast_var)
  print(tab_data)
  if (ast_var == "MDR_INT") {
    pcv_susceptible <- tab_data[2,1]
    nvt_susceptible <- tab_data[1,1]
    pcv_ns <- tab_data[2,2]
    npcv_ns <- tab_data[1,2]
  } else {
    pcv_susceptible <- tab_data[2,ncol(tab_data)]
    nvt_susceptible <- tab_data[1,ncol(tab_data)]
    pcv_ns <- rowSums(tab_data)[2] - pcv_susceptible
    npcv_ns <- rowSums(tab_data)[1] - nvt_susceptible
  }
  cont_table <- matrix(c(pcv_susceptible, pcv_ns, nvt_susceptible, npcv_ns),
                       nrow = 2, ncol = 2, byrow = T)
  stat_test <- fisher.test(cont_table)
  pvals_pcv13_ast <- c(pvals_pcv13_ast, stat_test$p.value)
  odds_pcv13_ast <- c(odds_pcv13_ast, stat_test$estimate)
  
}

adjps <- p.adjust(pvals_pcv13_ast, method = "bonferroni")

# AST Pneumosil stats
ast_pneumo <- select(df, c(4,5,7,9,11,13,17,19,21, 35))
pvals_pneumo_ast <- c()
odds_pneumo_ast <- c()
for (i in 1:9) {
  ast_var <- as.character(colnames(ast_pneumo)[i])
  tab_data <- table(ast_pneumo$PNEUMOSIL, df[ , ast_var])
  print(ast_var)
  print(tab_data)
  if (ast_var == "MDR_INT") {
    pcv_susceptible <- tab_data[2,1]
    nvt_susceptible <- tab_data[1,1]
    pcv_ns <- tab_data[2,2]
    npcv_ns <- tab_data[1,2]
  } else {
    pcv_susceptible <- tab_data[2,ncol(tab_data)]
    nvt_susceptible <- tab_data[1,ncol(tab_data)]
    pcv_ns <- rowSums(tab_data)[2] - pcv_susceptible
    npcv_ns <- rowSums(tab_data)[1] - nvt_susceptible
  }
  cont_table <- matrix(c(pcv_susceptible, pcv_ns, nvt_susceptible, npcv_ns),
                       nrow = 2, ncol = 2, byrow = T)
  stat_test <- fisher.test(cont_table)
  pvals_pneumo_ast <- c(pvals_pneumo_ast, stat_test$p.value)
  odds_pneumo_ast <- c(odds_pneumo_ast, stat_test$estimate)
  
}

adjps <- p.adjust(pvals_pneumo_ast, method = "bonferroni")

## common serotypes vs AMR - PRE and POST
predominant_serotypes <- c("14", "15B", "19A", "19F", 
                           "23F", "3", "6A", "6B")

df_amr_sero <- select(df, c(4,6,9,11, 13, 15, 19, 21, 23, 27, 34))
df_amr_common_sero <- df_amr_sero %>%
  filter(SEROTYPE %in% predominant_serotypes)
prevacc_df <- df_amr_common_sero[df_amr_common_sero$PCV_PERIOD == "PRE", ]
postvacc_df <- df_amr_common_sero[df_amr_common_sero$PCV_PERIOD == "POST", ]
pvals <- c()
ors <- c()
seros <- c()
ast_s_pre <- c()
ast_ns_pre <- c()
ast_s_post <- c()
ast_ns_post <- c()
amrs <- c()
for (i in 1:9) {
  amr_id <- colnames(df_amr_common_sero)[i]
  print(amr_id)
  for (sero in predominant_serotypes) {
    if (amr_id == "MDR_INT") {
      pre_s <- nrow(prevacc_df[prevacc_df[[amr_id]] == "N" & prevacc_df$SEROTYPE == sero, ])
      pre_ns <- nrow(prevacc_df[prevacc_df[[amr_id]] == "Y" & prevacc_df$SEROTYPE == sero, ])
      post_s <- nrow(postvacc_df[postvacc_df[[amr_id]] == "N" & postvacc_df$SEROTYPE == sero , ])
      post_ns <- nrow(postvacc_df[postvacc_df[[amr_id]] == "Y" & postvacc_df$SEROTYPE == sero , ])
    } else {
      pre_s <- nrow(prevacc_df[prevacc_df[[amr_id]] == "S" & prevacc_df$SEROTYPE == sero, ])
      pre_ns <- nrow(prevacc_df[prevacc_df[[amr_id]] == "NS" & prevacc_df$SEROTYPE == sero, ])
      post_s <- nrow(postvacc_df[postvacc_df[[amr_id]] == "S" & postvacc_df$SEROTYPE == sero, ])
      post_ns <- nrow(postvacc_df[postvacc_df[[amr_id]] == "NS" & postvacc_df$SEROTYPE == sero, ])
  
    }
    
    cont_table <- matrix(c(post_ns, post_s, pre_ns, pre_s), nrow = 2,
                         ncol = 2, byrow = T)
    row.names(cont_table) <- c("POST", "PRE")
    colnames(cont_table) <- c("NS", "S")
    print(cont_table)
    cont_table[cont_table == 0] <- 0.5
    fisher_stats <- fisher.test(cont_table)
    pvals <- c(pvals, fisher_stats$p.value)
    ors <- c(ors, fisher_stats$estimate)
    seros <- c(seros, sero)
    ast_s_pre <- c(ast_s_pre, pre_s)
    ast_ns_pre <- c(ast_ns_pre, pre_ns)
    ast_s_post <- c(ast_s_post, post_s)
    ast_ns_post <- c(ast_ns_post, post_ns)
    amrs <- c(amrs, amr_id)
    
  }
  
}

# final df
ast_common_sero_pre <- data.frame("SEROTYPE" = seros,
                                  "Antibiotic" = amrs,
                                  "PRE_SUS" = ast_s_pre,
                                  "PRE_NS" = ast_ns_pre,
                                  "POST_SUS" = ast_s_post,
                                  "POST_NS" = ast_ns_post,
                                  "Odds" = ors,
                                  "Raw P" =  pvals)
write.csv(ast_common_sero_pre, "common_serotypes_ast_pre_post.csv")
# Manually recalculate OR if any cell values = 0
ast_res <- read.csv("common_serotypes_ast_pre_post.csv", header = T)
adjpvals <- p.adjust(ast_res$Raw.P, method = "bonferroni")
ast_res["AdjP"] <- adjpvals
write.csv(ast_res, "top_sero_ast_pre.csv")
