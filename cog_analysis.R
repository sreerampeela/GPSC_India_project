library(xlsx)
library(reshape2)
library(ggplot2)
library(tidyr)
library(dplyr)

### CORE GENOME
df <- read.xlsx("EMAPPER_OUT.xlsx", sheetName = "TABLES")

df_prop <- select(df, c(1, 10:12))
df_prop <- transform(df_prop, sum_val = rowSums(df_prop[ , 2:4]))
df_clean <- df_prop[df_prop$sum_val != 0, 1:4]
df_prop_long <- melt(df_clean)
ggplot(data=df_prop_long, aes(x=variable, y=Category, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient(low="white", high = "red")


## SHELL GENOME
df_shell <- read.xlsx("EMAPPER_OUT.xlsx", sheetName = "SHELL_TABLE")
df_shell <- select(df_shell, c(1, 5:7))
df_shell <- na.omit(df_shell)
df_shell_long <- melt(df_shell)
ggplot(data=df_shell_long, aes(x=variable, y=Category, fill=value)) + 
  geom_tile(label = "percentage of genes") + 
  scale_fill_gradient(low="white", high = "red")+
  labs(x = "species")


## CLOUD GENOME
df_cloud <- read.xlsx("EMAPPER_OUT.xlsx", sheetName = "cloud_table")
df_cloud <- select(df_cloud, c(1, 5:7))
df_cloud <- na.omit(df_cloud)
df_cloud_long <- melt(df_cloud)
ggplot(data=df_cloud_long, aes(x=variable, y=Category, fill=value)) + 
  geom_tile(label = "percentage of genes") + 
  scale_fill_gradient(low="white", high = "red")+
  labs(x = "species")
