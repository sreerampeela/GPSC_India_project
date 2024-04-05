## heatmap for virulence genes
library(xlsx)
library(dplyr)
library(ggplot2)
library(reshape2)
df <- read.xlsx("data_analysis.xlsx", sheet = "vfs_new")
#df <- select(df, -"NUM_FOUND")
#rownames(df) <- df$GPSC
#df <- select(df, -"GPSC")
df_long <- melt(data = df)
heatmap_new <- ggplot(data=df_long, aes(x=variable, 
                                        y= VF)) +
  geom_tile(aes(fill = value), colour = "white")+
  scale_fill_gradient(low = "blue", high = "red", na.value ="white")+
  labs(title = "Virulence genes in core genomes")+
  xlab("GPSC type")+
  ylab("Virulence gene")


## script for enrichment table
df <- read.xlsx("enrichment_cleanedData.xlsx", header = T, sheetName = "commonGOs")
# scale data using MS EXCEL functions
df_data <- select(df, c("term.name","gpsc1_SCALED", "gpsc2_SCALED", "gpsc8_SCALED", 
                        "gpsc9_SCALED","gpsc10_SCALED", "ippg_SCALED"))
colnames(df_data) <- c("GO_Term", "GPSC1","GPSC2","GPSC8","GPSC9","GPSC10","IPPG")
df_long <- melt(df_data)

plt <- ggplot(data=df_long, aes(x=variable, y=GO_Term)) +
  geom_tile(aes(fill=value)) +
  scale_fill_gradient(low = "blue", high = "red", na.value ="white")+
  xlab("GPSC type")+
  ylab("GO Term")+
  guides(fill = guide_colourbar(title = "Z score"))

df <- read.xlsx("enrichment_cleanedData.xlsx", header = T, sheetName = "gpsc10_unique")
df_data <- select(df, c("vf_pres", "term.name", "description", "gpsc10_prop", "gpsc10_fdr"))

#df_highProp <- df_data[df_data$gpsc10_prop >= 0.9, ]
df_vf <- df_data[df_data$vf_pres == 1, ]
ggplot(data=df_vf, aes(x=-log10(gpsc10_fdr), y=description, size = gpsc10_prop)) + geom_point()

## plot for num virulence genes vs GO terms. colored on FDR values and bubble
## size based on proportion of genes in GPSC10df_data <- select(df, c("vf_pres", "term.name", "description", "gpsc10_prop", "gpsc10_fdr"))

df_data <- select(df, c("vf_pres", "term.name", "description", "gpsc10_prop", "gpsc10_fdr", "num_vfs"))
df_gpsc10_most <- df_data[!is.na(df_data$num_vfs), ] # removed all NAs

ggplot(data=df_gpsc10_most, aes(x=num_vfs, y=term.name, size = gpsc10_prop)) + 
  geom_point(data=df_gpsc10_most, aes(x=num_vfs, y=term.name, color=-log10(gpsc10_fdr))) + 
  xlab("Number of virulence genes") + ylab("GO Term") + scale_size(name="Proportion") +  
  scale_x_continuous(breaks=c(0:10), 
                     labels=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))


## for COG categories
df <- read.xlsx("eggNOG_ann.xlsx", header = T, sheetName = "COG_categories")
df_prop <- select(df, c("Category", "coreGenome_z",	"shellGenome_z",	"cloudGenome_z"))
colnames(df_prop) <- c("Category", "Core genome", "shell genome", "cloud genome")
dfmelt <- melt(df_prop)
dfmelt <- na.omit(dfmelt)
ggplot(data=dfmelt, aes(x=variable, y=Category, size = value, color = variable, 
                        transparency = 0.6))+geom_point()

ggplot(data=dfmelt, aes(x=variable, y=Category))+
  geom_tile(aes (fill = value), colour = "white")+
  scale_fill_gradient(low = "#c2bbae", high = "#1f06a0")+
  xlab("Genome compartment")
  
