## script to run stats for repeats
# load libraries
library(tidyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(FSA)

#main function
perform_analysis <- function(resDF, commonVFs, highGPSC, repType) {
  dataDF <- resDF %>%
    filter(Gene %in% commonVFs, GPSC %in% highGPSC) %>%
    select(Gene, GPSC, Upstream)
  colnames(dataDF) <- c("Gene", "GPSC", "value")
  if (length(unique(dataDF$value)) > 1) {
     shapiro_res <- dataDF %>%
      group_by(Gene) %>%
      summarise(test_statistic = shapiro.test(value)$statistic, 
                p_value = shapiro.test(value)$p.value,
                gene_mean = mean(value),
                gene_median = median(value),
                gene_q1 = quantile(value, 0.25),
                gene_q3 = quantile(value, 0.75))
    shapiro_res$adjP <- p.adjust(shapiro_res$p_value, method = "BH")
    print(shapiro_res)
    kruskall_res <- dataDF %>%
      group_by(Gene) %>%
      summarise(test_res = kruskal.test(value ~ GPSC)$statistic,
                test_p = kruskal.test(value ~ GPSC)$p.value)
    kruskall_res$adjP <- p.adjust(kruskall_res$test_p, method = "BH")
    print(kruskall_res)
    for (gene in commonVFs) {
    gene_res <- dunnTest(value ~ as.factor(GPSC), 
                         data = dataDF, 
                         method = "bh")
    gene_comp <- gene_res$res$Comparison
    gene_padj <- gene_res$res$P.adj
    gene_df <- data.frame("Comparison" = gene_comp, "adjP" = gene_padj)
    gene_df["Gene"] <- rep(gene, times = nrow(gene_df))
    fname <- paste(c(gene, repType, "csv"), collapse = ".")
    write.csv(gene_df, fname, row.names = FALSE)
  }
    
    
  }
  
}

## generating images
gen_plot <- function(repDF, repType,  commonVFs, highGPSC) {
  dataDF <- repDF %>%
    filter(Gene %in% commonVFs, GPSC %in% highGPSC) %>%
    select(Gene, GPSC, Upstream)
  colnames(dataDF) <- c("Gene", "GPSC", "value")
  plt_sprite <- ggplot(dataDF, aes(x = GPSC, y= value)) + geom_boxplot() + 
    facet_wrap(~ Gene, scales = "free") + labs(y = "Number of repeats upstream") + 
    theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8))
  plt_name <- paste(c(repType, "boxPlot", "png"), collapse = ".")
  ggsave(filename = plt_name, plot = plt_sprite, dpi = 600, width = 9, height = 5)
  
  
}
