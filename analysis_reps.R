library(dplyr)
library(reshape2)
# library(ggpubr)
library(ggplot2)
# theme_set(theme_pubr())
library(dunn.test)

df_sm <- read.csv("smitis_all.csv", header = T)

df_sps <- read.csv("sps_all.csv", header = T)

df_spn_jip <- read.csv("jip_2024_reps.csv", header = T)

######## PLOTS ##############
df_jip <- melt(df_spn_jip)
colnames(df_jip) <- c("sample", "repType", "number")
plt1 <- ggplot(data=df_jip, aes(y=number, color = repType)) + 
  geom_boxplot() 

df_sm_long <- melt(df_sm)
colnames(df_sm_long) <- c("sample", "repType", "number")
plt2 <- ggplot(data=df_sm_long, aes(y=number, color = repType)) + 
  geom_boxplot()

df_sps_long <- melt(df_sps)
colnames(df_sps_long) <- c("sample", "repType", "number")
plt3 <- ggplot(data=df_sps_long, aes(y=number, color = repType)) + 
  geom_boxplot()


plts <- list(plt1, plt2, plt3)
ggarrange(plotlist = plts,
          labels = c("pneumococcus (JIPMER)", "S.mitis", "S.pseudopneumoniae"),
          ncol = 3, nrow = 1, font.label = list(size = 18), align = "hv",
          common.legend = TRUE, legend = "bottom")

####### boxA repeats ############
jip_mean_boxA <- mean(df_spn_jip$boxA, na.rm = T)

sm_mean_boxA <- mean(df_sm$boxA, na.rm = T)

sps_mean_boxA <- mean(df_sps$boxA, na.rm = T)

## creating a combined dataframe
df_reps <- data.frame(matrix(data=NA, ncol = 6, nrow = sum(nrow(df_sm), 
                                                           nrow(df_sps), 
                                                           nrow(df_spn_jip))))

colnames(df_reps) <- c("species", "boxA", "boxB", "boxC", "RUP", "SPRITE")

df_reps$species <- c(rep(c("SPN"), times=197), rep(c("SM"), times = nrow(df_sm)), 
                     rep(c("SPS"), times = nrow(df_sps)))
df_reps$boxA <- c(df_spn_jip$boxA, df_sm$boxA, df_sps$boxA)
df_reps$boxB <- c(df_spn_jip$boxB, df_sm$boxB, df_sps$boxB)
df_reps$boxC <- c(df_spn_jip$boxC, df_sm$boxC, df_sps$boxC)
df_reps$RUP <- c(df_spn_jip$RUP, df_sm$RUP, df_sps$RUP)
df_reps$SPRITE <- c(df_spn_jip$SPRITE, df_sm$SPRITE, df_sps$SPRITE)

df_reps$species <- as.factor(df_reps$species)

## using for loop for shapiro test
for (i in 2:6) {
  print(colnames(df_reps)[i])
  pval <- df_reps %>%
    group_by(species) %>%
    summarise(pval = shapiro.test(.[[i]])$p.value)
  print(pval)
  normvars <- filter(pval, pval > 0.05)
  # print(summary(normvars))
 if (nrow(normvars) < nlevels(df_reps$species)) {
   statres <- kruskal.test(df_reps[, i] ~ df_reps[, 1])
   posthocRes <- dunn.test::dunn.test(df_reps[, i], g = df_reps[, 1],
                                   method = "bh")
 } else {
   #   # all groups have normal dist
     statres <- aov(df_reps[, i] ~ df_reps[, 1])
     posthocRes <- TukeyHSD(statres)
  print(statres)
  print(summary(posthocRes))

 }
}

### Plots colored by species type
df_reps_long <- melt(df_reps)

ggplot(df_reps_long, aes(x=variable, y=value, fill=species)) + geom_boxplot() + 
  labs(x = "Repeat type", y = "Number of repeats") + 
  theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),
        axis.title = element_text(size = 14), legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))
