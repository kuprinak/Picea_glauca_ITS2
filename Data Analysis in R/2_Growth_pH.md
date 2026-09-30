---
title: "Growth_pH"
author: "Kristina Kuprina"
date: "`r Sys.Date()`"
output:
  html_document: default
---

```{r Packages, message=FALSE, warning=FALSE}
library("ggplot2")
library("GGally")
library("Hmisc")
library("readr")
library("dplyr")
library("networkD3")
library("tidyr")
library('tibble')
library('pals')
library("rstatix")
library("cowplot")
library("ggrepel")
```

```{r session info}
sessionInfo()
```

## Correlation matrix (Tree growth)

```{r fig.height=8, fig.width=8, warning=FALSE}
sample_info_no_R <-read.table("Alaska_info_no_R.txt", row.names=1, header=T, check.names=F)

df_correlation <- sample_info_no_R[, c("BAI_10y", "detBAI_10y", "rcsBAI_10y", "Age", "DBH", "Height", "pH", "Site")]
df_correlation <- na.omit(df_correlation)

# Compute correlations + p-values
cor_results <- rcorr(as.matrix(df_correlation[, 1:7])) # column 9 (Site) is used for point colours
cor_matrix <- cor_results$r
p_matrix <- cor_results$P

site_colors <- c(Alaska_Range = "#EE7733", Brooks_Range = "#009988", Interior = "#5E4FA2")

lower_fun_sig <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point(aes(color = Site), alpha = 0.7, size = 1.5) +
    scale_color_manual(values = site_colors,
                       labels = c("Alaska Range", "Brooks Range", "Interior Alaska")) +
    geom_smooth(method = "lm", se = FALSE, color = "darkred", size = 0.7)}

ggpairs(df_correlation, columns = 1:7,
        lower = list(continuous = lower_fun_sig),
        upper = list(continuous = wrap("cor", size = 4)),
        diag = list(continuous = wrap("densityDiag", alpha = 0.5, fill = "lightblue")),
        legend = c(2, 1),
        progress = FALSE) +
theme_bw() + theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "lightgray"),
        axis.text = element_text(size = 8),
        legend.position = "bottom")
```

## Differences between Plots 

```{r fig.height=10, fig.width=10}
palette <- c("#EE7733","#CC6677", "#009988","#337538","#33BBEE", "#5E4FA2")

pH <- ggplot(sample_info_no_R, aes(x = Plot, y = pH, fill = Plot)) + geom_violin(trim = FALSE) +geom_point(size = 2, color = "gray10", alpha = 0.3,
             position = position_jitter(width = 0.1)) + scale_fill_manual(values = palette) + ylab("Soil pH") + theme_minimal() + theme(legend.position = "none")
Age_plot <- ggplot(sample_info_no_R, aes(x = Plot, y = Age, fill = Plot)) + geom_violin(trim = FALSE) +geom_point(size = 2, color = "gray10", alpha = 0.3,
             position = position_jitter(width = 0.1)) + scale_fill_manual(values = palette) + ylab("Tree age, years") + theme_minimal() + theme(legend.position = "none")
DBH_plot <- ggplot(sample_info_no_R, aes(x = Plot, y = DBH, fill = Plot)) + geom_violin(trim = FALSE) +geom_point(size = 2, color = "gray10", alpha = 0.3,
             position = position_jitter(width = 0.1)) + scale_fill_manual(values = palette) + ylab("DBH, cm") + theme_minimal() + theme(legend.position = "none")
Height_plot <- ggplot(sample_info_no_R, aes(x = Plot, y = Height, fill = Plot)) + geom_violin(trim = FALSE) +geom_point(size = 2, color = "gray10", alpha = 0.3, 
              position = position_jitter(width = 0.1)) + scale_fill_manual(values = palette) + ylab("Tree height, m") + theme_minimal() + theme(legend.position = "none")

plot_grid(pH, Age_plot, DBH_plot, Height_plot, labels = c("Soil pH", "Tree age", "Tree DBH", "Tree height" ), ncol = 2)
```

## Soil pH in BF (map)

```{r BF_pH_map, fig.height=4.5, fig.width=6.5}
BF_pH <- subset(sample_info_no_R, Plot == "BF" & !is.na(pH))
BF_pH$Sample <- rownames(BF_pH)

ggplot(BF_pH, aes(x = Longitude, y = Latitude)) +
  geom_point(aes(fill = pH), shape = 21, size = 4, color = "gray20") +
  geom_text_repel(aes(label = paste0(Sample, "\npH ", pH)), size = 2.8, point.size = 4, max.overlaps = Inf, lineheight = 0.9) +
  scale_fill_viridis_c(name = "Soil pH") +
  coord_quickmap() +
  theme_bw()
```

## Statistical tests

```{r pH vs habitat/Plot - test}
shapiro.test(sample_info_no_R$pH) # data are not normally distributed
bartlett.test(pH ~ Plot, data = sample_info_no_R) # variance is not equal

dunn_test(sample_info_no_R , pH ~ Habitat) #not significant 
dunn_test(sample_info_no_R , pH ~ Plot, p.adjust.method = "bonferroni") 
dunn_test(sample_info_no_R , pH ~ Site, p.adjust.method = "bonferroni") 
```

```{r Tree age vs habitat - test}
shapiro.test(sample_info_no_R$Age) # data are not normally distributed
bartlett.test(Age ~ Habitat, data = sample_info_no_R) # variance is equal
kruskal_test(sample_info_no_R , Age ~ Habitat) #significant
```

```{r Tree age vs Plot - test}
shapiro.test(sample_info_no_R$Age) # data are not normally distributed
bartlett.test(Age ~ Plot, data = sample_info_no_R) # variance is not equal

res.kruskal <-kruskal_test(sample_info_no_R , Age ~ Plot)
res.kruskal #significant

dunn_test(sample_info_no_R , Age ~ Plot, p.adjust.method = "bonferroni") 
```

## DBH vs Habitat/Plot 

```{r DBH vs Habitat - test}
shapiro.test(sample_info_no_R$DBH) # data are not normally distributed
bartlett.test(DBH ~ Habitat, data = sample_info_no_R) # variance is not equal
dunn_test(sample_info_no_R , DBH ~ Habitat) # not significant
```

```{r DBH vs Plot - test}
shapiro.test(sample_info_no_R$DBH) # data are not normally distributed
bartlett.test(DBH ~ Plot, data = sample_info_no_R) # variance is not equal

res.kruskal <-kruskal_test(sample_info_no_R , DBH ~ Plot)
res.kruskal #significant

dunn_test(sample_info_no_R , DBH ~ Plot, p.adjust.method = "bonferroni") 
```

## Height vs Habitat/Plot 

```{r Height vs Habitat - test}
shapiro.test(sample_info_no_R$Height) # data are not normally distributed
bartlett.test(Height ~ Habitat, data = sample_info_no_R) # variance is equal
dunn_test(sample_info_no_R , Height ~ Habitat) # significant
```

```{r Height vs Plot -test}
shapiro.test(sample_info_no_R$Height) # data are not normally distributed
bartlett.test(Height ~ Plot, data = sample_info_no_R) # variance is not equal

res.kruskal <-kruskal_test(sample_info_no_R , Height ~ Plot)
res.kruskal #significant
dunn_test(sample_info_no_R , Height ~ Plot, p.adjust.method = "bonferroni") 
```



















