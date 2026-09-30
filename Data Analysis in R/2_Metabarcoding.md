---
title: "Metabarcoding"
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
library("DESeq2")
library("phyloseq")
library('ggrepel')
library('ggalluvial')
library("metacoder")
library("microeco")
library("dplyr")
library("iNEXT")
library("networkD3")
library("tidyr")
library('tibble')
library('pals')
library("rstatix")
library("cowplot")
library("vegan")
```

```{r session info}
sessionInfo()
```

```{r sample info}
sample_info <-read.table("Alaska_info.txt", row.names=1, header=T, check.names=F)
sample_info_no_R <- read.table("Alaska_info_no_R.txt", row.names=1, header=T, check.names=F)

tax_tab <- read_delim("Alaska_taxonomy.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
count_tab <-read.table("Alaska_counts.txt", header=T, row.names=1, check.names=F)
count_tab <- count_tab[rownames(count_tab) %in% tax_tab$OTU, ] 
```

### Treatment of NTC

```{r Treatment of NTC - 1}
NTC <-rowSums(count_tab[, 1:8]) # sum of NTC counts for each OTU
sample_counts <- rowSums(count_tab[, 9:143]) # sum of sample counts for each OTU

# normalize sample counts by dividing the samples' total by 17 – as there are 17 as many samples (135) as there are NTCs (8) (135/8=17)
norm_sample_OTU_counts <- sample_counts/17

# which OTUs are deemed likely contaminants based on the threshold noted above:
blank_OTUs <- names(NTC[NTC *10 > norm_sample_OTU_counts]) # OTUs which have more counts than 1/10 of the sample counts
length(blank_OTUs) 

# How many counts are kept
colSums(count_tab[!rownames(count_tab) %in% blank_OTUs, ]) / colSums(count_tab) * 100

filt_count_tab <- count_tab[!rownames(count_tab) %in% blank_OTUs, -c(1:8)] # removing blank_OTUs and the blank samples from further analysis
filt_count_tab <- filt_count_tab[rowSums(filt_count_tab == 0) != ncol(filt_count_tab), ] # delete OTUs with 0 counts
```

### Normalizing for sampling depth (variance stabilizing transformation) 

```{r}
sample_info<-sample_info[order(rownames(sample_info)), ]
filt_count_tab1<- filt_count_tab[,order(names(filt_count_tab))]

deseq_counts <- DESeqDataSetFromMatrix(filt_count_tab1, colData = sample_info, design = ~ Site) 
deseq_counts <- estimateSizeFactors(deseq_counts,type = "poscounts")
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts,blind = TRUE)
vst_trans_count_tab <- assay(deseq_counts_vst) #save as matrix
vst_trans_count_tab_df <- as.data.frame(vst_trans_count_tab)
```

### PCoA (Plot)

```{r fig.height=5, fig.width=8}
palette <- c("#EE7733","#CC6677", "#009988","#337538","#33BBEE", "#5E4FA2")

# create phyloseq object with normalized data
vst_count_phy <- otu_table(vst_trans_count_tab, taxa_are_rows=T)
sample_info_tab_phy <- sample_data(sample_info)
vst_physeq <- phyloseq(vst_count_phy, sample_info_tab_phy)

#generating and visualizing the PCoA with phyloseq
vst_pcoa <- ordinate(vst_physeq, method="PCoA", distance="euclidean")

eigen_vals <- vst_pcoa$values$Eigenvalues # to scale the axes
total_var <- sum(eigen_vals)
axis1 <- 100 * eigen_vals[1] / total_var
axis2 <- 100 * eigen_vals[2] / total_var

plot_ordination(vst_physeq, vst_pcoa, color = "Plot", shape = "Habitat") +
  geom_point(size = 2) +
  geom_text_repel(aes(label = rownames(sample_info)), size = 2, force = 1, max.overlaps = 100) +
  stat_ellipse(aes(fill = Plot), geom = "polygon", color = NA, alpha = 0.1) +
  coord_fixed(sqrt(eigen_vals[2] / eigen_vals[1])) +
    labs(x = paste0("Principal Coordinate Axis 1 (", round(axis1, 2), "%)"),
    y = paste0("Principal Coordinate Axis 2 (", round(axis2, 2), "%)")) +
  theme_classic() +
  scale_shape_manual(values = c(20, 17)) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette)
```

## Without technical replicates

### Normalizing for sampling depth (variance stabilizing transformation) 

```{r Preparing dataset no R}
#reorder the row and column names
sample_info_no_R<-sample_info_no_R[order(rownames(sample_info_no_R)), ]
filt_count_tab<- filt_count_tab[,order(names(filt_count_tab))]
filt_count_tab_no_r <- filt_count_tab[colnames(filt_count_tab) %in% rownames(sample_info_no_R) ]
```

```{r create DESeq2 object: Growth no R}
# create a DESeq2 object
deseq_counts <- DESeqDataSetFromMatrix(filt_count_tab_no_r, colData = sample_info_no_R, design = ~ Habitat) 
deseq_counts <- estimateSizeFactors(deseq_counts,type = "poscounts")
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts,blind = TRUE)
vst_trans_count_tab <- assay(deseq_counts_vst) #save as matrix
```

### PCoA (Plot, All OTUs)

```{r fig.height=4, fig.width=8}
# create phyloseq object with normalized data
vst_count_phy <- otu_table(vst_trans_count_tab, taxa_are_rows=T)
sample_info_tab_phy <- sample_data(sample_info_no_R)
vst_physeq <- phyloseq(vst_count_phy, sample_info_tab_phy)

#generating and visualizing the PCoA with phyloseq
vst_pcoa <- ordinate(vst_physeq, method="PCoA", distance="euclidean")

eigen_vals <- vst_pcoa$values$Eigenvalues # ato scale the axes
total_var <- sum(eigen_vals)
axis1 <- 100 * eigen_vals[1] / total_var
axis2 <- 100 * eigen_vals[2] / total_var

plot_ordination(vst_physeq, vst_pcoa, color = "Plot", shape = "Habitat") +
  geom_point(size = 2, alpha = 0.5) +
  stat_ellipse(aes(fill = Plot), geom = "polygon", color = NA, alpha = 0.2) +
  coord_fixed(sqrt(eigen_vals[2] / eigen_vals[1])) +
    labs(x = paste0("Principal Coordinate Axis 1 (", round(axis1, 2), "%)"),
    y = paste0("Principal Coordinate Axis 2 (", round(axis2, 2), "%)")) +
  theme_classic() +
  scale_shape_manual(values = c(20, 17)) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette)
```

### PCoA (Plot, only ECM)

```{r fig.height=4, fig.width=8}
ECM_OTUs <- tax_tab$OTU[tax_tab$Guild %in% c("Ectomycorrhizal")]
length(ECM_OTUs)

Sapro_OTU<-tax_tab$OTU[tax_tab$Guild %in% c("Wood Saprotroph","Plant Saprotroph","Undefined Saprotroph", "Plant Pathogen")]
length(Sapro_OTU)

ECM_OTUs_present <- intersect(ECM_OTUs, rownames(vst_trans_count_tab))

# create phyloseq object with normalized data
vst_count_phy_ECM <- otu_table(vst_trans_count_tab[ECM_OTUs_present, ], taxa_are_rows=T)
sample_info_tab_phy_ECM <- sample_data(sample_info_no_R)
vst_physeq_ECM <- phyloseq(vst_count_phy_ECM, sample_info_tab_phy)

#generating and visualizing the PCoA with phyloseq
vst_pcoa_ECM <- ordinate(vst_physeq_ECM, method="PCoA", distance="euclidean")

eigen_vals_ECM <- vst_pcoa_ECM$values$Eigenvalues # to scale the axes
total_var <- sum(eigen_vals_ECM)
axis1 <- 100 * eigen_vals_ECM[1] / total_var
axis2 <- 100 * eigen_vals_ECM[2] / total_var

plot_ordination(vst_physeq_ECM, vst_pcoa_ECM, color = "Plot", shape = "Habitat") +
  geom_point(size = 2, alpha = 0.5) +
  stat_ellipse(aes(fill = Plot), geom = "polygon", color = NA, alpha = 0.2) +
  coord_fixed(sqrt(eigen_vals_ECM[2] / eigen_vals_ECM[1])) +
    labs(x = paste0("Principal Coordinate Axis 1 (", round(axis1, 2), "%)"),
    y = paste0("Principal Coordinate Axis 2 (", round(axis2, 2), "%)")) +
  theme_classic() +
  scale_shape_manual(values = c(20, 17)) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette)
```


## Relative Abundancy

```{r Relative Abundancy - preparing data, fig.height=16, fig.width=8}

tax_table<-as.data.frame(tax_tab)
row.names(tax_table) <- tax_table$OTU

tax_table <- tax_table[, -1]
tax_table<- tax_table[order(rownames(tax_table)),]

# create microeco file
meco_fungi <- microtable$new(sample_table = sample_info_no_R, otu_table = filt_count_tab, tax_table = tax_table)
 
mycol <- c("#A6CEE3" ,"#1F78B4" ,"#B2DF8A" ,"#33A02C" ,"#FB9A99", "#E31A1C" ,"#FDBF6F" ,"#FF7F00" ,"#CAB2D6" ,
           "#6A3D9A" ,"#FFFF99" ,"#B15928" , "#6fdc8c", "green4", "#004144")

t1 <- trans_abund$new(dataset = meco_fungi, taxrank = "Phylum", ntaxa = 6)
t2 <- trans_abund$new(dataset = meco_fungi, taxrank = "Class", ntaxa = 11)
t3 <- trans_abund$new(dataset = meco_fungi, taxrank = "Order", ntaxa = 11)
t4 <- trans_abund$new(dataset = meco_fungi, taxrank = "Family", ntaxa = 11)
t5 <- trans_abund$new(dataset = meco_fungi, taxrank = "Genus", ntaxa = 12)

p1<-t1$plot_bar(xtext_keep = TRUE, xtext_angle = 90,facet = "Plot",xtext_size = 6, barwidth = 1, color_values = mycol)
p2<-t2$plot_bar(xtext_keep = TRUE, xtext_angle = 90,facet = "Plot",xtext_size = 6, barwidth = 1,color_values = mycol)
p3<-t3$plot_bar(xtext_keep = TRUE, xtext_angle = 90,facet = "Plot",xtext_size = 6, barwidth = 1,color_values = mycol)
p4<-t4$plot_bar(xtext_keep = TRUE, xtext_angle = 90,facet = "Plot",xtext_size = 6, barwidth = 1,color_values = mycol)

library(gridExtra)
grid.arrange(p1, p2, p3, p4, ncol = 1)
```

### Relative abundane - Genus

```{r}
t5$plot_bar(xtext_keep = TRUE, xtext_angle = 90,facet = "Plot",xtext_size = 6, barwidth = 1,color_values = mycol)
```

# Prediction of the RAF composition (all OTUs): Variance partitioning

```{r}
meta <- as(sample_data(vst_physeq), "data.frame")
keep <- complete.cases(meta[, c("Age","pH","Habitat","Site")])
vst_mat <- as.data.frame(t(otu_table(vst_physeq)))
vst_mat_clean <- vst_mat[keep, ]
meta_clean <- meta[keep, ]

# List of predictors
env <- meta_clean[, c("pH", "Habitat")]
env$Habitat <- as.factor(env$Habitat)
space <- data.frame(Site = as.factor(meta_clean$Site))
growth <- meta_clean$Age

# Variance partitioning

vp <- varpart(vst_mat_clean, ~ pH, ~ Habitat, ~ Site, ~ Age, data = meta_clean)
vals <- vp$part$indfract$Adj.R.square[1:4]

plot(vp, Xnames = c("pH", "Habitat","Site",  "Age"), bg = c("#EE7733","#5E4FA2","#337538" ,"#33BBEE"), cex = 1, digits = 2)
```

## Redundance analyses - All OTUs

```{r}
rda_env <- rda(vst_mat_clean ~ pH + Habitat + Condition(Site), data = meta_clean)
anova.cca(rda_env, by = "term", permutations = 9999)

rda_space <- rda(vst_mat_clean ~ Site, data = meta_clean)
anova.cca(rda_space, permutations = 9999)

rda_age <- rda(vst_mat_clean ~ Age + Condition(pH + Habitat + Site), data = meta_clean)
anova.cca(rda_age, by = "term", permutations = 9999)
```

# Prediction of the RAF composition (only ECM OTUs): Variance partitioning

```{r}
meta <- as(sample_data(vst_physeq), "data.frame")
keep <- complete.cases(meta[, c("Age","pH","Habitat","Site")])
vst_mat <- as.data.frame(t(otu_table(vst_physeq)))
vst_mat_clean <- vst_mat[keep, ]
meta_clean <- meta[keep, ]

tax_data <- tax_tab
tax_data_ECM <- tax_data[tax_data$Guild == "Ectomycorrhizal", ]
vst_mat_clean_ECM <- vst_mat_clean[, colnames(vst_mat_clean) %in% tax_data_ECM$OTU]


# List of predictors
env <- meta_clean[, c("pH", "Habitat")]
env$Habitat <- as.factor(env$Habitat)
space <- data.frame(Site = as.factor(meta_clean$Site))
growth <- meta_clean$Age

# Variance partitioning

vp <- varpart(vst_mat_clean_ECM, ~ pH, ~ Habitat, ~ Site, ~ Age, data = meta_clean)
vals <- vp$part$indfract$Adj.R.square[1:4]

plot(vp, Xnames = c("pH", "Habitat","Site",  "Age"), bg = c("#EE7733","#5E4FA2","#337538" ,"#33BBEE"), cex = 1, digits = 2)
```

## Redundance analyses - Only ECM

```{r}
rda_env <- rda(vst_mat_clean_ECM ~ pH + Habitat + Condition(Site), data = meta_clean)
anova.cca(rda_env, by = "term", permutations = 9999)

rda_space <- rda(vst_mat_clean_ECM ~ Site, data = meta_clean)
anova.cca(rda_space, permutations = 9999)

rda_age <- rda(vst_mat_clean_ECM ~ Age + Condition(pH + Habitat + Site), data = meta_clean)
anova.cca(rda_age, by = "term", permutations = 9999)
```

## Mantel test (all OTUs)

```{r Mantel tests within plots - All OTUs}
# test whether RAF community dissimilarity increases with the distance between
# trees inside a plot (distance decay), and whether such a pattern is explained
# by differences in soil pH between the trees

plots <- sort(unique(as.character(meta_clean$Plot)))   # AF AT BF BT IF IT

set.seed(50)

# ---------- AF
sel <- which(as.character(meta_clean$Plot) == "AF" &
               !is.na(meta_clean$Latitude) & !is.na(meta_clean$Longitude))
meta_plot <- meta_clean[sel, ]
lat_mean  <- mean(meta_plot$Latitude)
x <- (meta_plot$Longitude - mean(meta_plot$Longitude)) * 111320 * cos(lat_mean * pi / 180)
y <- (meta_plot$Latitude - lat_mean) * 110540
dist_geo <- dist(cbind(x, y))
dist_com <- dist(vst_mat_clean[sel, ])
dist_pH  <- dist(as.matrix(meta_plot$pH))
mt_geo    <- mantel(dist_com, dist_geo, method = "pearson", permutations = 9999)
mt_geo_pH <- mantel.partial(dist_com, dist_geo, dist_pH, method = "pearson", permutations = 9999)
mt_pH_geo <- mantel(dist_pH, dist_geo, method = "pearson", permutations = 9999)
m_af <- data.frame(Plot = "AF", n = length(sel),
                r_geo    = mt_geo$statistic,    p_geo    = mt_geo$signif,
                r_geo_pH = mt_geo_pH$statistic, p_geo_pH = mt_geo_pH$signif,
                r_pH_geo = mt_pH_geo$statistic, p_pH_geo = mt_pH_geo$signif)

# ---------- AT
sel <- which(as.character(meta_clean$Plot) == "AT" &
               !is.na(meta_clean$Latitude) & !is.na(meta_clean$Longitude))
meta_plot <- meta_clean[sel, ]
lat_mean  <- mean(meta_plot$Latitude)
x <- (meta_plot$Longitude - mean(meta_plot$Longitude)) * 111320 * cos(lat_mean * pi / 180)
y <- (meta_plot$Latitude - lat_mean) * 110540
dist_geo <- dist(cbind(x, y))
dist_com <- dist(vst_mat_clean[sel, ])
dist_pH  <- dist(as.matrix(meta_plot$pH))
mt_geo    <- mantel(dist_com, dist_geo, method = "pearson", permutations = 9999)
mt_geo_pH <- mantel.partial(dist_com, dist_geo, dist_pH, method = "pearson", permutations = 9999)
mt_pH_geo <- mantel(dist_pH, dist_geo, method = "pearson", permutations = 9999)
m_at <- data.frame(Plot = "AT", n = length(sel),
                r_geo    = mt_geo$statistic,    p_geo    = mt_geo$signif,
                r_geo_pH = mt_geo_pH$statistic, p_geo_pH = mt_geo_pH$signif,
                r_pH_geo = mt_pH_geo$statistic, p_pH_geo = mt_pH_geo$signif)

# ---------- BF
sel <- which(as.character(meta_clean$Plot) == "BF" &
               !is.na(meta_clean$Latitude) & !is.na(meta_clean$Longitude))
meta_plot <- meta_clean[sel, ]
lat_mean  <- mean(meta_plot$Latitude)
x <- (meta_plot$Longitude - mean(meta_plot$Longitude)) * 111320 * cos(lat_mean * pi / 180)
y <- (meta_plot$Latitude - lat_mean) * 110540
dist_geo <- dist(cbind(x, y))
dist_com <- dist(vst_mat_clean[sel, ])
dist_pH  <- dist(as.matrix(meta_plot$pH))
mt_geo    <- mantel(dist_com, dist_geo, method = "pearson", permutations = 9999)
mt_geo_pH <- mantel.partial(dist_com, dist_geo, dist_pH, method = "pearson", permutations = 9999)
mt_pH_geo <- mantel(dist_pH, dist_geo, method = "pearson", permutations = 9999)
m_bf <- data.frame(Plot = "BF", n = length(sel),
                r_geo    = mt_geo$statistic,    p_geo    = mt_geo$signif,
                r_geo_pH = mt_geo_pH$statistic, p_geo_pH = mt_geo_pH$signif,
                r_pH_geo = mt_pH_geo$statistic, p_pH_geo = mt_pH_geo$signif)

# ---------- BT
sel <- which(as.character(meta_clean$Plot) == "BT" &
               !is.na(meta_clean$Latitude) & !is.na(meta_clean$Longitude))
meta_plot <- meta_clean[sel, ]
lat_mean  <- mean(meta_plot$Latitude)
x <- (meta_plot$Longitude - mean(meta_plot$Longitude)) * 111320 * cos(lat_mean * pi / 180)
y <- (meta_plot$Latitude - lat_mean) * 110540
dist_geo <- dist(cbind(x, y))
dist_com <- dist(vst_mat_clean[sel, ])
dist_pH  <- dist(as.matrix(meta_plot$pH))
mt_geo    <- mantel(dist_com, dist_geo, method = "pearson", permutations = 9999)
mt_geo_pH <- mantel.partial(dist_com, dist_geo, dist_pH, method = "pearson", permutations = 9999)
mt_pH_geo <- mantel(dist_pH, dist_geo, method = "pearson", permutations = 9999)
m_bt <- data.frame(Plot = "BT", n = length(sel),
                r_geo    = mt_geo$statistic,    p_geo    = mt_geo$signif,
                r_geo_pH = mt_geo_pH$statistic, p_geo_pH = mt_geo_pH$signif,
                r_pH_geo = mt_pH_geo$statistic, p_pH_geo = mt_pH_geo$signif)

# ---------- IF
sel <- which(as.character(meta_clean$Plot) == "IF" &
               !is.na(meta_clean$Latitude) & !is.na(meta_clean$Longitude))
meta_plot <- meta_clean[sel, ]
lat_mean  <- mean(meta_plot$Latitude)
x <- (meta_plot$Longitude - mean(meta_plot$Longitude)) * 111320 * cos(lat_mean * pi / 180)
y <- (meta_plot$Latitude - lat_mean) * 110540
dist_geo <- dist(cbind(x, y))
dist_com <- dist(vst_mat_clean[sel, ])
dist_pH  <- dist(as.matrix(meta_plot$pH))
mt_geo    <- mantel(dist_com, dist_geo, method = "pearson", permutations = 9999)
mt_geo_pH <- mantel.partial(dist_com, dist_geo, dist_pH, method = "pearson", permutations = 9999)
mt_pH_geo <- mantel(dist_pH, dist_geo, method = "pearson", permutations = 9999)
m_if <- data.frame(Plot = "IF", n = length(sel),
                r_geo    = mt_geo$statistic,    p_geo    = mt_geo$signif,
                r_geo_pH = mt_geo_pH$statistic, p_geo_pH = mt_geo_pH$signif,
                r_pH_geo = mt_pH_geo$statistic, p_pH_geo = mt_pH_geo$signif)

# ---------- IT
sel <- which(as.character(meta_clean$Plot) == "IT" &
               !is.na(meta_clean$Latitude) & !is.na(meta_clean$Longitude))
meta_plot <- meta_clean[sel, ]
lat_mean  <- mean(meta_plot$Latitude)
x <- (meta_plot$Longitude - mean(meta_plot$Longitude)) * 111320 * cos(lat_mean * pi / 180)
y <- (meta_plot$Latitude - lat_mean) * 110540
dist_geo <- dist(cbind(x, y))
dist_com <- dist(vst_mat_clean[sel, ])
dist_pH  <- dist(as.matrix(meta_plot$pH))
mt_geo    <- mantel(dist_com, dist_geo, method = "pearson", permutations = 9999)
mt_geo_pH <- mantel.partial(dist_com, dist_geo, dist_pH, method = "pearson", permutations = 9999)
mt_pH_geo <- mantel(dist_pH, dist_geo, method = "pearson", permutations = 9999)
m_it <- data.frame(Plot = "IT", n = length(sel),
                r_geo    = mt_geo$statistic,    p_geo    = mt_geo$signif,
                r_geo_pH = mt_geo_pH$statistic, p_geo_pH = mt_geo_pH$signif,
                r_pH_geo = mt_pH_geo$statistic, p_pH_geo = mt_pH_geo$signif)

mantel_res <- rbind(m_af, m_at, m_bf, m_bt, m_if, m_it)

# adjust the p-values across the six plots
mantel_res$p_geo_adj    <- p.adjust(mantel_res$p_geo,    method = "BH")
mantel_res$p_geo_pH_adj <- p.adjust(mantel_res$p_geo_pH, method = "BH")
mantel_res$p_pH_geo_adj <- p.adjust(mantel_res$p_pH_geo, method = "BH")

mantel_res
```

## Mantel test (ECM OTUs)

```{r Mantel tests within plots - ECM}
# test whether RAF community dissimilarity increases with the distance between
# trees inside a plot (distance decay), and whether such a pattern is explained
# by differences in soil pH between the trees

plots <- sort(unique(as.character(meta_clean$Plot)))   # AF AT BF BT IF IT

set.seed(50)

# ---------- AF
sel <- which(as.character(meta_clean$Plot) == "AF" &
               !is.na(meta_clean$Latitude) & !is.na(meta_clean$Longitude))
meta_plot <- meta_clean[sel, ]
lat_mean  <- mean(meta_plot$Latitude)
x <- (meta_plot$Longitude - mean(meta_plot$Longitude)) * 111320 * cos(lat_mean * pi / 180)
y <- (meta_plot$Latitude - lat_mean) * 110540
dist_geo <- dist(cbind(x, y))
dist_com <- dist(vst_mat_clean_ECM[sel, ])
dist_pH  <- dist(as.matrix(meta_plot$pH))
mt_geo    <- mantel(dist_com, dist_geo, method = "pearson", permutations = 9999)
mt_geo_pH <- mantel.partial(dist_com, dist_geo, dist_pH, method = "pearson", permutations = 9999)
mt_pH_geo <- mantel(dist_pH, dist_geo, method = "pearson", permutations = 9999)
m_af <- data.frame(Plot = "AF", n = length(sel),
                r_geo    = mt_geo$statistic,    p_geo    = mt_geo$signif,
                r_geo_pH = mt_geo_pH$statistic, p_geo_pH = mt_geo_pH$signif,
                r_pH_geo = mt_pH_geo$statistic, p_pH_geo = mt_pH_geo$signif)

# ---------- AT
sel <- which(as.character(meta_clean$Plot) == "AT" &
               !is.na(meta_clean$Latitude) & !is.na(meta_clean$Longitude))
meta_plot <- meta_clean[sel, ]
lat_mean  <- mean(meta_plot$Latitude)
x <- (meta_plot$Longitude - mean(meta_plot$Longitude)) * 111320 * cos(lat_mean * pi / 180)
y <- (meta_plot$Latitude - lat_mean) * 110540
dist_geo <- dist(cbind(x, y))
dist_com <- dist(vst_mat_clean_ECM[sel, ])
dist_pH  <- dist(as.matrix(meta_plot$pH))
mt_geo    <- mantel(dist_com, dist_geo, method = "pearson", permutations = 9999)
mt_geo_pH <- mantel.partial(dist_com, dist_geo, dist_pH, method = "pearson", permutations = 9999)
mt_pH_geo <- mantel(dist_pH, dist_geo, method = "pearson", permutations = 9999)
m_at <- data.frame(Plot = "AT", n = length(sel),
                r_geo    = mt_geo$statistic,    p_geo    = mt_geo$signif,
                r_geo_pH = mt_geo_pH$statistic, p_geo_pH = mt_geo_pH$signif,
                r_pH_geo = mt_pH_geo$statistic, p_pH_geo = mt_pH_geo$signif)

# ---------- BF
sel <- which(as.character(meta_clean$Plot) == "BF" &
               !is.na(meta_clean$Latitude) & !is.na(meta_clean$Longitude))
meta_plot <- meta_clean[sel, ]
lat_mean  <- mean(meta_plot$Latitude)
x <- (meta_plot$Longitude - mean(meta_plot$Longitude)) * 111320 * cos(lat_mean * pi / 180)
y <- (meta_plot$Latitude - lat_mean) * 110540
dist_geo <- dist(cbind(x, y))
dist_com <- dist(vst_mat_clean_ECM[sel, ])
dist_pH  <- dist(as.matrix(meta_plot$pH))
mt_geo    <- mantel(dist_com, dist_geo, method = "pearson", permutations = 9999)
mt_geo_pH <- mantel.partial(dist_com, dist_geo, dist_pH, method = "pearson", permutations = 9999)
mt_pH_geo <- mantel(dist_pH, dist_geo, method = "pearson", permutations = 9999)
m_bf <- data.frame(Plot = "BF", n = length(sel),
                r_geo    = mt_geo$statistic,    p_geo    = mt_geo$signif,
                r_geo_pH = mt_geo_pH$statistic, p_geo_pH = mt_geo_pH$signif,
                r_pH_geo = mt_pH_geo$statistic, p_pH_geo = mt_pH_geo$signif)

# ---------- BT
sel <- which(as.character(meta_clean$Plot) == "BT" &
               !is.na(meta_clean$Latitude) & !is.na(meta_clean$Longitude))
meta_plot <- meta_clean[sel, ]
lat_mean  <- mean(meta_plot$Latitude)
x <- (meta_plot$Longitude - mean(meta_plot$Longitude)) * 111320 * cos(lat_mean * pi / 180)
y <- (meta_plot$Latitude - lat_mean) * 110540
dist_geo <- dist(cbind(x, y))
dist_com <- dist(vst_mat_clean_ECM[sel, ])
dist_pH  <- dist(as.matrix(meta_plot$pH))
mt_geo    <- mantel(dist_com, dist_geo, method = "pearson", permutations = 9999)
mt_geo_pH <- mantel.partial(dist_com, dist_geo, dist_pH, method = "pearson", permutations = 9999)
mt_pH_geo <- mantel(dist_pH, dist_geo, method = "pearson", permutations = 9999)
m_bt <- data.frame(Plot = "BT", n = length(sel),
                r_geo    = mt_geo$statistic,    p_geo    = mt_geo$signif,
                r_geo_pH = mt_geo_pH$statistic, p_geo_pH = mt_geo_pH$signif,
                r_pH_geo = mt_pH_geo$statistic, p_pH_geo = mt_pH_geo$signif)

# ---------- IF
sel <- which(as.character(meta_clean$Plot) == "IF" &
               !is.na(meta_clean$Latitude) & !is.na(meta_clean$Longitude))
meta_plot <- meta_clean[sel, ]
lat_mean  <- mean(meta_plot$Latitude)
x <- (meta_plot$Longitude - mean(meta_plot$Longitude)) * 111320 * cos(lat_mean * pi / 180)
y <- (meta_plot$Latitude - lat_mean) * 110540
dist_geo <- dist(cbind(x, y))
dist_com <- dist(vst_mat_clean_ECM[sel, ])
dist_pH  <- dist(as.matrix(meta_plot$pH))
mt_geo    <- mantel(dist_com, dist_geo, method = "pearson", permutations = 9999)
mt_geo_pH <- mantel.partial(dist_com, dist_geo, dist_pH, method = "pearson", permutations = 9999)
mt_pH_geo <- mantel(dist_pH, dist_geo, method = "pearson", permutations = 9999)
m_if <- data.frame(Plot = "IF", n = length(sel),
                r_geo    = mt_geo$statistic,    p_geo    = mt_geo$signif,
                r_geo_pH = mt_geo_pH$statistic, p_geo_pH = mt_geo_pH$signif,
                r_pH_geo = mt_pH_geo$statistic, p_pH_geo = mt_pH_geo$signif)

# ---------- IT
sel <- which(as.character(meta_clean$Plot) == "IT" &
               !is.na(meta_clean$Latitude) & !is.na(meta_clean$Longitude))
meta_plot <- meta_clean[sel, ]
lat_mean  <- mean(meta_plot$Latitude)
x <- (meta_plot$Longitude - mean(meta_plot$Longitude)) * 111320 * cos(lat_mean * pi / 180)
y <- (meta_plot$Latitude - lat_mean) * 110540
dist_geo <- dist(cbind(x, y))
dist_com <- dist(vst_mat_clean_ECM[sel, ])
dist_pH  <- dist(as.matrix(meta_plot$pH))
mt_geo    <- mantel(dist_com, dist_geo, method = "pearson", permutations = 9999)
mt_geo_pH <- mantel.partial(dist_com, dist_geo, dist_pH, method = "pearson", permutations = 9999)
mt_pH_geo <- mantel(dist_pH, dist_geo, method = "pearson", permutations = 9999)
m_it <- data.frame(Plot = "IT", n = length(sel),
                r_geo    = mt_geo$statistic,    p_geo    = mt_geo$signif,
                r_geo_pH = mt_geo_pH$statistic, p_geo_pH = mt_geo_pH$signif,
                r_pH_geo = mt_pH_geo$statistic, p_pH_geo = mt_pH_geo$signif)

mantel_res_ECM <- rbind(m_af, m_at, m_bf, m_bt, m_if, m_it)

# adjust the p-values across the six plots
mantel_res_ECM$p_geo_adj    <- p.adjust(mantel_res_ECM$p_geo,    method = "BH")
mantel_res_ECM$p_geo_pH_adj <- p.adjust(mantel_res_ECM$p_geo_pH, method = "BH")
mantel_res_ECM$p_pH_geo_adj <- p.adjust(mantel_res_ECM$p_pH_geo, method = "BH")

mantel_res_ECM
```


## Prediction of Tree growth (rcsBAI) - all OTUs 

```{r}
meta2 <- as(sample_data(vst_physeq), "data.frame")
keep2 <- complete.cases(meta2[, c("Age","pH","Habitat","Site", "rcsBAI_5y", "rcsBAI_10y","rcsBAI_15y","rcsBAI_30y")])
vst_mat <- as.data.frame(t(otu_table(vst_physeq)))
vst_mat_clean2 <- vst_mat[keep2, ]
meta_clean2 <- meta2[keep2, ]

meta_clean2$rcsBAI_5y_trans <- log(meta_clean2$rcsBAI_5y)
meta_clean2$rcsBAI_10y_trans <- log(meta_clean2$rcsBAI_10y)
meta_clean2$rcsBAI_15y_trans <- log(meta_clean2$rcsBAI_15y)
meta_clean2$rcsBAI_30y_trans <- log(meta_clean2$rcsBAI_30y)

# PCA on OTU table
pca <- prcomp(vst_mat_clean2, scale. = TRUE)
pcs_df <- as.data.frame(pca$x[, 1:5]) 

# Add BAI and environmental predictors
pcs_df$rcsBAI_5y_trans <- meta_clean2$rcsBAI_5y_trans
pcs_df$rcsBAI_10y_trans <- meta_clean2$rcsBAI_10y_trans
pcs_df$rcsBAI_15y_trans <- meta_clean2$rcsBAI_15y_trans
pcs_df$rcsBAI_30y_trans <- meta_clean2$rcsBAI_30y_trans
pcs_df$pH <- meta_clean2$pH
pcs_df$Habitat <- meta_clean2$Habitat
pcs_df$Site <- meta_clean2$Site

#  RDA using PCs + environmental variables
rda_comm <-  rda(pcs_df$rcsBAI_5y_trans ~ PC1 + PC2 + PC3 + PC4 + PC5 + Condition(Site + Habitat + pH),data = pcs_df)
anova(rda_comm, permutations = 9999)
RsquareAdj(rda_comm)

rda_comm <- rda(pcs_df$rcsBAI_10y_trans ~ PC1 + PC2 + PC3 + PC4 + PC5 + Condition(Site + Habitat + pH),data = pcs_df)
anova(rda_comm, permutations = 9999)
RsquareAdj(rda_comm)

rda_comm <- rda(pcs_df$rcsBAI_15y_trans ~ PC1 + PC2 + PC3 + PC4 + PC5 + Condition(Site + Habitat + pH),data = pcs_df)
anova(rda_comm, permutations = 9999)
RsquareAdj(rda_comm)

rda_comm <- rda(pcs_df$rcsBAI_30y_trans ~ PC1 + PC2 + PC3 + PC4 + PC5 + Condition(Site + Habitat + pH),data = pcs_df)
anova(rda_comm, permutations = 9999)
RsquareAdj(rda_comm)
```

## Prediction of Tree growth (rcsBAI) - only ECM

```{r}
tax_data <- tax_tab
tax_data_ECM <- tax_data[tax_data$Guild == "Ectomycorrhizal", ]
vst_mat_clean2_ECM <- vst_mat_clean2[, colnames(vst_mat_clean2) %in% tax_data_ECM$OTU]

# PCA on OTU table
pca_ECM <- prcomp(vst_mat_clean2_ECM, scale. = TRUE)
pcs_df_ECM <- as.data.frame(pca_ECM$x[, 1:5])

# Add BAI and environmental predictors
pcs_df_ECM$rcsBAI_5y_trans <- meta_clean2$rcsBAI_5y_trans
pcs_df_ECM$rcsBAI_10y_trans <- meta_clean2$rcsBAI_10y_trans
pcs_df_ECM$rcsBAI_15y_trans <- meta_clean2$rcsBAI_15y_trans
pcs_df_ECM$rcsBAI_30y_trans <- meta_clean2$rcsBAI_30y_trans
pcs_df_ECM$pH <- meta_clean2$pH
pcs_df_ECM$Habitat <- meta_clean2$Habitat
pcs_df_ECM$Site <- meta_clean2$Site

#  RDA using PCs + environmental variables
rda_comm <- rda(pcs_df_ECM$rcsBAI_5y_trans ~ PC1 + PC2 + PC3 + PC4 + PC5 + Condition(Site + Habitat + pH),data = pcs_df_ECM)
anova(rda_comm, permutations = 9999)
RsquareAdj(rda_comm)

rda_comm <- rda(pcs_df_ECM$rcsBAI_10y_trans ~ PC1 + PC2 + PC3 + PC4 + PC5 + Condition(Site + Habitat + pH),data = pcs_df_ECM)
anova(rda_comm, permutations = 9999)
RsquareAdj(rda_comm)

rda_comm <- rda(pcs_df_ECM$rcsBAI_15y_trans ~ PC1 + PC2 + PC3 + PC4 + PC5 + Condition(Site + Habitat + pH),data = pcs_df_ECM)
anova(rda_comm, permutations = 9999)
RsquareAdj(rda_comm)

rda_comm <- rda(pcs_df_ECM$rcsBAI_30y_trans ~ PC1 + PC2 + PC3 + PC4 + PC5 + Condition(Site + Habitat + pH),data = pcs_df_ECM)
anova(rda_comm, permutations = 9999)
RsquareAdj(rda_comm)
```

## OTU abundancy per Taxa per Site - treeline vs forest (age-adjusted): Genus

```{r warning=FALSE}

# genus-level count matrix
rank <- "Genus"
tax_map <- as.data.frame(tax_tab[, c("OTU", rank)])
names(tax_map) <- c("OTU", "taxon")
tax_map <- tax_map[tax_map$taxon != "-" & !is.na(tax_map$taxon), ]
tax_map <- tax_map[tax_map$OTU %in% rownames(filt_count_tab_no_r), ]
taxon_counts <- rowsum(as.matrix(filt_count_tab_no_r[tax_map$OTU, ]),
                       group = tax_map$taxon)

# ---------- Alaska_Range
samples <- rownames(sample_info_no_R[sample_info_no_R$Site == "Alaska_Range", ])
samples <- samples[!is.na(sample_info_no_R[samples, "Age"])]
info <- sample_info_no_R[samples, ]
info$Habitat <- factor(info$Habitat, levels = c("forest", "treeline"))
info$Age_s   <- as.numeric(scale(info$Age))

mat    <- taxon_counts[, samples]
mat    <- mat[rowSums(mat > 0) >= 3, ]                          # taxa in at least 3 trees
mat_ab <- mat[rowSums(mat > 0) >= ceiling(0.30 * ncol(mat)), ]  # taxa in at least 30% of trees

dds <- DESeqDataSetFromMatrix(mat_ab, colData = info, design = ~ Age_s + Habitat)
dds <- estimateSizeFactors(dds, type = "poscounts")
res_hab <- results(DESeq(dds, test = "LRT", reduced = ~ Age_s,   fitType = "local", quiet = TRUE),
                   pAdjustMethod = "BH", independentFiltering = FALSE, cooksCutoff = FALSE)
res_age <- results(DESeq(dds, test = "LRT", reduced = ~ Habitat, fitType = "local", quiet = TRUE),
                   pAdjustMethod = "BH", independentFiltering = FALSE, cooksCutoff = FALSE)
da_ar <- data.frame(Site = "Alaska_Range", taxon = rownames(res_hab),
                  log2FC = res_hab$log2FoldChange, p_adj = res_hab$padj,
                  log2FC_age = res_age$log2FoldChange, p_adj_age = res_age$padj)

occ_ar <- data.frame(Site = "Alaska_Range", taxon = rownames(mat),
                   n_forest = NA, n_treeline = NA, p = NA)

# one Fisher test per taxon
for (i in 1:nrow(mat)) {   
  tb <- table(factor(mat[i, ] > 0, levels = c(FALSE, TRUE)), info$Habitat)
  occ_ar$n_forest[i]   <- tb["TRUE", "forest"]
  occ_ar$n_treeline[i] <- tb["TRUE", "treeline"]
  occ_ar$p[i]          <- fisher.test(tb)$p.value}

occ_ar$p_adj <- p.adjust(occ_ar$p, method = "BH")
occ_ar$p     <- NULL

# ---------- Brooks_Range
samples <- rownames(sample_info_no_R[sample_info_no_R$Site == "Brooks_Range", ])
samples <- samples[!is.na(sample_info_no_R[samples, "Age"])]
info <- sample_info_no_R[samples, ]
info$Habitat <- factor(info$Habitat, levels = c("forest", "treeline"))
info$Age_s   <- as.numeric(scale(info$Age))

mat    <- taxon_counts[, samples]
mat    <- mat[rowSums(mat > 0) >= 3, ]                          # taxa in at least 3 trees
mat_ab <- mat[rowSums(mat > 0) >= ceiling(0.30 * ncol(mat)), ]  # taxa in at least 30% of trees

dds <- DESeqDataSetFromMatrix(mat_ab, colData = info, design = ~ Age_s + Habitat)
dds <- estimateSizeFactors(dds, type = "poscounts")
res_hab <- results(DESeq(dds, test = "LRT", reduced = ~ Age_s,   fitType = "local", quiet = TRUE),
                   pAdjustMethod = "BH", independentFiltering = FALSE, cooksCutoff = FALSE)
res_age <- results(DESeq(dds, test = "LRT", reduced = ~ Habitat, fitType = "local", quiet = TRUE),
                   pAdjustMethod = "BH", independentFiltering = FALSE, cooksCutoff = FALSE)
da_br <- data.frame(Site = "Brooks_Range", taxon = rownames(res_hab),
                  log2FC = res_hab$log2FoldChange, p_adj = res_hab$padj,
                  log2FC_age = res_age$log2FoldChange, p_adj_age = res_age$padj)

occ_br <- data.frame(Site = "Brooks_Range", taxon = rownames(mat),
                   n_forest = NA, n_treeline = NA, p = NA)

# one Fisher test per taxon
for (i in 1:nrow(mat)) {   
  tb <- table(factor(mat[i, ] > 0, levels = c(FALSE, TRUE)), info$Habitat)
  occ_br$n_forest[i]   <- tb["TRUE", "forest"]
  occ_br$n_treeline[i] <- tb["TRUE", "treeline"]
  occ_br$p[i]          <- fisher.test(tb)$p.value}

occ_br$p_adj <- p.adjust(occ_br$p, method = "BH")
occ_br$p     <- NULL

# ---------- Interior
samples <- rownames(sample_info_no_R[sample_info_no_R$Site == "Interior", ])
samples <- samples[!is.na(sample_info_no_R[samples, "Age"])]
info <- sample_info_no_R[samples, ]
info$Habitat <- factor(info$Habitat, levels = c("forest", "treeline"))
info$Age_s   <- as.numeric(scale(info$Age))

mat    <- taxon_counts[, samples]
mat    <- mat[rowSums(mat > 0) >= 3, ]                          # taxa in at least 3 trees
mat_ab <- mat[rowSums(mat > 0) >= ceiling(0.30 * ncol(mat)), ]  # taxa in at least 30% of trees

dds <- DESeqDataSetFromMatrix(mat_ab, colData = info, design = ~ Age_s + Habitat)
dds <- estimateSizeFactors(dds, type = "poscounts")
res_hab <- results(DESeq(dds, test = "LRT", reduced = ~ Age_s,   fitType = "local", quiet = TRUE),
                   pAdjustMethod = "BH", independentFiltering = FALSE, cooksCutoff = FALSE)
res_age <- results(DESeq(dds, test = "LRT", reduced = ~ Habitat, fitType = "local", quiet = TRUE),
                   pAdjustMethod = "BH", independentFiltering = FALSE, cooksCutoff = FALSE)
da_in <- data.frame(Site = "Interior", taxon = rownames(res_hab),
                  log2FC = res_hab$log2FoldChange, p_adj = res_hab$padj,
                  log2FC_age = res_age$log2FoldChange, p_adj_age = res_age$padj)

occ_in <- data.frame(Site = "Interior", taxon = rownames(mat),
                   n_forest = NA, n_treeline = NA, p = NA)

# one Fisher test per taxon
for (i in 1:nrow(mat)) {   
  tb <- table(factor(mat[i, ] > 0, levels = c(FALSE, TRUE)), info$Habitat)
  occ_in$n_forest[i]   <- tb["TRUE", "forest"]
  occ_in$n_treeline[i] <- tb["TRUE", "treeline"]
  occ_in$p[i]          <- fisher.test(tb)$p.value}

occ_in$p_adj <- p.adjust(occ_in$p, method = "BH")
occ_in$p     <- NULL

da_all  <- rbind(da_ar, da_br, da_in)
occ_all <- rbind(occ_ar, occ_br, occ_in)

da_all[!is.na(da_all$p_adj)     & da_all$p_adj     < 0.05, ]
occ_all[occ_all$p_adj < 0.05, ]
da_all[!is.na(da_all$p_adj_age) & da_all$p_adj_age < 0.05, ]   # taxa driven by age

res_all <- merge(da_all, occ_all, by = c("Site", "taxon"),
                 all = TRUE, suffixes = c("_abund", "_occ"))
res_all$log2FC      <- signif(res_all$log2FC, 3)      # signif, not round: keeps small p-values
res_all$p_adj_abund <- signif(res_all$p_adj_abund, 3)
res_all$log2FC_age  <- signif(res_all$log2FC_age, 3)
res_all$p_adj_age   <- signif(res_all$p_adj_age, 3)
res_all$p_adj_occ   <- signif(res_all$p_adj_occ, 3)
write.csv(res_all, "Habitat_tests.csv", row.names = FALSE)
```

## OTU abundancy per Taxa per Site - treeline vs forest (age-adjusted): Species

```{r warning=FALSE}

# genus-level count matrix
rank <- "Species...14"
tax_map <- as.data.frame(tax_tab[, c("OTU", rank)])
names(tax_map) <- c("OTU", "taxon")
tax_map <- tax_map[tax_map$taxon != "-" & !is.na(tax_map$taxon), ]
tax_map <- tax_map[tax_map$OTU %in% rownames(filt_count_tab_no_r), ]
taxon_counts <- rowsum(as.matrix(filt_count_tab_no_r[tax_map$OTU, ]),
                       group = tax_map$taxon)

# ---------- Alaska_Range
samples <- rownames(sample_info_no_R[sample_info_no_R$Site == "Alaska_Range", ])
samples <- samples[!is.na(sample_info_no_R[samples, "Age"])]
info <- sample_info_no_R[samples, ]
info$Habitat <- factor(info$Habitat, levels = c("forest", "treeline"))
info$Age_s   <- as.numeric(scale(info$Age))

mat    <- taxon_counts[, samples]
mat    <- mat[rowSums(mat > 0) >= 3, ]                          # taxa in at least 3 trees
mat_ab <- mat[rowSums(mat > 0) >= ceiling(0.30 * ncol(mat)), ]  # taxa in at least 30% of trees

dds <- DESeqDataSetFromMatrix(mat_ab, colData = info, design = ~ Age_s + Habitat)
dds <- estimateSizeFactors(dds, type = "poscounts")
res_hab <- results(DESeq(dds, test = "LRT", reduced = ~ Age_s,   fitType = "local", quiet = TRUE),
                   pAdjustMethod = "BH", independentFiltering = FALSE, cooksCutoff = FALSE)
res_age <- results(DESeq(dds, test = "LRT", reduced = ~ Habitat, fitType = "local", quiet = TRUE),
                   pAdjustMethod = "BH", independentFiltering = FALSE, cooksCutoff = FALSE)
da_ar <- data.frame(Site = "Alaska_Range", taxon = rownames(res_hab),
                  log2FC = res_hab$log2FoldChange, p_adj = res_hab$padj,
                  log2FC_age = res_age$log2FoldChange, p_adj_age = res_age$padj)

occ_ar <- data.frame(Site = "Alaska_Range", taxon = rownames(mat),
                   n_forest = NA, n_treeline = NA, p = NA)

# one Fisher test per taxon
for (i in 1:nrow(mat)) {   
  tb <- table(factor(mat[i, ] > 0, levels = c(FALSE, TRUE)), info$Habitat)
  occ_ar$n_forest[i]   <- tb["TRUE", "forest"]
  occ_ar$n_treeline[i] <- tb["TRUE", "treeline"]
  occ_ar$p[i]          <- fisher.test(tb)$p.value}

occ_ar$p_adj <- p.adjust(occ_ar$p, method = "BH")
occ_ar$p     <- NULL

# ---------- Brooks_Range
samples <- rownames(sample_info_no_R[sample_info_no_R$Site == "Brooks_Range", ])
samples <- samples[!is.na(sample_info_no_R[samples, "Age"])]
info <- sample_info_no_R[samples, ]
info$Habitat <- factor(info$Habitat, levels = c("forest", "treeline"))
info$Age_s   <- as.numeric(scale(info$Age))

mat    <- taxon_counts[, samples]
mat    <- mat[rowSums(mat > 0) >= 3, ]                          # taxa in at least 3 trees
mat_ab <- mat[rowSums(mat > 0) >= ceiling(0.30 * ncol(mat)), ]  # taxa in at least 30% of trees

dds <- DESeqDataSetFromMatrix(mat_ab, colData = info, design = ~ Age_s + Habitat)
dds <- estimateSizeFactors(dds, type = "poscounts")
res_hab <- results(DESeq(dds, test = "LRT", reduced = ~ Age_s,   fitType = "local", quiet = TRUE),
                   pAdjustMethod = "BH", independentFiltering = FALSE, cooksCutoff = FALSE)
res_age <- results(DESeq(dds, test = "LRT", reduced = ~ Habitat, fitType = "local", quiet = TRUE),
                   pAdjustMethod = "BH", independentFiltering = FALSE, cooksCutoff = FALSE)
da_br <- data.frame(Site = "Brooks_Range", taxon = rownames(res_hab),
                  log2FC = res_hab$log2FoldChange, p_adj = res_hab$padj,
                  log2FC_age = res_age$log2FoldChange, p_adj_age = res_age$padj)

occ_br <- data.frame(Site = "Brooks_Range", taxon = rownames(mat),
                   n_forest = NA, n_treeline = NA, p = NA)

# one Fisher test per taxon
for (i in 1:nrow(mat)) {  
  tb <- table(factor(mat[i, ] > 0, levels = c(FALSE, TRUE)), info$Habitat)
  occ_br$n_forest[i]   <- tb["TRUE", "forest"]
  occ_br$n_treeline[i] <- tb["TRUE", "treeline"]
  occ_br$p[i]          <- fisher.test(tb)$p.value}

occ_br$p_adj <- p.adjust(occ_br$p, method = "BH")
occ_br$p     <- NULL

# ---------- Interior
samples <- rownames(sample_info_no_R[sample_info_no_R$Site == "Interior", ])
samples <- samples[!is.na(sample_info_no_R[samples, "Age"])]
info <- sample_info_no_R[samples, ]
info$Habitat <- factor(info$Habitat, levels = c("forest", "treeline"))
info$Age_s   <- as.numeric(scale(info$Age))

mat    <- taxon_counts[, samples]
mat    <- mat[rowSums(mat > 0) >= 3, ]                          # taxa in at least 3 trees
mat_ab <- mat[rowSums(mat > 0) >= ceiling(0.30 * ncol(mat)), ]  # taxa in at least 30% of trees

dds <- DESeqDataSetFromMatrix(mat_ab, colData = info, design = ~ Age_s + Habitat)
dds <- estimateSizeFactors(dds, type = "poscounts")
res_hab <- results(DESeq(dds, test = "LRT", reduced = ~ Age_s,   fitType = "local", quiet = TRUE),
                   pAdjustMethod = "BH", independentFiltering = FALSE, cooksCutoff = FALSE)
res_age <- results(DESeq(dds, test = "LRT", reduced = ~ Habitat, fitType = "local", quiet = TRUE),
                   pAdjustMethod = "BH", independentFiltering = FALSE, cooksCutoff = FALSE)
da_in <- data.frame(Site = "Interior", taxon = rownames(res_hab),
                  log2FC = res_hab$log2FoldChange, p_adj = res_hab$padj,
                  log2FC_age = res_age$log2FoldChange, p_adj_age = res_age$padj)

occ_in <- data.frame(Site = "Interior", taxon = rownames(mat),
                   n_forest = NA, n_treeline = NA, p = NA)

# one Fisher test per taxon
for (i in 1:nrow(mat)) {   
  tb <- table(factor(mat[i, ] > 0, levels = c(FALSE, TRUE)), info$Habitat)
  occ_in$n_forest[i]   <- tb["TRUE", "forest"]
  occ_in$n_treeline[i] <- tb["TRUE", "treeline"]
  occ_in$p[i]          <- fisher.test(tb)$p.value}

occ_in$p_adj <- p.adjust(occ_in$p, method = "BH")
occ_in$p     <- NULL

da_all  <- rbind(da_ar, da_br, da_in)
occ_all <- rbind(occ_ar, occ_br, occ_in)

da_all[!is.na(da_all$p_adj)     & da_all$p_adj     < 0.05, ]
occ_all[occ_all$p_adj < 0.05, ]
da_all[!is.na(da_all$p_adj_age) & da_all$p_adj_age < 0.05, ]   # taxa driven by age

res_all <- merge(da_all, occ_all, by = c("Site", "taxon"),
                 all = TRUE, suffixes = c("_abund", "_occ"))
res_all$log2FC      <- signif(res_all$log2FC, 3)      # signif, not round: keeps small p-values
res_all$p_adj_abund <- signif(res_all$p_adj_abund, 3)
res_all$log2FC_age  <- signif(res_all$log2FC_age, 3)
res_all$p_adj_age   <- signif(res_all$p_adj_age, 3)
res_all$p_adj_occ   <- signif(res_all$p_adj_occ, 3)
write.csv(res_all, "Habitat_tests_species.csv", row.names = FALSE)
```


## OTU abundancy per Taxa per Site

```{r}
tax_data <-tax_tab
tax_data <- tax_data[tax_data$OTU %in% rownames(filt_count_tab_no_r), ] # I used row counts
counts_tree <- as.data.frame(filt_count_tab_no_r)

counts_tree$OTU<-rownames(counts_tree)
rownames(counts_tree) <- NULL

# Combine tables sorting by OTU
Tax_otu_data <- left_join(tax_data, counts_tree,  by = "OTU" ) 
Sites <- unique(sample_info_no_R$Site)

set.seed(50)

samples <- rownames(sample_info_no_R[sample_info_no_R$Site == "Alaska_Range", ])
obj_ar <- parse_tax_data(Tax_otu_data, class_cols = "taxonomy", class_sep = ";")
obj_ar <- filter_taxa(obj_ar, taxon_names != "-")
obj_ar$data$tax_abund <- calc_taxon_abund(obj_ar, data = "tax_data", cols = samples)
obj_ar$data$tax_abund$total_counts <- rowSums(obj_ar$data$tax_abund[, samples])

samples <- rownames(sample_info_no_R[sample_info_no_R$Site == "Brooks_Range", ])
obj_br <- parse_tax_data(Tax_otu_data, class_cols = "taxonomy", class_sep = ";")
obj_br <- filter_taxa(obj_br, taxon_names != "-")
obj_br$data$tax_abund <- calc_taxon_abund(obj_br, data = "tax_data", cols = samples)
obj_br$data$tax_abund$total_counts <- rowSums(obj_br$data$tax_abund[, samples])

samples <- rownames(sample_info_no_R[sample_info_no_R$Site == "Interior", ])
obj_in <- parse_tax_data(Tax_otu_data, class_cols = "taxonomy", class_sep = ";")
obj_in <- filter_taxa(obj_in, taxon_names != "-")
obj_in$data$tax_abund <- calc_taxon_abund(obj_in, data = "tax_data", cols = samples)
obj_in$data$tax_abund$total_counts <- rowSums(obj_in$data$tax_abund[, samples])

results <- list(Alaska_Range = obj_ar, Brooks_Range = obj_br, Interior = obj_in)

obj_loc <- results[["Alaska_Range"]]
obj_loc_D <- filter_taxa(obj_loc, n_obs > 1)
set.seed(50)

hD<-heat_tree(obj_loc_D, node_label = gsub(pattern = "\\[|\\]", replacement = "", taxon_names),
            node_size = log(total_counts),
            node_size_range = c(0.0005, 0.025), 
            node_color = log(total_counts), # one taxa with 0 counts
            node_color_range = c("white","#33BBEE","#EE7733"),
            node_label_size_range = c(0.01, 0.025), 
            node_color_axis_label = "Abundance log(count number)",
            layout = "davidson-harel", initial_layout = "reingold-tilford")
plot(hD)
```

```{r}
obj_loc <- results[["Brooks_Range"]]
obj_loc_B <- filter_taxa(obj_loc, n_obs > 1)

set.seed(50)

hB<-heat_tree(obj_loc_B, node_label = gsub(pattern = "\\[|\\]", replacement = "", taxon_names),
            node_size = log(total_counts),
            node_size_range = c(0.0005, 0.025), 
            node_color = log(total_counts), # no zeros
            node_color_range = c("white","#33BBEE","#EE7733"),
            node_label_size_range = c(0.01, 0.025), 
            node_color_axis_label = "Abundance log(count number)",
            layout = "davidson-harel", initial_layout = "reingold-tilford")
plot(hB)
```

```{r}
obj_loc <- results[["Interior"]]
obj_loc_I <- filter_taxa(obj_loc, n_obs > 1)

set.seed(50)

hI<-heat_tree(obj_loc_I, node_label = gsub(pattern = "\\[|\\]", replacement = "", taxon_names),
            node_size = log(total_counts),
            node_size_range = c(0.0005, 0.025), 
            node_color = log(total_counts), #three taxa with 0 counts
            node_color_range = c("white","#33BBEE","#EE7733"),
            node_label_size_range = c(0.01, 0.025), 
            node_color_axis_label = "Abundance log(count number)",
            layout = "davidson-harel", initial_layout = "reingold-tilford")
plot(hI)
```


## Species diversty (All OTUs): treeline vs forest

```{r}
#Hill number = effective number of species.  Hill numbers determines the measures’ sensitivity to species relative abundances.  Hill numbers include the three most widely used species diversity measures as special cases: species richness (q = 0), Shannon diversity (q = 1, as the effective number of common species in the assemblage) and Simpson diversity (q = 2, effective number of dominant species in the assemblage)

#I used row counts of reads after removing of OTUs found in NTCs and technical replicates

count_tab <-filt_count_tab_no_r
count_tab <- as.data.frame(count_tab > 0) # to make presence/absence matrix

fD_info_rows <- sample_info_no_R[sample_info_no_R$Habitat == "forest" & sample_info_no_R$Site == "Alaska_Range",]
fD_counts_col <- count_tab[colnames(count_tab) %in% rownames(fD_info_rows) ]
tD_info_rows <- sample_info_no_R[sample_info_no_R$Habitat == "treeline" & sample_info_no_R$Site == "Alaska_Range", ]
tD_counts_col <- count_tab[colnames(count_tab) %in% rownames(tD_info_rows) ]

fI_info_rows <- sample_info_no_R[sample_info_no_R$Habitat == "forest" & sample_info_no_R$Site == "Interior",]
fI_counts_col <- count_tab[colnames(count_tab) %in% rownames(fI_info_rows) ]
tI_info_rows <- sample_info_no_R[sample_info_no_R$Habitat == "treeline" & sample_info_no_R$Site == "Interior", ]
tI_counts_col <- count_tab[colnames(count_tab) %in% rownames(tI_info_rows) ]

fB_info_rows <- sample_info_no_R[sample_info_no_R$Habitat == "forest" & sample_info_no_R$Site == "Brooks_Range",]
fB_counts_col <- count_tab[colnames(count_tab) %in% rownames(fB_info_rows) ]
tB_info_rows <- sample_info_no_R[sample_info_no_R$Habitat == "treeline" & sample_info_no_R$Site == "Brooks_Range", ]
tB_counts_col <- count_tab[colnames(count_tab) %in% rownames(tB_info_rows) ]


# create a dataframe with incidence data summed by rows within plots
inc_bytype<- cbind(rowSums(fD_counts_col), rowSums(tD_counts_col), rowSums(fI_counts_col), rowSums(tI_counts_col), rowSums(fB_counts_col), rowSums(tB_counts_col))

colnames(inc_bytype) <- c("A_Forest", "A_Treeline", "I_Forest", "I_Treeline","B_Forest", "B_Treeline")
inc_bytype<-as.data.frame(inc_bytype)
n_samples <-c(ncol(fD_counts_col), ncol(tD_counts_col), ncol(fI_counts_col), ncol(tI_counts_col), ncol(fB_counts_col), ncol(tB_counts_col))
rownames(inc_bytype)<-tax_data$OTU
inc_bytype <-rbind(n_samples, inc_bytype)

# Get the diversity estimates
out <- iNEXT(inc_bytype, q=c(0,1,2), datatype ="incidence_freq", se = TRUE)

# Plot the diversity estimates
# Sample-size-based rarefaction and extrapolation (R/E) curve (type=1). This curve plots diversity estimates with 95% confidence intervals (if se=TRUE) as a function of sample size up to double the reference sample size, by default, or a user‐specified endpoint:

ggiNEXT(out, type=1, facet.var="Order.q", se = TRUE, color.var="Assemblage")
```

```{r}
ggiNEXT(out, type=2, facet.var="None", color.var="Assemblage", se = TRUE)
```

## Species diversty (ECM OTUs): treeline vs forest

```{r}
#I used row counts of reads after removing of OTUs found in NTCs and technical replicates
count_tab <-filt_count_tab_no_r

# filter OTU from tax_data:
tax_data     <- tax_tab[tax_tab$OTU %in% rownames(count_tab), ]
tax_data_ECM <- tax_data[tax_data$Guild %in% "Ectomycorrhizal", ]

count_tab <- count_tab[rownames(count_tab) %in% tax_data_ECM$OTU, ]
count_tab <- as.data.frame(count_tab > 0)   # presence/absence matrix

fD_info_rows <- sample_info_no_R[sample_info_no_R$Habitat == "forest" & sample_info_no_R$Site == "Alaska_Range",]
fD_counts_col <- count_tab[colnames(count_tab) %in% rownames(fD_info_rows) ]
tD_info_rows <- sample_info_no_R[sample_info_no_R$Habitat == "treeline" & sample_info_no_R$Site == "Alaska_Range", ]
tD_counts_col <- count_tab[colnames(count_tab) %in% rownames(tD_info_rows) ]

fI_info_rows <- sample_info_no_R[sample_info_no_R$Habitat == "forest" & sample_info_no_R$Site == "Interior",]
fI_counts_col <- count_tab[colnames(count_tab) %in% rownames(fI_info_rows) ]
tI_info_rows <- sample_info_no_R[sample_info_no_R$Habitat == "treeline" & sample_info_no_R$Site == "Interior", ]
tI_counts_col <- count_tab[colnames(count_tab) %in% rownames(tI_info_rows) ]

fB_info_rows <- sample_info_no_R[sample_info_no_R$Habitat == "forest" & sample_info_no_R$Site == "Brooks_Range",]
fB_counts_col <- count_tab[colnames(count_tab) %in% rownames(fB_info_rows) ]
tB_info_rows <- sample_info_no_R[sample_info_no_R$Habitat == "treeline" & sample_info_no_R$Site == "Brooks_Range", ]
tB_counts_col <- count_tab[colnames(count_tab) %in% rownames(tB_info_rows) ]

# create a dataframe with incidence data summed by rows within plots
inc_bytype<- cbind(rowSums(fD_counts_col), rowSums(tD_counts_col), rowSums(fI_counts_col), rowSums(tI_counts_col), rowSums(fB_counts_col), rowSums(tB_counts_col))

colnames(inc_bytype) <- c("A_Forest", "A_Treeline", "I_Forest", "I_Treeline","B_Forest", "B_Treeline")
inc_bytype<-as.data.frame(inc_bytype)
n_samples <-c(ncol(fD_counts_col), ncol(tD_counts_col), ncol(fI_counts_col), ncol(tI_counts_col), ncol(fB_counts_col), ncol(tB_counts_col))
inc_bytype <-rbind(n_samples, inc_bytype)

# Get the diversity estimates
out <- iNEXT(inc_bytype, q=c(0,1,2), datatype ="incidence_freq", se = TRUE)

# Plot the diversity estimates
# Sample-size-based rarefaction and extrapolation (R/E) curve (type=1). This curve plots diversity estimates with 95% confidence intervals (if se=TRUE) as a function of sample size up to double the reference sample size, by default, or a user‐specified endpoint:

ggiNEXT(out, type=1, facet.var="Order.q", se = TRUE, color.var="Assemblage")
```



