---
title: "Metabarcoding: Trophic Guilds"
author: "Kristina Kuprina"
date: "`r Sys.Date()`"
output: html_document
---

```{r Packages, message=FALSE, warning=FALSE}
library("ggplot2")
library("readr")
library("DESeq2")
library("phyloseq")
library('ggrepel')
library('ggalluvial')
library("microeco")
library("dplyr")
library("tidyr")
library('tibble')
library('pals')
library("rstatix")
library("cowplot")
library("vegan")
library("car")
library("mgcv")
```

```{r session info}
sessionInfo()
```

## Data preparation

```{r sample info}
sample_info_no_R <- read.table("Alaska_info_no_R.txt", row.names=1, header=T, check.names=F) # output of BAI_detBAI_RCS.Rmd: raw BAI, detrended BAI (detBAI) and RCS BAI (rcsBAI)
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
filt_count_tab <- count_tab[!rownames(count_tab) %in% blank_OTUs, -c(1:8)] # removing blank_OTUs and the blank samples from further analysis
filt_count_tab <- filt_count_tab[rowSums(filt_count_tab == 0) != ncol(filt_count_tab), ] # delete OTUs with 0 counts
```

### Normalizing for sampling depth (variance stabilizing transformation) 

```{r Preparing dataset no R}
#reorder the row and column names
sample_info_no_R<-sample_info_no_R[order(rownames(sample_info_no_R)), ]
filt_count_tab<- filt_count_tab[,order(names(filt_count_tab))]
filt_count_tab_no_r <- filt_count_tab[colnames(filt_count_tab) %in% rownames(sample_info_no_R) ]

# create a DESeq2 object
deseq_counts <- DESeqDataSetFromMatrix(filt_count_tab_no_r, colData = sample_info_no_R, design = ~ Habitat) 
deseq_counts <- estimateSizeFactors(deseq_counts,type = "poscounts")
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts,blind = TRUE)
vst_trans_count_tab <- assay(deseq_counts_vst) #save as matrix
vst_trans_count_tab_df0 <- as.data.frame(vst_trans_count_tab)
#transform lowest values (absent counts) to 0
vst_trans_count_tab_df <- sweep(vst_trans_count_tab_df0, 2, apply(vst_trans_count_tab, 2, min), "-") 
```

## Functional guilds relative abundancy

```{r}
tax_table<-as.data.frame(tax_tab)
row.names(tax_table) <- tax_table$OTU
tax_table <- tax_table[, -1]
tax_table<- tax_table[order(rownames(tax_table)),]

# create microeco file
meco_fungi <- microtable$new(sample_table = sample_info_no_R, otu_table = filt_count_tab_no_r, tax_table = tax_table)
mycol <- c("#A6CEE3" ,"#1F78B4" ,"#B2DF8A" ,"#33A02C" ,"#FB9A99", "#E31A1C" ,"#FDBF6F" ,"#FF7F00" ,"#CAB2D6" ,"#6A3D9A" ,"#FFFF99" ,"#B15928" , "#6fdc8c", "green4", "#004144")

t1 <- trans_abund$new(dataset = meco_fungi, taxrank = "Guild", ntaxa = 15)
t1$plot_bar(bar_type = "notfull", use_alluvium = TRUE,xtext_keep = TRUE, xtext_angle = 90,facet = "Plot",xtext_size = 6, color_values = mycol)
```

```{r}
palette <- c("#EE7733","#CC6677", "#009988","#337538","#33BBEE", "#5E4FA2")
t1$plot_box(group = "Plot",  xtext_angle = 90, color_values = palette) + ylab("Relative abundance (%)")
```


## Setup for the guild comparisons

```{r}
abund_rel <- t1$data_abund  # relative abundance data

sapro_guilds <- c("Wood Saprotroph","Plant Saprotroph","Undefined Saprotroph")

# sum the saprotroph guilds within each tree -> one row per tree
sapro_rel_hab <- aggregate(Abundance ~ Sample + Site + Habitat, data = abund_rel[abund_rel$Taxonomy %in% sapro_guilds,], FUN = sum)

anyDuplicated(sapro_rel_hab$Sample) # should be 0
```

## Treeline vs forest - relative abundance data

```{r}
# ECM
ECM_hab_IN <- wilcox_test(abund_rel[abund_rel$Taxonomy == "Ectomycorrhizal" & abund_rel$Site == "Interior",], Abundance ~ Habitat, detailed = TRUE)
ECM_hab_AR <- wilcox_test(abund_rel[abund_rel$Taxonomy == "Ectomycorrhizal" & abund_rel$Site == "Alaska_Range",], Abundance ~ Habitat, detailed = TRUE)
ECM_hab_BR <- wilcox_test(abund_rel[abund_rel$Taxonomy == "Ectomycorrhizal" & abund_rel$Site == "Brooks_Range",], Abundance ~ Habitat, detailed = TRUE)

ECM_hab <- rbind(ECM_hab_IN, ECM_hab_AR, ECM_hab_BR)
ECM_hab$Site <- c("Interior","Alaska_Range","Brooks_Range")
ECM_hab$p.adj <- p.adjust(ECM_hab$p, method = "BH")
ECM_hab

# Saprotrophs
SAP_hab_IN <- wilcox_test(sapro_rel_hab[sapro_rel_hab$Site == "Interior",], Abundance ~ Habitat, detailed = TRUE)
SAP_hab_AR <- wilcox_test(sapro_rel_hab[sapro_rel_hab$Site == "Alaska_Range",], Abundance ~ Habitat, detailed = TRUE)
SAP_hab_BR <- wilcox_test(sapro_rel_hab[sapro_rel_hab$Site == "Brooks_Range",], Abundance ~ Habitat, detailed = TRUE)

SAP_hab <- rbind(SAP_hab_IN, SAP_hab_AR, SAP_hab_BR)
SAP_hab$Site <- c("Interior","Alaska_Range","Brooks_Range")
SAP_hab$p.adj <- p.adjust(SAP_hab$p, method = "BH")
SAP_hab
```


# Guild abundance vs tree growth

## Data preparation

```{r}
# one value per tree, from the relative abundance data
ECM_tree <- abund_rel[abund_rel$Taxonomy == "Ectomycorrhizal", c("Sample","Abundance")]
names(ECM_tree)[2] <- "ECM"

SAP_tree <- aggregate(Abundance ~ Sample, data = abund_rel[abund_rel$Taxonomy %in% sapro_guilds,], FUN = sum)
names(SAP_tree)[2] <- "SAP"

meta_bai <- sample_info_no_R
meta_bai$Sample <- rownames(meta_bai)

guild_bai <- merge(ECM_tree, SAP_tree, by = "Sample")
guild_bai <- merge(guild_bai, meta_bai, by = "Sample")
guild_bai <- guild_bai[is.finite(guild_bai$BAI_5y), ] # remove samples with NA (same 5 trees for all BAI variables)
nrow(guild_bai) # should be 106
```

## Model check: Gaussian vs Gamma (ECM, BAI 10 years)

```{r}
#Check for multicollinearity
vif(lm(log(BAI_10y) ~ ECM + Habitat + Site + log(DBH), data = guild_bai))
vif(lm(log(BAI_10y) ~ SAP + Habitat + Site + log(DBH), data = guild_bai))

# same model, two error distributions
m_gaus <- gam(BAI_10y ~ ECM + Habitat + Site + log(DBH),
  family = gaussian(link = "log"),
  method = "REML",
  data = guild_bai)

m_gamma <- gam(BAI_10y ~ ECM + Habitat + Site + log(DBH),
  family = Gamma(link = "log"),
  method = "REML",
  data = guild_bai)

AIC(m_gaus, m_gamma)  # lower = better
gam.check(m_gaus)     # "Resids vs. linear pred." fans out 
gam.check(m_gamma)

summary(m_gamma)

# prediction line at the median DBH of each plot, so it shows the ECM effect only
nd <- expand.grid(ECM = seq(min(guild_bai$ECM), max(guild_bai$ECM), length.out = 50),
                  Habitat = unique(guild_bai$Habitat),
                  Site = unique(guild_bai$Site))
med_dbh <- aggregate(DBH ~ Habitat + Site, data = guild_bai, FUN = median)
nd <- merge(nd, med_dbh, by = c("Habitat", "Site"))
nd$pred <- predict(m_gamma, newdata = nd, type = "response")

ggplot(guild_bai, aes(x = ECM, y = BAI_10y)) +
  geom_point(alpha = 0.6, color = "#EE7733") +
  geom_line(data = nd, aes(y = pred), color = "#009988", linewidth = 1) +
  facet_grid(Habitat ~ Site) +
  scale_y_log10() +
  theme_bw() + labs(x = "Relative abundance of ECM fungi (%)", y = "Mean BAI over 10 years (cm²/year, log scale)")
```

## Does the ECM effect differ between sites? -> no

```{r}
# both models: same family and same covariates, only the interaction differs
m_ecm <- gam(rcsBAI_10y ~ ECM + Habitat + Site , family = Gamma(link = "log"), method = "REML", data = guild_bai)
m_ecm_int <- gam(rcsBAI_10y ~ ECM * Site + Habitat , family = Gamma(link = "log"), method = "REML", data = guild_bai) #-> no
anova(m_ecm, m_ecm_int, test = "F")

# slope within each site
guild_bai$Site_IN <- relevel(factor(guild_bai$Site), ref = "Interior")
guild_bai$Site_AR <- relevel(factor(guild_bai$Site), ref = "Alaska_Range")
guild_bai$Site_BR <- relevel(factor(guild_bai$Site), ref = "Brooks_Range")

summary(gam(rcsBAI_10y ~ ECM * Site_IN + Habitat , family = Gamma(link = "log"), method = "REML", data = guild_bai))$p.table["ECM",] #-> ns
summary(gam(rcsBAI_10y ~ ECM * Site_AR + Habitat , family = Gamma(link = "log"), method = "REML", data = guild_bai))$p.table["ECM",] #-> ns
summary(gam(rcsBAI_10y ~ ECM * Site_BR + Habitat , family = Gamma(link = "log"), method = "REML", data = guild_bai))$p.table["ECM",] #-> ns
```

## Guild abundance vs raw BAI

```{r}
## Choice of error distribution: Gamma vs Gaussian, both with a log link.
## Positive dAIC = the Gamma fit is better. Checked once here instead of refitting
## a Gaussian model beside every model below.

m_chk_bai <- gam(BAI_10y ~ ECM + Habitat + Site + log(DBH), family = Gamma(link = "log"), method = "REML", data = guild_bai)
g_chk_bai <- gam(BAI_10y ~ ECM + Habitat + Site + log(DBH), family = gaussian(link = "log"), method = "REML", data = guild_bai)

m_chk_det <- gam(detBAI_10y ~ ECM + Habitat + Site, family = Gamma(link = "log"), method = "REML", data = guild_bai)
g_chk_det <- gam(detBAI_10y ~ ECM + Habitat + Site, family = gaussian(link = "log"), method = "REML", data = guild_bai)

m_chk_rcs <- gam(rcsBAI_10y ~ ECM + Plot, family = Gamma(link = "log"), method = "REML", data = guild_bai)
g_chk_rcs <- gam(rcsBAI_10y ~ ECM + Plot, family = gaussian(link = "log"), method = "REML", data = guild_bai)

data.frame(Model = c("BAI_10y", "detBAI_10y", "rcsBAI_10y"),
           dAIC  = c(AIC(g_chk_bai) - AIC(m_chk_bai),
                     AIC(g_chk_det) - AIC(m_chk_det),
                     AIC(g_chk_rcs) - AIC(m_chk_rcs)))
```

```{r}
# Gamma error distribution with a log link
m_ecm_5y <- gam(BAI_5y ~ ECM + Habitat + Site + log(DBH), family = Gamma(link = "log"), method = "REML", data = guild_bai)
c_ecm_5y <- summary(m_ecm_5y)$p.table["ECM", ]

m_ecm_10y <- gam(BAI_10y ~ ECM + Habitat + Site + log(DBH), family = Gamma(link = "log"), method = "REML", data = guild_bai)
c_ecm_10y <- summary(m_ecm_10y)$p.table["ECM", ]

m_ecm_15y <- gam(BAI_15y ~ ECM + Habitat + Site + log(DBH), family = Gamma(link = "log"), method = "REML", data = guild_bai)
c_ecm_15y <- summary(m_ecm_15y)$p.table["ECM", ]

m_ecm_30y <- gam(BAI_30y ~ ECM + Habitat + Site + log(DBH), family = Gamma(link = "log"), method = "REML", data = guild_bai)
c_ecm_30y <- summary(m_ecm_30y)$p.table["ECM", ]

m_sap_5y <- gam(BAI_5y ~ SAP + Habitat + Site + log(DBH), family = Gamma(link = "log"), method = "REML", data = guild_bai)
c_sap_5y <- summary(m_sap_5y)$p.table["SAP", ]

m_sap_10y <- gam(BAI_10y ~ SAP + Habitat + Site + log(DBH), family = Gamma(link = "log"), method = "REML", data = guild_bai)
c_sap_10y <- summary(m_sap_10y)$p.table["SAP", ]

m_sap_15y <- gam(BAI_15y ~ SAP + Habitat + Site + log(DBH), family = Gamma(link = "log"), method = "REML", data = guild_bai)
c_sap_15y <- summary(m_sap_15y)$p.table["SAP", ]

m_sap_30y <- gam(BAI_30y ~ SAP + Habitat + Site + log(DBH), family = Gamma(link = "log"), method = "REML", data = guild_bai)
c_sap_30y <- summary(m_sap_30y)$p.table["SAP", ]

res <- data.frame(
Guild = c("ECM", "ECM", "ECM", "ECM", "SAP", "SAP", "SAP", "SAP"),
Period = c("BAI_5y", "BAI_10y", "BAI_15y", "BAI_30y", "BAI_5y", "BAI_10y", "BAI_15y", "BAI_30y"),
Estimate = c(c_ecm_5y["Estimate"], c_ecm_10y["Estimate"], c_ecm_15y["Estimate"], c_ecm_30y["Estimate"], c_sap_5y["Estimate"], c_sap_10y["Estimate"], c_sap_15y["Estimate"], c_sap_30y["Estimate"]),
SE = c(c_ecm_5y["Std. Error"], c_ecm_10y["Std. Error"], c_ecm_15y["Std. Error"], c_ecm_30y["Std. Error"], c_sap_5y["Std. Error"], c_sap_10y["Std. Error"], c_sap_15y["Std. Error"], c_sap_30y["Std. Error"]),
p = c(c_ecm_5y["Pr(>|t|)"], c_ecm_10y["Pr(>|t|)"], c_ecm_15y["Pr(>|t|)"], c_ecm_30y["Pr(>|t|)"], c_sap_5y["Pr(>|t|)"], c_sap_10y["Pr(>|t|)"], c_sap_15y["Pr(>|t|)"], c_sap_30y["Pr(>|t|)"]),
R2_adj= c(summary(m_ecm_5y)$r.sq, summary(m_ecm_10y)$r.sq, summary(m_ecm_15y)$r.sq, summary(m_ecm_30y)$r.sq, summary(m_sap_5y)$r.sq, summary(m_sap_10y)$r.sq, summary(m_sap_15y)$r.sq, summary(m_sap_30y)$r.sq),
df= c(df.residual(m_ecm_5y), df.residual(m_ecm_10y), df.residual(m_ecm_15y), df.residual(m_ecm_30y), df.residual(m_sap_5y), df.residual(m_sap_10y), df.residual(m_sap_15y), df.residual(m_sap_30y)))

# BH correction for the 8 tests
res$p_adj <- p.adjust(res$p, method = "BH")   
res$Effect_10pct <- exp(10 * res$Estimate)   # change in BAI per 10 percentage points of guild abundance
res
```

## Guild abundance vs detrended BAI

```{r}
# detrended BAI is already size-free, so no log(DBH) needed
m_ecm_5y <- gam(detBAI_5y ~ ECM + Habitat + Site, family = Gamma(link = "log"), method = "REML", data = guild_bai)
c_ecm_5y <- summary(m_ecm_5y)$p.table["ECM", ]

m_ecm_10y <- gam(detBAI_10y ~ ECM + Habitat + Site, family = Gamma(link = "log"), method = "REML", data = guild_bai)
c_ecm_10y <- summary(m_ecm_10y)$p.table["ECM", ]

m_ecm_15y <- gam(detBAI_15y ~ ECM + Habitat + Site, family = Gamma(link = "log"), method = "REML", data = guild_bai)
c_ecm_15y <- summary(m_ecm_15y)$p.table["ECM", ]

m_ecm_30y <- gam(detBAI_30y ~ ECM + Habitat + Site, family = Gamma(link = "log"), method = "REML", data = guild_bai)
c_ecm_30y <- summary(m_ecm_30y)$p.table["ECM", ]

m_sap_5y <- gam(detBAI_5y ~ SAP + Habitat + Site, family = Gamma(link = "log"), method = "REML", data = guild_bai)
c_sap_5y <- summary(m_sap_5y)$p.table["SAP", ]

m_sap_10y <- gam(detBAI_10y ~ SAP + Habitat + Site, family = Gamma(link = "log"), method = "REML", data = guild_bai)
c_sap_10y <- summary(m_sap_10y)$p.table["SAP", ]

m_sap_15y <- gam(detBAI_15y ~ SAP + Habitat + Site, family = Gamma(link = "log"), method = "REML", data = guild_bai)
c_sap_15y <- summary(m_sap_15y)$p.table["SAP", ]

m_sap_30y <- gam(detBAI_30y ~ SAP + Habitat + Site, family = Gamma(link = "log"), method = "REML", data = guild_bai)
c_sap_30y <- summary(m_sap_30y)$p.table["SAP", ]

res <- data.frame(
Guild = c("ECM", "ECM", "ECM", "ECM", "SAP", "SAP", "SAP", "SAP"),
Period = c("detBAI_5y", "detBAI_10y", "detBAI_15y", "detBAI_30y", "detBAI_5y", "detBAI_10y", "detBAI_15y", "detBAI_30y"),
Estimate = c(c_ecm_5y["Estimate"], c_ecm_10y["Estimate"], c_ecm_15y["Estimate"], c_ecm_30y["Estimate"], c_sap_5y["Estimate"], c_sap_10y["Estimate"], c_sap_15y["Estimate"], c_sap_30y["Estimate"]),
SE = c(c_ecm_5y["Std. Error"], c_ecm_10y["Std. Error"], c_ecm_15y["Std. Error"], c_ecm_30y["Std. Error"], c_sap_5y["Std. Error"], c_sap_10y["Std. Error"], c_sap_15y["Std. Error"], c_sap_30y["Std. Error"]),
p  = c(c_ecm_5y["Pr(>|t|)"], c_ecm_10y["Pr(>|t|)"], c_ecm_15y["Pr(>|t|)"], c_ecm_30y["Pr(>|t|)"], c_sap_5y["Pr(>|t|)"], c_sap_10y["Pr(>|t|)"], c_sap_15y["Pr(>|t|)"], c_sap_30y["Pr(>|t|)"]),
R2_adj   = c(summary(m_ecm_5y)$r.sq, summary(m_ecm_10y)$r.sq, summary(m_ecm_15y)$r.sq, summary(m_ecm_30y)$r.sq, summary(m_sap_5y)$r.sq, summary(m_sap_10y)$r.sq, summary(m_sap_15y)$r.sq, summary(m_sap_30y)$r.sq),
df = c(df.residual(m_ecm_5y), df.residual(m_ecm_10y), df.residual(m_ecm_15y), df.residual(m_ecm_30y), df.residual(m_sap_5y), df.residual(m_sap_10y), df.residual(m_sap_15y), df.residual(m_sap_30y)))

# BH correction for the 8 tests
res$p_adj <- p.adjust(res$p, method = "BH")   
res$Effect_10pct <- exp(10 * res$Estimate)   # change in BAI per 10 percentage points of guild abundance
res
```

## Guild abundance vs RCS BAI

```{r}
# RCS BAI: growth relative to trees of the same cambial age in the same plot (no log(DBH)).
# Each plot has its own regional curve, so the index level differs between plots:
# Plot (= Site x Habitat) is used instead of Habitat + Site

#Check for multicollinearity
vif(lm(log(rcsBAI_10y) ~ ECM + Plot, data = guild_bai))
vif(lm(log(rcsBAI_10y) ~ SAP + Plot, data = guild_bai))

m_ecm_5y <- gam(rcsBAI_5y ~ ECM + Plot, family = Gamma(link = "log"), method = "REML", data = guild_bai)
c_ecm_5y <- summary(m_ecm_5y)$p.table["ECM", ]

m_ecm_10y <- gam(rcsBAI_10y ~ ECM + Plot, family = Gamma(link = "log"), method = "REML", data = guild_bai)
c_ecm_10y <- summary(m_ecm_10y)$p.table["ECM", ]

m_ecm_15y <- gam(rcsBAI_15y ~ ECM + Plot, family = Gamma(link = "log"), method = "REML", data = guild_bai)
c_ecm_15y <- summary(m_ecm_15y)$p.table["ECM", ]

m_ecm_30y <- gam(rcsBAI_30y ~ ECM + Plot, family = Gamma(link = "log"), method = "REML", data = guild_bai)
c_ecm_30y <- summary(m_ecm_30y)$p.table["ECM", ]

m_sap_5y <- gam(rcsBAI_5y ~ SAP + Plot, family = Gamma(link = "log"), method = "REML", data = guild_bai)
c_sap_5y <- summary(m_sap_5y)$p.table["SAP", ]

m_sap_10y <- gam(rcsBAI_10y ~ SAP + Plot, family = Gamma(link = "log"), method = "REML", data = guild_bai)
c_sap_10y <- summary(m_sap_10y)$p.table["SAP", ]

m_sap_15y <- gam(rcsBAI_15y ~ SAP + Plot, family = Gamma(link = "log"), method = "REML", data = guild_bai)
c_sap_15y <- summary(m_sap_15y)$p.table["SAP", ]

m_sap_30y <- gam(rcsBAI_30y ~ SAP + Plot, family = Gamma(link = "log"), method = "REML", data = guild_bai)
c_sap_30y <- summary(m_sap_30y)$p.table["SAP", ]

res <- data.frame(
Guild = c("ECM", "ECM", "ECM", "ECM", "SAP", "SAP", "SAP", "SAP"),
Period = c("rcsBAI_5y", "rcsBAI_10y", "rcsBAI_15y", "rcsBAI_30y", "rcsBAI_5y", "rcsBAI_10y", "rcsBAI_15y", "rcsBAI_30y"),
Estimate = c(c_ecm_5y["Estimate"], c_ecm_10y["Estimate"], c_ecm_15y["Estimate"], c_ecm_30y["Estimate"], c_sap_5y["Estimate"], c_sap_10y["Estimate"], c_sap_15y["Estimate"], c_sap_30y["Estimate"]),
SE  = c(c_ecm_5y["Std. Error"], c_ecm_10y["Std. Error"], c_ecm_15y["Std. Error"], c_ecm_30y["Std. Error"], c_sap_5y["Std. Error"], c_sap_10y["Std. Error"], c_sap_15y["Std. Error"], c_sap_30y["Std. Error"]),
p = c(c_ecm_5y["Pr(>|t|)"], c_ecm_10y["Pr(>|t|)"], c_ecm_15y["Pr(>|t|)"], c_ecm_30y["Pr(>|t|)"], c_sap_5y["Pr(>|t|)"], c_sap_10y["Pr(>|t|)"], c_sap_15y["Pr(>|t|)"], c_sap_30y["Pr(>|t|)"]),
R2_adj = c(summary(m_ecm_5y)$r.sq, summary(m_ecm_10y)$r.sq, summary(m_ecm_15y)$r.sq, summary(m_ecm_30y)$r.sq, summary(m_sap_5y)$r.sq, summary(m_sap_10y)$r.sq, summary(m_sap_15y)$r.sq, summary(m_sap_30y)$r.sq),
df = c(df.residual(m_ecm_5y), df.residual(m_ecm_10y), df.residual(m_ecm_15y), df.residual(m_ecm_30y), df.residual(m_sap_5y), df.residual(m_sap_10y), df.residual(m_sap_15y), df.residual(m_sap_30y)))

# BH correction for the 8 tests
res$p_adj <- p.adjust(res$p, method = "BH")   
res$Effect_10pct <- exp(10 * res$Estimate)   # change in BAI per 10 percentage points of guild abundance
res
```

