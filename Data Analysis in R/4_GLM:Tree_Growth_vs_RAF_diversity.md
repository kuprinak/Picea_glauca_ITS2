---
title: "GAM: Tree growth vs RAF diversity"
author: "Kristina Kuprina"
date: "`r Sys.Date()`"
output: html_document
---

```{r Packages, message=FALSE, warning=FALSE}
library("ggplot2")
library("readr")
library("phyloseq")
library("dplyr")
library("tidyr")
library('tibble')
library("rstatix")
library("vegan")
library("car")
library("mgcv")
```

```{r session info}
sessionInfo()
```

## Read data

```{r sample info, message=FALSE, warning=FALSE}
# output of BAI_detBAI_RCS.Rmd: raw BAI, detrended BAI (detBAI) and RCS BAI (rcsBAI)
sample_info_no_R <-read.table("Alaska_info_no_R.txt", row.names=1, header=T, check.names=F)
tax_tab <- read_delim("Alaska_taxonomy.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
count_tab <-read.table("Alaska_counts.txt", header=T, row.names=1, check.names=F)
count_tab <- count_tab[rownames(count_tab) %in% tax_tab$OTU, ] 
```

## NTC treatment

```{r Treatment of NTC -1}
NTC <-rowSums(count_tab[, 1:8]) # sum of NTC counts for each OTU
sample_counts <- rowSums(count_tab[, 9:143]) # sum of sample counts for each OTU

# normalize sample counts by dividing the samples' total by 17 – as there are 17 as many samples (135) as there are NTCs (8) (135/8=17)
norm_sample_OTU_counts <- sample_counts/17

# which OTUs are deemed likely contaminants based on the threshold noted above:
blank_OTUs <- names(NTC[NTC *10 > norm_sample_OTU_counts]) # OTUs which have more counts than 1/10 of the sample counts
length(blank_OTUs) 

filt_count_tab <- count_tab[!rownames(count_tab) %in% blank_OTUs, -c(1:8)] # removing blank_OTUs and the blank samples from further analysis
filt_count_tab <- filt_count_tab[rowSums(filt_count_tab == 0) != ncol(filt_count_tab), ] # delete OTUs with 0 counts

sample_info_no_R<-sample_info_no_R[order(rownames(sample_info_no_R)), ]
filt_count_tab<- filt_count_tab[,order(names(filt_count_tab))]
filt_count_tab_no_r <- filt_count_tab[colnames(filt_count_tab) %in% rownames(sample_info_no_R) ]
```


# All OTUs

## Data preparation 

```{r}
meta <- as(sample_info_no_R, "data.frame")
keep <- complete.cases(meta[, c("BAI_5y","Habitat","Site")])
meta_clean <- meta[keep, ]

# OTU table + diversity (raw counts)
otu <- otu_table(filt_count_tab_no_r, taxa_are_rows = TRUE)
physeq_tmp <- phyloseq(otu)

div_df <- estimate_richness(physeq_tmp, measures = c("Shannon", "Observed", "InvSimpson"))
rownames(div_df) <- sub("^X", "", rownames(div_df))
rownames(div_df) <- gsub("\\.", "-", rownames(div_df))
div_df$Sample <- rownames(div_df)

meta_clean$Sample <- rownames(meta_clean)

common <- intersect(div_df$Sample, meta_clean$Sample)
div_only <- setdiff(div_df$Sample, meta_clean$Sample)
div_df <- div_df[div_df$Sample %in% common, ]
meta2   <- meta_clean[meta_clean$Sample %in% common, ]
div_df2 <- merge(div_df, meta2, by = "Sample")

div_df2_clean <- div_df2[is.finite(div_df2$BAI_5y), ] # remove samples with NA
```

## Model check: Gaussian vs Gamma (Shannon, BAI 5 years)
 
```{r}
#Check for multicollinearity (should be less than 1.1)
vif(lm(log(BAI_5y) ~ Observed + Habitat + Site + log(DBH), data = div_df2))
vif(lm(log(BAI_5y) ~ Shannon + Habitat + Site + log(DBH), data = div_df2))
vif(lm(log(BAI_5y) ~ InvSimpson + Habitat + Site + log(DBH), data = div_df2))

# same model, two error distributions
m_gaus <- gam(BAI_5y ~ Shannon + Habitat + Site + log(DBH),
  family = gaussian(link = "log"),
  method = "REML",
  data = div_df2_clean)

m_gamma <- gam(BAI_5y ~ Shannon + Habitat + Site + log(DBH),
  family = Gamma(link = "log"),
  method = "REML",
  data = div_df2_clean)

AIC(m_gaus, m_gamma)  # lower = better
gam.check(m_gaus)     # "Resids vs. linear pred." fans out = variance grows with BAI
gam.check(m_gamma)

summary(m_gamma)

# prediction line at the median DBH of each plot, so it shows the Shannon effect only
nd <- expand.grid(Shannon = seq(min(div_df2_clean$Shannon), max(div_df2_clean$Shannon), length.out = 50),
                  Habitat = unique(div_df2_clean$Habitat),
                  Site = unique(div_df2_clean$Site))
med_dbh <- aggregate(DBH ~ Habitat + Site, data = div_df2_clean, FUN = median)
nd <- merge(nd, med_dbh, by = c("Habitat", "Site"))
nd$pred <- predict(m_gamma, newdata = nd, type = "response")

ggplot(div_df2_clean, aes(x = Shannon, y = BAI_5y)) + 
  geom_point(alpha = 0.6, color = "#EE7733") +
  geom_line(data = nd, aes(y = pred), color = "#009988", linewidth = 1) +
  facet_grid(Habitat ~ Site) +
  scale_y_log10() + labs(x = "Shannon index", y = "Mean BAI over 5 years (cm²/year, log scale)") +
  theme_bw()
  
```

## RAF diversity vs raw BAI

```{r}
## Choice of error distribution: Gamma vs Gaussian, both with a log link.
## Positive dAIC = the Gamma fit is better. Checked once here, for one model per
## growth measure, instead of refitting a Gaussian model beside every model below.

m_chk_bai <- gam(BAI_10y ~ Observed + Habitat + Site + log(DBH), family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
g_chk_bai <- gam(BAI_10y ~ Observed + Habitat + Site + log(DBH), family = gaussian(link = "log"), method = "REML", data = div_df2_clean)

m_chk_det <- gam(detBAI_10y ~ Observed + Habitat + Site, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
g_chk_det <- gam(detBAI_10y ~ Observed + Habitat + Site, family = gaussian(link = "log"), method = "REML", data = div_df2_clean)

m_chk_rcs <- gam(rcsBAI_10y ~ Observed + Plot, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
g_chk_rcs <- gam(rcsBAI_10y ~ Observed + Plot, family = gaussian(link = "log"), method = "REML", data = div_df2_clean)

data.frame(Model = c("BAI_10y", "detBAI_10y", "rcsBAI_10y"),
           dAIC  = c(AIC(g_chk_bai) - AIC(m_chk_bai),AIC(g_chk_det) - AIC(m_chk_det), AIC(g_chk_rcs) - AIC(m_chk_rcs)))
```

```{r}
# Gamma error distribution with a log link (see the family check above)
m_obs_5y <- gam(BAI_5y ~ Observed + Habitat + Site + log(DBH), family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_obs_5y <- summary(m_obs_5y)$p.table["Observed", ]

m_obs_10y <- gam(BAI_10y ~ Observed + Habitat + Site + log(DBH), family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_obs_10y <- summary(m_obs_10y)$p.table["Observed", ]

m_obs_15y <- gam(BAI_15y ~ Observed + Habitat + Site + log(DBH), family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_obs_15y <- summary(m_obs_15y)$p.table["Observed", ]

m_obs_30y <- gam(BAI_30y ~ Observed + Habitat + Site + log(DBH), family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_obs_30y <- summary(m_obs_30y)$p.table["Observed", ]

m_sha_5y <- gam(BAI_5y ~ Shannon + Habitat + Site + log(DBH), family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_sha_5y <- summary(m_sha_5y)$p.table["Shannon", ]

m_sha_10y <- gam(BAI_10y ~ Shannon + Habitat + Site + log(DBH), family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_sha_10y <- summary(m_sha_10y)$p.table["Shannon", ]

m_sha_15y <- gam(BAI_15y ~ Shannon + Habitat + Site + log(DBH), family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_sha_15y <- summary(m_sha_15y)$p.table["Shannon", ]

m_sha_30y <- gam(BAI_30y ~ Shannon + Habitat + Site + log(DBH), family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_sha_30y <- summary(m_sha_30y)$p.table["Shannon", ]

m_inv_5y <- gam(BAI_5y ~ InvSimpson + Habitat + Site + log(DBH), family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_inv_5y <- summary(m_inv_5y)$p.table["InvSimpson", ]

m_inv_10y <- gam(BAI_10y ~ InvSimpson + Habitat + Site + log(DBH), family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_inv_10y <- summary(m_inv_10y)$p.table["InvSimpson", ]

m_inv_15y <- gam(BAI_15y ~ InvSimpson + Habitat + Site + log(DBH), family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_inv_15y <- summary(m_inv_15y)$p.table["InvSimpson", ]

m_inv_30y <- gam(BAI_30y ~ InvSimpson + Habitat + Site + log(DBH), family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_inv_30y <- summary(m_inv_30y)$p.table["InvSimpson", ]

res <- data.frame(
Index = c("Observed", "Observed", "Observed", "Observed", "Shannon", "Shannon", "Shannon", "Shannon", "InvSimpson", "InvSimpson", "InvSimpson", "InvSimpson"),
Period = c("BAI_5y", "BAI_10y", "BAI_15y", "BAI_30y", "BAI_5y", "BAI_10y", "BAI_15y", "BAI_30y", "BAI_5y", "BAI_10y", "BAI_15y", "BAI_30y"),
Estimate = c(c_obs_5y["Estimate"], c_obs_10y["Estimate"], c_obs_15y["Estimate"], c_obs_30y["Estimate"], c_sha_5y["Estimate"], c_sha_10y["Estimate"], c_sha_15y["Estimate"], c_sha_30y["Estimate"], c_inv_5y["Estimate"], c_inv_10y["Estimate"], c_inv_15y["Estimate"], c_inv_30y["Estimate"]),
SE  = c(c_obs_5y["Std. Error"], c_obs_10y["Std. Error"], c_obs_15y["Std. Error"], c_obs_30y["Std. Error"], c_sha_5y["Std. Error"], c_sha_10y["Std. Error"], c_sha_15y["Std. Error"], c_sha_30y["Std. Error"], c_inv_5y["Std. Error"], c_inv_10y["Std. Error"], c_inv_15y["Std. Error"], c_inv_30y["Std. Error"]),
p = c(c_obs_5y["Pr(>|t|)"], c_obs_10y["Pr(>|t|)"], c_obs_15y["Pr(>|t|)"], c_obs_30y["Pr(>|t|)"], c_sha_5y["Pr(>|t|)"], c_sha_10y["Pr(>|t|)"], c_sha_15y["Pr(>|t|)"], c_sha_30y["Pr(>|t|)"], c_inv_5y["Pr(>|t|)"], c_inv_10y["Pr(>|t|)"], c_inv_15y["Pr(>|t|)"], c_inv_30y["Pr(>|t|)"]),
R2_adj = c(summary(m_obs_5y)$r.sq, summary(m_obs_10y)$r.sq, summary(m_obs_15y)$r.sq, summary(m_obs_30y)$r.sq, summary(m_sha_5y)$r.sq, summary(m_sha_10y)$r.sq, summary(m_sha_15y)$r.sq, summary(m_sha_30y)$r.sq, summary(m_inv_5y)$r.sq, summary(m_inv_10y)$r.sq, summary(m_inv_15y)$r.sq, summary(m_inv_30y)$r.sq),
df = c(df.residual(m_obs_5y), df.residual(m_obs_10y), df.residual(m_obs_15y), df.residual(m_obs_30y), df.residual(m_sha_5y), df.residual(m_sha_10y), df.residual(m_sha_15y), df.residual(m_sha_30y), df.residual(m_inv_5y), df.residual(m_inv_10y), df.residual(m_inv_15y), df.residual(m_inv_30y)))

# BH correction for the 12 tests
res$p_adj <- p.adjust(res$p, method = "BH")   

## Writing the table of results
sd_obs <- sd(div_df2_clean$Observed)
sd_sha <- sd(div_df2_clean$Shannon)
sd_inv <- sd(div_df2_clean$InvSimpson)

res$CI_low  <- res$Estimate - qt(0.975, res$df) * res$SE
res$CI_high <- res$Estimate + qt(0.975, res$df) * res$SE
res$SD_index   <- c(sd_obs, sd_obs, sd_obs, sd_obs, sd_sha, sd_sha, sd_sha, sd_sha, sd_inv, sd_inv, sd_inv, sd_inv)
res$per_SD_pct <- 100 * (exp(res$Estimate * res$SD_index) - 1)
write.csv(res, "TableST5_All_OTU_BAI.csv", row.names = FALSE)
res
```

## RAF diversity vs detrended BAI

```{r}
# detrended BAI is already size-free, so no log(DBH) term
m_obs_5y <- gam(detBAI_5y ~ Observed + Habitat + Site, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_obs_5y <- summary(m_obs_5y)$p.table["Observed", ]

m_obs_10y <- gam(detBAI_10y ~ Observed + Habitat + Site, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_obs_10y <- summary(m_obs_10y)$p.table["Observed", ]

m_obs_15y <- gam(detBAI_15y ~ Observed + Habitat + Site, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_obs_15y <- summary(m_obs_15y)$p.table["Observed", ]

m_obs_30y <- gam(detBAI_30y ~ Observed + Habitat + Site, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_obs_30y <- summary(m_obs_30y)$p.table["Observed", ]

m_sha_5y <- gam(detBAI_5y ~ Shannon + Habitat + Site, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_sha_5y <- summary(m_sha_5y)$p.table["Shannon", ]

m_sha_10y <- gam(detBAI_10y ~ Shannon + Habitat + Site, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_sha_10y <- summary(m_sha_10y)$p.table["Shannon", ]

m_sha_15y <- gam(detBAI_15y ~ Shannon + Habitat + Site, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_sha_15y <- summary(m_sha_15y)$p.table["Shannon", ]

m_sha_30y <- gam(detBAI_30y ~ Shannon + Habitat + Site, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_sha_30y <- summary(m_sha_30y)$p.table["Shannon", ]

m_inv_5y <- gam(detBAI_5y ~ InvSimpson + Habitat + Site, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_inv_5y <- summary(m_inv_5y)$p.table["InvSimpson", ]

m_inv_10y <- gam(detBAI_10y ~ InvSimpson + Habitat + Site, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_inv_10y <- summary(m_inv_10y)$p.table["InvSimpson", ]

m_inv_15y <- gam(detBAI_15y ~ InvSimpson + Habitat + Site, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_inv_15y <- summary(m_inv_15y)$p.table["InvSimpson", ]

m_inv_30y <- gam(detBAI_30y ~ InvSimpson + Habitat + Site, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_inv_30y <- summary(m_inv_30y)$p.table["InvSimpson", ]

res <- data.frame(
Index = c("Observed", "Observed", "Observed", "Observed", "Shannon", "Shannon", "Shannon", "Shannon", "InvSimpson", "InvSimpson", "InvSimpson", "InvSimpson"),
Period = c("detBAI_5y", "detBAI_10y", "detBAI_15y", "detBAI_30y", "detBAI_5y", "detBAI_10y", "detBAI_15y", "detBAI_30y", "detBAI_5y", "detBAI_10y", "detBAI_15y", "detBAI_30y"),
Estimate = c(c_obs_5y["Estimate"], c_obs_10y["Estimate"], c_obs_15y["Estimate"], c_obs_30y["Estimate"], c_sha_5y["Estimate"], c_sha_10y["Estimate"], c_sha_15y["Estimate"], c_sha_30y["Estimate"], c_inv_5y["Estimate"], c_inv_10y["Estimate"], c_inv_15y["Estimate"], c_inv_30y["Estimate"]),
SE = c(c_obs_5y["Std. Error"], c_obs_10y["Std. Error"], c_obs_15y["Std. Error"], c_obs_30y["Std. Error"], c_sha_5y["Std. Error"], c_sha_10y["Std. Error"], c_sha_15y["Std. Error"], c_sha_30y["Std. Error"], c_inv_5y["Std. Error"], c_inv_10y["Std. Error"], c_inv_15y["Std. Error"], c_inv_30y["Std. Error"]),
p = c(c_obs_5y["Pr(>|t|)"], c_obs_10y["Pr(>|t|)"], c_obs_15y["Pr(>|t|)"], c_obs_30y["Pr(>|t|)"], c_sha_5y["Pr(>|t|)"], c_sha_10y["Pr(>|t|)"], c_sha_15y["Pr(>|t|)"], c_sha_30y["Pr(>|t|)"], c_inv_5y["Pr(>|t|)"], c_inv_10y["Pr(>|t|)"], c_inv_15y["Pr(>|t|)"], c_inv_30y["Pr(>|t|)"]),
R2_adj = c(summary(m_obs_5y)$r.sq, summary(m_obs_10y)$r.sq, summary(m_obs_15y)$r.sq, summary(m_obs_30y)$r.sq, summary(m_sha_5y)$r.sq, summary(m_sha_10y)$r.sq, summary(m_sha_15y)$r.sq, summary(m_sha_30y)$r.sq, summary(m_inv_5y)$r.sq, summary(m_inv_10y)$r.sq, summary(m_inv_15y)$r.sq, summary(m_inv_30y)$r.sq),
df = c(df.residual(m_obs_5y), df.residual(m_obs_10y), df.residual(m_obs_15y), df.residual(m_obs_30y), df.residual(m_sha_5y), df.residual(m_sha_10y), df.residual(m_sha_15y), df.residual(m_sha_30y), df.residual(m_inv_5y), df.residual(m_inv_10y), df.residual(m_inv_15y), df.residual(m_inv_30y)))

# BH correction for the 12 tests
res$p_adj <- p.adjust(res$p, method = "BH")   

## Writing the table of results
sd_obs <- sd(div_df2_clean$Observed)
sd_sha <- sd(div_df2_clean$Shannon)
sd_inv <- sd(div_df2_clean$InvSimpson)

res$CI_low  <- res$Estimate - qt(0.975, res$df) * res$SE
res$CI_high <- res$Estimate + qt(0.975, res$df) * res$SE
res$SD_index   <- c(sd_obs, sd_obs, sd_obs, sd_obs, sd_sha, sd_sha, sd_sha, sd_sha, sd_inv, sd_inv, sd_inv, sd_inv)
res$per_SD_pct <- 100 * (exp(res$Estimate * res$SD_index) - 1)
write.csv(res, "TableST5_All_OTU_detBAI.csv", row.names = FALSE)
res
```

## RAF diversity vs RCS BAI

```{r}
# RCS BAI: growth relative to trees of the same cambial age in the same plot (no log(DBH) term).
# Each plot has its own regional curve, so the index level differs between plots:
# Plot (= Site x Habitat) is used instead of Habitat + Site

#Check for multicollinearity
vif(lm(log(rcsBAI_5y) ~ Observed + Plot, data = div_df2))
vif(lm(log(rcsBAI_5y) ~ Shannon + Plot, data = div_df2))
vif(lm(log(rcsBAI_5y) ~ InvSimpson + Plot, data = div_df2))

m_obs_5y <- gam(rcsBAI_5y ~ Observed + Plot, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_obs_5y <- summary(m_obs_5y)$p.table["Observed", ]

m_obs_10y <- gam(rcsBAI_10y ~ Observed + Plot, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_obs_10y <- summary(m_obs_10y)$p.table["Observed", ]

m_obs_15y <- gam(rcsBAI_15y ~ Observed + Plot, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_obs_15y <- summary(m_obs_15y)$p.table["Observed", ]

m_obs_30y <- gam(rcsBAI_30y ~ Observed + Plot, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_obs_30y <- summary(m_obs_30y)$p.table["Observed", ]

m_sha_5y <- gam(rcsBAI_5y ~ Shannon + Plot, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_sha_5y <- summary(m_sha_5y)$p.table["Shannon", ]

m_sha_10y <- gam(rcsBAI_10y ~ Shannon + Plot, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_sha_10y <- summary(m_sha_10y)$p.table["Shannon", ]

m_sha_15y <- gam(rcsBAI_15y ~ Shannon + Plot, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_sha_15y <- summary(m_sha_15y)$p.table["Shannon", ]

m_sha_30y <- gam(rcsBAI_30y ~ Shannon + Plot, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_sha_30y <- summary(m_sha_30y)$p.table["Shannon", ]

m_inv_5y <- gam(rcsBAI_5y ~ InvSimpson + Plot, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_inv_5y <- summary(m_inv_5y)$p.table["InvSimpson", ]

m_inv_10y <- gam(rcsBAI_10y ~ InvSimpson + Plot, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_inv_10y <- summary(m_inv_10y)$p.table["InvSimpson", ]

m_inv_15y <- gam(rcsBAI_15y ~ InvSimpson + Plot, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_inv_15y <- summary(m_inv_15y)$p.table["InvSimpson", ]

m_inv_30y <- gam(rcsBAI_30y ~ InvSimpson + Plot, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_inv_30y <- summary(m_inv_30y)$p.table["InvSimpson", ]

res <- data.frame(
Index = c("Observed", "Observed", "Observed", "Observed", "Shannon", "Shannon", "Shannon", "Shannon", "InvSimpson", "InvSimpson", "InvSimpson", "InvSimpson"),
Period = c("rcsBAI_5y", "rcsBAI_10y", "rcsBAI_15y", "rcsBAI_30y", "rcsBAI_5y", "rcsBAI_10y", "rcsBAI_15y", "rcsBAI_30y", "rcsBAI_5y", "rcsBAI_10y", "rcsBAI_15y", "rcsBAI_30y"),
Estimate = c(c_obs_5y["Estimate"], c_obs_10y["Estimate"], c_obs_15y["Estimate"], c_obs_30y["Estimate"], c_sha_5y["Estimate"], c_sha_10y["Estimate"], c_sha_15y["Estimate"], c_sha_30y["Estimate"], c_inv_5y["Estimate"], c_inv_10y["Estimate"], c_inv_15y["Estimate"], c_inv_30y["Estimate"]),
SE = c(c_obs_5y["Std. Error"], c_obs_10y["Std. Error"], c_obs_15y["Std. Error"], c_obs_30y["Std. Error"], c_sha_5y["Std. Error"], c_sha_10y["Std. Error"], c_sha_15y["Std. Error"], c_sha_30y["Std. Error"], c_inv_5y["Std. Error"], c_inv_10y["Std. Error"], c_inv_15y["Std. Error"], c_inv_30y["Std. Error"]),
p = c(c_obs_5y["Pr(>|t|)"], c_obs_10y["Pr(>|t|)"], c_obs_15y["Pr(>|t|)"], c_obs_30y["Pr(>|t|)"], c_sha_5y["Pr(>|t|)"], c_sha_10y["Pr(>|t|)"], c_sha_15y["Pr(>|t|)"], c_sha_30y["Pr(>|t|)"], c_inv_5y["Pr(>|t|)"], c_inv_10y["Pr(>|t|)"], c_inv_15y["Pr(>|t|)"], c_inv_30y["Pr(>|t|)"]),
R2_adj = c(summary(m_obs_5y)$r.sq, summary(m_obs_10y)$r.sq, summary(m_obs_15y)$r.sq, summary(m_obs_30y)$r.sq, summary(m_sha_5y)$r.sq, summary(m_sha_10y)$r.sq, summary(m_sha_15y)$r.sq, summary(m_sha_30y)$r.sq, summary(m_inv_5y)$r.sq, summary(m_inv_10y)$r.sq, summary(m_inv_15y)$r.sq, summary(m_inv_30y)$r.sq),
df = c(df.residual(m_obs_5y), df.residual(m_obs_10y), df.residual(m_obs_15y), df.residual(m_obs_30y), df.residual(m_sha_5y), df.residual(m_sha_10y), df.residual(m_sha_15y), df.residual(m_sha_30y), df.residual(m_inv_5y), df.residual(m_inv_10y), df.residual(m_inv_15y), df.residual(m_inv_30y)))

# BH correction for the 12 tests
res$p_adj <- p.adjust(res$p, method = "BH")   

## Writing the table of results
sd_obs <- sd(div_df2_clean$Observed)
sd_sha <- sd(div_df2_clean$Shannon)
sd_inv <- sd(div_df2_clean$InvSimpson)

res$CI_low  <- res$Estimate - qt(0.975, res$df) * res$SE
res$CI_high <- res$Estimate + qt(0.975, res$df) * res$SE
res$SD_index   <- c(sd_obs, sd_obs, sd_obs, sd_obs, sd_sha, sd_sha, sd_sha, sd_sha, sd_inv, sd_inv, sd_inv, sd_inv)
res$per_SD_pct <- 100 * (exp(res$Estimate * res$SD_index) - 1)
write.csv(res, "TableST5_All_OTU_rcsBAI.csv", row.names = FALSE)
res
```

## RAF diversity vs DBH and age

```{r}
# DBH and age as response: tests whether tree size/age differ with RAF diversity (association only)
m_obs_dbh <- gam(DBH ~ Observed + Habitat + Site, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_obs_dbh <- summary(m_obs_dbh)$p.table["Observed", ]

m_obs_age <- gam(Age ~ Observed + Habitat + Site, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_obs_age <- summary(m_obs_age)$p.table["Observed", ]

m_sha_dbh <- gam(DBH ~ Shannon + Habitat + Site, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_sha_dbh <- summary(m_sha_dbh)$p.table["Shannon", ]

m_sha_age <- gam(Age ~ Shannon + Habitat + Site, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_sha_age <- summary(m_sha_age)$p.table["Shannon", ]

m_inv_dbh <- gam(DBH ~ InvSimpson + Habitat + Site, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_inv_dbh <- summary(m_inv_dbh)$p.table["InvSimpson", ]

m_inv_age <- gam(Age ~ InvSimpson + Habitat + Site, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_inv_age <- summary(m_inv_age)$p.table["InvSimpson", ]

res <- data.frame(
Index = c("Observed", "Observed", "Shannon", "Shannon", "InvSimpson", "InvSimpson"),
Period = c("DBH", "Age", "DBH", "Age", "DBH", "Age"),
Estimate = c(c_obs_dbh["Estimate"], c_obs_age["Estimate"], c_sha_dbh["Estimate"], c_sha_age["Estimate"], c_inv_dbh["Estimate"], c_inv_age["Estimate"]),
SE = c(c_obs_dbh["Std. Error"], c_obs_age["Std. Error"], c_sha_dbh["Std. Error"], c_sha_age["Std. Error"], c_inv_dbh["Std. Error"], c_inv_age["Std. Error"]),
p  = c(c_obs_dbh["Pr(>|t|)"], c_obs_age["Pr(>|t|)"], c_sha_dbh["Pr(>|t|)"], c_sha_age["Pr(>|t|)"], c_inv_dbh["Pr(>|t|)"], c_inv_age["Pr(>|t|)"]),
R2_adj = c(summary(m_obs_dbh)$r.sq, summary(m_obs_age)$r.sq, summary(m_sha_dbh)$r.sq, summary(m_sha_age)$r.sq, summary(m_inv_dbh)$r.sq, summary(m_inv_age)$r.sq),
df = c(df.residual(m_obs_dbh), df.residual(m_obs_age), df.residual(m_sha_dbh), df.residual(m_sha_age), df.residual(m_inv_dbh), df.residual(m_inv_age)))

# BH correction for the 6 tests
res$p_adj <- p.adjust(res$p, method = "BH")  
res

```

# ECM OTUs

## Data preparation (ECM)

```{r}
# OTU table + diversity
# Raw count phyloseq object - select ECM
otu <- otu_table(filt_count_tab_no_r, taxa_are_rows = TRUE)
otu_ECM <- tax_tab[tax_tab$Guild %in% c("Ectomycorrhizal"),]
ECM_OTUs_present <- intersect(otu_ECM$OTU, rownames(filt_count_tab_no_r))

# create phyloseq object with raw ECM counts
otu_phy_ECM <- otu_table(filt_count_tab_no_r[ECM_OTUs_present,], taxa_are_rows=T)
otu_phy_ECM_data <- sample_data(sample_info_no_R)
phy_ECM <- phyloseq(otu_phy_ECM, otu_phy_ECM_data)

div_df <- estimate_richness(phy_ECM, measures = c("Shannon", "Observed", "InvSimpson"))
rownames(div_df) <- sub("^X", "", rownames(div_df))
rownames(div_df) <- gsub("\\.", "-", rownames(div_df))
div_df$Sample <- rownames(div_df)

meta_clean$Sample <- rownames(meta_clean)

common <- intersect(div_df$Sample, meta_clean$Sample)
div_only <- setdiff(div_df$Sample, meta_clean$Sample)
div_df <- div_df[div_df$Sample %in% common, ]
meta2   <- meta_clean[meta_clean$Sample %in% common, ]
div_df2 <- merge(div_df, meta2, by = "Sample")

div_df2_clean <- div_df2[is.finite(div_df2$BAI_5y), ] # remove samples with NA
```

## RAF diversity vs raw BAI (ECM)

```{r}
#Check for multicollinearity
vif(lm(log(BAI_5y) ~ Observed + Habitat + Site + log(DBH), data = div_df2))
vif(lm(log(BAI_5y) ~ Shannon + Habitat + Site + log(DBH), data = div_df2))
vif(lm(log(BAI_5y) ~ InvSimpson + Habitat + Site + log(DBH), data = div_df2))

m_obs_5y <- gam(BAI_5y ~ Observed + Habitat + Site + log(DBH), family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_obs_5y <- summary(m_obs_5y)$p.table["Observed", ]

m_obs_10y <- gam(BAI_10y ~ Observed + Habitat + Site + log(DBH), family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_obs_10y <- summary(m_obs_10y)$p.table["Observed", ]

m_obs_15y <- gam(BAI_15y ~ Observed + Habitat + Site + log(DBH), family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_obs_15y <- summary(m_obs_15y)$p.table["Observed", ]

m_obs_30y <- gam(BAI_30y ~ Observed + Habitat + Site + log(DBH), family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_obs_30y <- summary(m_obs_30y)$p.table["Observed", ]

m_sha_5y <- gam(BAI_5y ~ Shannon + Habitat + Site + log(DBH), family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_sha_5y <- summary(m_sha_5y)$p.table["Shannon", ]

m_sha_10y <- gam(BAI_10y ~ Shannon + Habitat + Site + log(DBH), family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_sha_10y <- summary(m_sha_10y)$p.table["Shannon", ]

m_sha_15y <- gam(BAI_15y ~ Shannon + Habitat + Site + log(DBH), family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_sha_15y <- summary(m_sha_15y)$p.table["Shannon", ]

m_sha_30y <- gam(BAI_30y ~ Shannon + Habitat + Site + log(DBH), family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_sha_30y <- summary(m_sha_30y)$p.table["Shannon", ]

m_inv_5y <- gam(BAI_5y ~ InvSimpson + Habitat + Site + log(DBH), family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_inv_5y <- summary(m_inv_5y)$p.table["InvSimpson", ]

m_inv_10y <- gam(BAI_10y ~ InvSimpson + Habitat + Site + log(DBH), family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_inv_10y <- summary(m_inv_10y)$p.table["InvSimpson", ]

m_inv_15y <- gam(BAI_15y ~ InvSimpson + Habitat + Site + log(DBH), family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_inv_15y <- summary(m_inv_15y)$p.table["InvSimpson", ]

m_inv_30y <- gam(BAI_30y ~ InvSimpson + Habitat + Site + log(DBH), family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_inv_30y <- summary(m_inv_30y)$p.table["InvSimpson", ]

res <- data.frame(
Index = c("Observed", "Observed", "Observed", "Observed", "Shannon", "Shannon", "Shannon", "Shannon", "InvSimpson", "InvSimpson", "InvSimpson", "InvSimpson"),
Period = c("BAI_5y", "BAI_10y", "BAI_15y", "BAI_30y", "BAI_5y", "BAI_10y", "BAI_15y", "BAI_30y", "BAI_5y", "BAI_10y", "BAI_15y", "BAI_30y"),
Estimate = c(c_obs_5y["Estimate"], c_obs_10y["Estimate"], c_obs_15y["Estimate"], c_obs_30y["Estimate"], c_sha_5y["Estimate"], c_sha_10y["Estimate"], c_sha_15y["Estimate"], c_sha_30y["Estimate"], c_inv_5y["Estimate"], c_inv_10y["Estimate"], c_inv_15y["Estimate"], c_inv_30y["Estimate"]),
SE = c(c_obs_5y["Std. Error"], c_obs_10y["Std. Error"], c_obs_15y["Std. Error"], c_obs_30y["Std. Error"], c_sha_5y["Std. Error"], c_sha_10y["Std. Error"], c_sha_15y["Std. Error"], c_sha_30y["Std. Error"], c_inv_5y["Std. Error"], c_inv_10y["Std. Error"], c_inv_15y["Std. Error"], c_inv_30y["Std. Error"]),
p = c(c_obs_5y["Pr(>|t|)"], c_obs_10y["Pr(>|t|)"], c_obs_15y["Pr(>|t|)"], c_obs_30y["Pr(>|t|)"], c_sha_5y["Pr(>|t|)"], c_sha_10y["Pr(>|t|)"], c_sha_15y["Pr(>|t|)"], c_sha_30y["Pr(>|t|)"], c_inv_5y["Pr(>|t|)"], c_inv_10y["Pr(>|t|)"], c_inv_15y["Pr(>|t|)"], c_inv_30y["Pr(>|t|)"]),
R2_adj = c(summary(m_obs_5y)$r.sq, summary(m_obs_10y)$r.sq, summary(m_obs_15y)$r.sq, summary(m_obs_30y)$r.sq, summary(m_sha_5y)$r.sq, summary(m_sha_10y)$r.sq, summary(m_sha_15y)$r.sq, summary(m_sha_30y)$r.sq, summary(m_inv_5y)$r.sq, summary(m_inv_10y)$r.sq, summary(m_inv_15y)$r.sq, summary(m_inv_30y)$r.sq),
df = c(df.residual(m_obs_5y), df.residual(m_obs_10y), df.residual(m_obs_15y), df.residual(m_obs_30y), df.residual(m_sha_5y), df.residual(m_sha_10y), df.residual(m_sha_15y), df.residual(m_sha_30y), df.residual(m_inv_5y), df.residual(m_inv_10y), df.residual(m_inv_15y), df.residual(m_inv_30y)))

# BH correction for the 12 tests
res$p_adj <- p.adjust(res$p, method = "BH")  

## Writing the table of results
sd_obs <- sd(div_df2_clean$Observed)
sd_sha <- sd(div_df2_clean$Shannon)
sd_inv <- sd(div_df2_clean$InvSimpson)

res$CI_low  <- res$Estimate - qt(0.975, res$df) * res$SE
res$CI_high <- res$Estimate + qt(0.975, res$df) * res$SE
res$SD_index   <- c(sd_obs, sd_obs, sd_obs, sd_obs, sd_sha, sd_sha, sd_sha, sd_sha, sd_inv, sd_inv, sd_inv, sd_inv)
res$per_SD_pct <- 100 * (exp(res$Estimate * res$SD_index) - 1)
write.csv(res, "TableST5_ECM_OTU_BAI.csv", row.names = FALSE)
res
```

## RAF diversity vs detrended BAI (ECM)

```{r}
m_obs_5y <- gam(detBAI_5y ~ Observed + Habitat + Site, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_obs_5y <- summary(m_obs_5y)$p.table["Observed", ]

m_obs_10y <- gam(detBAI_10y ~ Observed + Habitat + Site, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_obs_10y <- summary(m_obs_10y)$p.table["Observed", ]

m_obs_15y <- gam(detBAI_15y ~ Observed + Habitat + Site, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_obs_15y <- summary(m_obs_15y)$p.table["Observed", ]

m_obs_30y <- gam(detBAI_30y ~ Observed + Habitat + Site, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_obs_30y <- summary(m_obs_30y)$p.table["Observed", ]

m_sha_5y <- gam(detBAI_5y ~ Shannon + Habitat + Site, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_sha_5y <- summary(m_sha_5y)$p.table["Shannon", ]

m_sha_10y <- gam(detBAI_10y ~ Shannon + Habitat + Site, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_sha_10y <- summary(m_sha_10y)$p.table["Shannon", ]

m_sha_15y <- gam(detBAI_15y ~ Shannon + Habitat + Site, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_sha_15y <- summary(m_sha_15y)$p.table["Shannon", ]

m_sha_30y <- gam(detBAI_30y ~ Shannon + Habitat + Site, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_sha_30y <- summary(m_sha_30y)$p.table["Shannon", ]

m_inv_5y <- gam(detBAI_5y ~ InvSimpson + Habitat + Site, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_inv_5y <- summary(m_inv_5y)$p.table["InvSimpson", ]

m_inv_10y <- gam(detBAI_10y ~ InvSimpson + Habitat + Site, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_inv_10y <- summary(m_inv_10y)$p.table["InvSimpson", ]

m_inv_15y <- gam(detBAI_15y ~ InvSimpson + Habitat + Site, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_inv_15y <- summary(m_inv_15y)$p.table["InvSimpson", ]

m_inv_30y <- gam(detBAI_30y ~ InvSimpson + Habitat + Site, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_inv_30y <- summary(m_inv_30y)$p.table["InvSimpson", ]

res <- data.frame(
Index = c("Observed", "Observed", "Observed", "Observed", "Shannon", "Shannon", "Shannon", "Shannon", "InvSimpson", "InvSimpson", "InvSimpson", "InvSimpson"),
Period = c("detBAI_5y", "detBAI_10y", "detBAI_15y", "detBAI_30y", "detBAI_5y", "detBAI_10y", "detBAI_15y", "detBAI_30y", "detBAI_5y", "detBAI_10y", "detBAI_15y", "detBAI_30y"),
Estimate = c(c_obs_5y["Estimate"], c_obs_10y["Estimate"], c_obs_15y["Estimate"], c_obs_30y["Estimate"], c_sha_5y["Estimate"], c_sha_10y["Estimate"], c_sha_15y["Estimate"], c_sha_30y["Estimate"], c_inv_5y["Estimate"], c_inv_10y["Estimate"], c_inv_15y["Estimate"], c_inv_30y["Estimate"]),
SE = c(c_obs_5y["Std. Error"], c_obs_10y["Std. Error"], c_obs_15y["Std. Error"], c_obs_30y["Std. Error"], c_sha_5y["Std. Error"], c_sha_10y["Std. Error"], c_sha_15y["Std. Error"], c_sha_30y["Std. Error"], c_inv_5y["Std. Error"], c_inv_10y["Std. Error"], c_inv_15y["Std. Error"], c_inv_30y["Std. Error"]),
p  = c(c_obs_5y["Pr(>|t|)"], c_obs_10y["Pr(>|t|)"], c_obs_15y["Pr(>|t|)"], c_obs_30y["Pr(>|t|)"], c_sha_5y["Pr(>|t|)"], c_sha_10y["Pr(>|t|)"], c_sha_15y["Pr(>|t|)"], c_sha_30y["Pr(>|t|)"], c_inv_5y["Pr(>|t|)"], c_inv_10y["Pr(>|t|)"], c_inv_15y["Pr(>|t|)"], c_inv_30y["Pr(>|t|)"]),
R2_adj  = c(summary(m_obs_5y)$r.sq, summary(m_obs_10y)$r.sq, summary(m_obs_15y)$r.sq, summary(m_obs_30y)$r.sq, summary(m_sha_5y)$r.sq, summary(m_sha_10y)$r.sq, summary(m_sha_15y)$r.sq, summary(m_sha_30y)$r.sq, summary(m_inv_5y)$r.sq, summary(m_inv_10y)$r.sq, summary(m_inv_15y)$r.sq, summary(m_inv_30y)$r.sq),
df = c(df.residual(m_obs_5y), df.residual(m_obs_10y), df.residual(m_obs_15y), df.residual(m_obs_30y), df.residual(m_sha_5y), df.residual(m_sha_10y), df.residual(m_sha_15y), df.residual(m_sha_30y), df.residual(m_inv_5y), df.residual(m_inv_10y), df.residual(m_inv_15y), df.residual(m_inv_30y)))

res$p_adj <- p.adjust(res$p, method = "BH")   # BH correction for the 12 tests

## Writing the table of results
sd_obs <- sd(div_df2_clean$Observed)
sd_sha <- sd(div_df2_clean$Shannon)
sd_inv <- sd(div_df2_clean$InvSimpson)

res$CI_low  <- res$Estimate - qt(0.975, res$df) * res$SE
res$CI_high <- res$Estimate + qt(0.975, res$df) * res$SE
res$SD_index   <- c(sd_obs, sd_obs, sd_obs, sd_obs, sd_sha, sd_sha, sd_sha, sd_sha, sd_inv, sd_inv, sd_inv, sd_inv)
res$per_SD_pct <- 100 * (exp(res$Estimate * res$SD_index) - 1)
write.csv(res, "TableST2_ECM_OTU_detBAI.csv", row.names = FALSE)
res
```

## RAF diversity vs rcsBAI (ECM)

```{r}
# RCS BAI: growth relative to trees of the same cambial age in the same plot (no log(DBH) term).
# Each plot has its own regional curve, so the index level differs between plots:
# Plot (= Site x Habitat) is used instead of Habitat + Site

#Check for multicollinearity
vif(lm(log(rcsBAI_5y) ~ Observed + Plot, data = div_df2))
vif(lm(log(rcsBAI_5y) ~ Shannon + Plot, data = div_df2))
vif(lm(log(rcsBAI_5y) ~ InvSimpson + Plot, data = div_df2))

m_obs_5y <- gam(rcsBAI_5y ~ Observed + Plot, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_obs_5y <- summary(m_obs_5y)$p.table["Observed", ]

m_obs_10y <- gam(rcsBAI_10y ~ Observed + Plot, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_obs_10y <- summary(m_obs_10y)$p.table["Observed", ]

m_obs_15y <- gam(rcsBAI_15y ~ Observed + Plot, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_obs_15y <- summary(m_obs_15y)$p.table["Observed", ]

m_obs_30y <- gam(rcsBAI_30y ~ Observed + Plot, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_obs_30y <- summary(m_obs_30y)$p.table["Observed", ]

m_sha_5y <- gam(rcsBAI_5y ~ Shannon + Plot, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_sha_5y <- summary(m_sha_5y)$p.table["Shannon", ]

m_sha_10y <- gam(rcsBAI_10y ~ Shannon + Plot, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_sha_10y <- summary(m_sha_10y)$p.table["Shannon", ]

m_sha_15y <- gam(rcsBAI_15y ~ Shannon + Plot, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_sha_15y <- summary(m_sha_15y)$p.table["Shannon", ]

m_sha_30y <- gam(rcsBAI_30y ~ Shannon + Plot, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_sha_30y <- summary(m_sha_30y)$p.table["Shannon", ]

m_inv_5y <- gam(rcsBAI_5y ~ InvSimpson + Plot, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_inv_5y <- summary(m_inv_5y)$p.table["InvSimpson", ]

m_inv_10y <- gam(rcsBAI_10y ~ InvSimpson + Plot, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_inv_10y <- summary(m_inv_10y)$p.table["InvSimpson", ]

m_inv_15y <- gam(rcsBAI_15y ~ InvSimpson + Plot, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_inv_15y <- summary(m_inv_15y)$p.table["InvSimpson", ]

m_inv_30y <- gam(rcsBAI_30y ~ InvSimpson + Plot, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
c_inv_30y <- summary(m_inv_30y)$p.table["InvSimpson", ]

res <- data.frame(
Index = c("Observed", "Observed", "Observed", "Observed", "Shannon", "Shannon", "Shannon", "Shannon", "InvSimpson", "InvSimpson", "InvSimpson", "InvSimpson"),
Period = c("rcsBAI_5y", "rcsBAI_10y", "rcsBAI_15y", "rcsBAI_30y", "rcsBAI_5y", "rcsBAI_10y", "rcsBAI_15y", "rcsBAI_30y", "rcsBAI_5y", "rcsBAI_10y", "rcsBAI_15y", "rcsBAI_30y"),
Estimate = c(c_obs_5y["Estimate"], c_obs_10y["Estimate"], c_obs_15y["Estimate"], c_obs_30y["Estimate"], c_sha_5y["Estimate"], c_sha_10y["Estimate"], c_sha_15y["Estimate"], c_sha_30y["Estimate"], c_inv_5y["Estimate"], c_inv_10y["Estimate"], c_inv_15y["Estimate"], c_inv_30y["Estimate"]),
SE = c(c_obs_5y["Std. Error"], c_obs_10y["Std. Error"], c_obs_15y["Std. Error"], c_obs_30y["Std. Error"], c_sha_5y["Std. Error"], c_sha_10y["Std. Error"], c_sha_15y["Std. Error"], c_sha_30y["Std. Error"], c_inv_5y["Std. Error"], c_inv_10y["Std. Error"], c_inv_15y["Std. Error"], c_inv_30y["Std. Error"]),
p = c(c_obs_5y["Pr(>|t|)"], c_obs_10y["Pr(>|t|)"], c_obs_15y["Pr(>|t|)"], c_obs_30y["Pr(>|t|)"], c_sha_5y["Pr(>|t|)"], c_sha_10y["Pr(>|t|)"], c_sha_15y["Pr(>|t|)"], c_sha_30y["Pr(>|t|)"], c_inv_5y["Pr(>|t|)"], c_inv_10y["Pr(>|t|)"], c_inv_15y["Pr(>|t|)"], c_inv_30y["Pr(>|t|)"]),
R2_adj= c(summary(m_obs_5y)$r.sq, summary(m_obs_10y)$r.sq, summary(m_obs_15y)$r.sq, summary(m_obs_30y)$r.sq, summary(m_sha_5y)$r.sq, summary(m_sha_10y)$r.sq, summary(m_sha_15y)$r.sq, summary(m_sha_30y)$r.sq, summary(m_inv_5y)$r.sq, summary(m_inv_10y)$r.sq, summary(m_inv_15y)$r.sq, summary(m_inv_30y)$r.sq),
df = c(df.residual(m_obs_5y), df.residual(m_obs_10y), df.residual(m_obs_15y), df.residual(m_obs_30y), df.residual(m_sha_5y), df.residual(m_sha_10y), df.residual(m_sha_15y), df.residual(m_sha_30y), df.residual(m_inv_5y), df.residual(m_inv_10y), df.residual(m_inv_15y), df.residual(m_inv_30y)))

# BH correction for the 12 tests
res$p_adj <- p.adjust(res$p, method = "BH")  

## Writing the table of results
sd_obs <- sd(div_df2_clean$Observed)
sd_sha <- sd(div_df2_clean$Shannon)
sd_inv <- sd(div_df2_clean$InvSimpson)

res$CI_low  <- res$Estimate - qt(0.975, res$df) * res$SE
res$CI_high <- res$Estimate + qt(0.975, res$df) * res$SE
res$SD_index   <- c(sd_obs, sd_obs, sd_obs, sd_obs, sd_sha, sd_sha, sd_sha, sd_sha, sd_inv, sd_inv, sd_inv, sd_inv)
res$per_SD_pct <- 100 * (exp(res$Estimate * res$SD_index) - 1)
write.csv(res, "TableST5_ECM_OTU_rcsBAI.csv", row.names = FALSE)
res
```


# ECM diversity vs tree age

## Distributions and transformations

```{r}

par(mfrow = c(2, 3))
hist(div_df2_clean$Age, main = "Age")
hist(log(div_df2_clean$Age), main = "log(Age)")
hist(div_df2_clean$Observed, main = "OTU richness")   # counts, right-skewed
hist(div_df2_clean$Shannon, main = "Shannon")         # close to symmetric
hist(div_df2_clean$InvSimpson, main = "InvSimpson")   # positive, right-skewed
hist(log(div_df2_clean$InvSimpson), main = "log(InvSimpson)")
par(mfrow = c(1, 1))

# tree age is right-skewed (skewness), log(Age) is more symmetric
```

## OTU richness vs age

```{r}
# counts and overdispersed: negative binomial, or Gaussian on log(richness)
m_nb <- gam(Observed ~ log(Age) + Habitat + Site, family = nb(), method = "REML", data = div_df2_clean)
m_ln <- gam(log(Observed) ~ log(Age) + Habitat + Site, family = gaussian(), method = "REML", data = div_df2_clean)

# the two AICs are on different response scales; the Jacobian term makes them comparable
AIC(m_nb)
AIC(m_ln) + 2 * sum(log(div_df2_clean$Observed))

# is the age effect linear? edf close to 1 = yes, a linear term is enough (we have 1.001)
m_s <- gam(Observed ~ s(log(Age), k = 5) + Habitat + Site, family = nb(), method = "REML", data = div_df2_clean)
summary(m_s)

summary(m_ln)
gam.check(m_ln)
```

## Shannon index vs age

```{r}
m_sh <- gam(Shannon ~ log(Age) + Habitat + Site, family = gaussian(), method = "REML", data = div_df2_clean)
m_sh_alt <- gam(Shannon ~ log(Age) + Habitat + Site, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
AIC(m_sh, m_sh_alt) # same response scale, directly comparable

m_sh_s <- gam(Shannon ~ s(log(Age), k = 5) + Habitat + Site, family = gaussian(), method = "REML", data = div_df2_clean)
summary(m_sh_s)

summary(m_sh)
gam.check(m_sh)
```

## Inverted Simpson index vs age

```{r}
m_is <- gam(InvSimpson ~ log(Age) + Habitat + Site, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
m_is_ln <- gam(log(InvSimpson) ~ log(Age) + Habitat + Site, family = gaussian(), method = "REML", data = div_df2_clean)

AIC(m_is)
AIC(m_is_ln) + 2 * sum(log(div_df2_clean$InvSimpson)) # Jacobian term 

m_is_s <- gam(InvSimpson ~ s(log(Age), k = 5) + Habitat + Site, family = Gamma(link = "log"), method = "REML", data = div_df2_clean)
summary(m_is_s)

summary(m_is)
gam.check(m_is)
```

## ECM relative abundance vs age

```{r}
# ECM relative abundance: share of ECM reads per sample
ecm_reads <- colSums(filt_count_tab_no_r[ECM_OTUs_present, ])
all_reads <- colSums(filt_count_tab_no_r)
ecm_abund <- data.frame(Sample = names(ecm_reads), ECM_prop = ecm_reads / all_reads)
div_df2_clean <- merge(div_df2_clean, ecm_abund, by = "Sample")
summary(div_df2_clean$ECM_prop)

par(mfrow = c(1, 2))
hist(div_df2_clean$ECM_prop, main = "ECM share of reads")                     # left-skewed
hist(qlogis(div_df2_clean$ECM_prop), main = "logit(ECM share)")               # more symmetric
par(mfrow = c(1, 1))

# proportions between 0 and 1: beta regression, or Gaussian on the logit
m_ab <- gam(ECM_prop ~ log(Age) + Habitat + Site, family = betar(link = "logit"), method = "REML", data = div_df2_clean)
m_ab_lg <- gam(qlogis(ECM_prop) ~ log(Age) + Habitat + Site, family = gaussian(), method = "REML", data = div_df2_clean)

# comparable AICs: the logit model needs the Jacobian term
AIC(m_ab)
AIC(m_ab_lg) + 2 * sum(log(div_df2_clean$ECM_prop * (1 - div_df2_clean$ECM_prop)))

# is the age effect linear? edf close to 1 = yes
m_ab_s <- gam(ECM_prop ~ s(log(Age), k = 5) + Habitat + Site, family = betar(link = "logit"), method = "REML", data = div_df2_clean)
summary(m_ab_s)

summary(m_ab)
gam.check(m_ab)

# age effect at the plot level
nd_ab <- expand.grid(Age = seq(min(div_df2_clean$Age), max(div_df2_clean$Age), length.out = 50),
                     Habitat = unique(div_df2_clean$Habitat),
                     Site = unique(div_df2_clean$Site))
nd_ab$pred <- predict(m_ab, newdata = nd_ab, type = "response")

ggplot(div_df2_clean, aes(x = Age, y = ECM_prop)) +
  geom_point(alpha = 0.6, color = "#EE7733") +
  geom_line(data = nd_ab, aes(y = pred), color = "#009988", linewidth = 1) +
  facet_grid(Habitat ~ Site) +
  scale_x_log10() +
  theme_bw() + labs(x = "Tree age (years, log scale)", y = "ECM share of reads")
  
```

## Summary of the four models

```{r}
# estimates are on the scale of each model's link: log for richness, identity for Shannon,
# log for InvSimpson and logit for the ECM share

res <- data.frame(Response = c("Observed (log)", "Shannon", "InvSimpson", "ECM share (logit)"),
                  Estimate = c(summary(m_ln)$p.table["log(Age)", 1],
                               summary(m_sh)$p.table["log(Age)", 1],
                               summary(m_is)$p.table["log(Age)", 1],
                               summary(m_ab)$p.table["log(Age)", 1]),
                  SE = c(summary(m_ln)$p.table["log(Age)", 2],
                         summary(m_sh)$p.table["log(Age)", 2],
                         summary(m_is)$p.table["log(Age)", 2],
                         summary(m_ab)$p.table["log(Age)", 2]),
                  p = c(summary(m_ln)$p.table["log(Age)", 4],
                        summary(m_sh)$p.table["log(Age)", 4],
                        summary(m_is)$p.table["log(Age)", 4],
                        summary(m_ab)$p.table["log(Age)", 4]),
                  R2_adj = c(summary(m_ln)$r.sq, summary(m_sh)$r.sq, summary(m_is)$r.sq, summary(m_ab)$r.sq))

res$CI_low  <- res$Estimate - 1.96 * res$SE
res$CI_high <- res$Estimate + 1.96 * res$SE
res$p_adj   <- p.adjust(res$p, method = "BH")   # four tests, as elsewhere
res <- res[, c("Response", "Estimate", "SE", "CI_low", "CI_high", "R2_adj", "p", "p_adj")]

# age effect for OTU richness, at the plot level
af <- div_df2_clean[div_df2_clean$Site == "Alaska_Range" & div_df2_clean$Habitat == "forest", ]
at <- div_df2_clean[div_df2_clean$Site == "Alaska_Range" & div_df2_clean$Habitat == "treeline", ]
bf <- div_df2_clean[div_df2_clean$Site == "Brooks_Range" & div_df2_clean$Habitat == "forest", ]
bt <- div_df2_clean[div_df2_clean$Site == "Brooks_Range" & div_df2_clean$Habitat == "treeline", ]
iff <- div_df2_clean[div_df2_clean$Site == "Interior" & div_df2_clean$Habitat == "forest", ]
it <- div_df2_clean[div_df2_clean$Site == "Interior" & div_df2_clean$Habitat == "treeline", ]

# prediction range per panel, so the line is not extrapolated beyond the ages sampled there
nd <- rbind(
  data.frame(Age = seq(min(af$Age), max(af$Age), length.out = 50), Habitat = "forest",   Site = "Alaska_Range"),
  data.frame(Age = seq(min(at$Age), max(at$Age), length.out = 50), Habitat = "treeline", Site = "Alaska_Range"),
  data.frame(Age = seq(min(bf$Age), max(bf$Age), length.out = 50), Habitat = "forest",   Site = "Brooks_Range"),
  data.frame(Age = seq(min(bt$Age), max(bt$Age), length.out = 50), Habitat = "treeline", Site = "Brooks_Range"),
  data.frame(Age = seq(min(iff$Age), max(iff$Age), length.out = 50), Habitat = "forest",   Site = "Interior"),
  data.frame(Age = seq(min(it$Age), max(it$Age), length.out = 50), Habitat = "treeline", Site = "Interior"))
nd$pred <- exp(predict(m_ln, newdata = nd)) # log(richness) -> richness

# one prediction grid per plot, so the line is not drawn beyond the ages sampled there
ggplot(div_df2_clean, aes(x = Age, y = Observed)) +
  geom_point(alpha = 0.6, color = "#EE7733") +
  geom_line(data = nd, aes(y = pred), color = "#009988", linewidth = 1) +
  facet_grid(Habitat ~ Site) +
  scale_x_log10() +
  theme_bw() +
  labs(x = "Tree age (years, log scale)", y = "ECM OTU richness")

res

## Writing the table of results
res$df      <- c(df.residual(m_ln), df.residual(m_sh), df.residual(m_is), df.residual(m_ab))
res$CI_low  <- res$Estimate - qt(0.975, res$df) * res$SE
res$CI_high <- res$Estimate + qt(0.975, res$df) * res$SE
write.csv(res, "Table-ST6_div_vs_age.csv", row.names = FALSE)
```

