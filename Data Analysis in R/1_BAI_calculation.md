---
title: "BAI, detrended BAI and RCS"
author: "Saroj Basnet, Kristina Kuprina"
date: "`r Sys.Date()`"
output: html_document
---

## Three growth measures for the 2022 trees
#### 1. BAI: raw basal area increment (bai.out with the measured DBH)
#### 2. detBAI: BAI detrended with a 30-year spline fitted to log(BAI)
#### 3. rcsBAI: BAI divided by the regional curve (RCS) of its plot, built from all 2012, 2015 and 2022 cores

### Packages

```{r}
#install.packages("dplR")
library(dplR)
```

### Read tree core data: 2022 cores (breast height)

```{r}
FB <- read.rwl("IF_IT_2022.rwl", format="tucson")
# core 45201 B is misdated after ~1990 (runs to 2024): removed until re-crossdated
drop <- grepl("^45201.?B$", colnames(FB))
cat("Removed core:", colnames(FB)[drop], "\n")
FB <- FB[, !drop]
FB[!is.na(FB) & FB < 0] <- NA # negative values are missing-value codes
IDs1 <- read.ids(FB, stc=c(0,5,1))
F_B <-  treeMean(rwl=FB, ids = IDs1 , na.rm=TRUE)

DF <- read.rwl("AF_2022.rwl", format="tucson")
DF[!is.na(DF) & DF < 0] <- NA
IDs2 <- read.ids(DF, stc=c(0,5,1))
D_F <-  treeMean(rwl=DF, ids = IDs2 , na.rm=TRUE)

DT <- read.rwl("AT_2022.rwl", format="tucson")
DT[!is.na(DT) & DT < 0] <- NA
IDs3 <- read.ids(DT, stc=c(2,3,1))
D_T <-  treeMean(rwl=DT, ids = IDs3 , na.rm=TRUE)
colnames(D_T) <- paste0("PE", colnames(D_T))

BF <- read.rwl("BF_2022.rwl", format="tucson")
BF[!is.na(BF) & BF < 0] <- NA
IDs4 <- read.ids(BF, stc=c(2,3,1))
B_F <-  treeMean(rwl=BF, ids = IDs4 , na.rm=TRUE)
colnames(B_F) <- paste0("BF", colnames(B_F))

BT <- read.rwl("BT_2022.rwl", format="tucson")
BT[!is.na(BT) & BT < 0] <- NA
IDs5 <- read.ids(BT, stc=c(2,3,1))
B_T <-  treeMean(rwl=BT, ids = IDs5 , na.rm=TRUE)
colnames(B_T) <- paste0("BT", colnames(B_T))
```

### Read older cores (2012 and 2015; taken below breast height)

```{r}
Bluff_Fairbanks <- read.rwl("IF_IT_2015.rwl", format="tucson")
Denali_Forest <- read.rwl("AF_2012.rwl", format="tucson")
Denali_Treeline <- read.rwl("AT_2012.rwl", format="tucson")
Brooks_Range_Forest <- read.rwl("BF_2012.rwl", format="tucson")
Brooks_Range_Treeline <- read.rwl("BT_2012.rwl", format="tucson")
```

### Read tree DBH data (2022, cm)

```{r}
fb <- read.csv("IF_IT.csv")
df <- read.csv("AF.csv")
dt <- read.csv("AT.csv")
bf <- read.csv("BF.csv")
bt <- read.csv("BT.csv")
```

### Load tree metadata (older cores)

```{r}
metadata <- read.csv("Metadata_Alaska.csv")
selected_data <- metadata[, c("ID", "Site", "DBH_cm", "diameter_at_coring_height_A_cm", "years_to_pith")]

## data for each site
filtered_data_DF <- selected_data[grepl("Denali NP, forest", selected_data$Site), ]
filtered_data_DT <- selected_data[grepl("Denali NP, treeline", selected_data$Site), ]
filtered_data_FB <- selected_data[grepl("Fairbanks bluff", selected_data$Site), ]
filtered_data_BF <- selected_data[grepl("Nutirwik Creek, S-facing slope, forest", selected_data$Site), ]
filtered_data_BT <- selected_data[grepl("Nutirwik Creek, S-facing slope, treeline", selected_data$Site), ]
```

# BAI

## BAI calculation: Interior Alaska (AKA Fairbanks)

```{r}
# DBH in mm, in the same order as the tree columns
Bluff_Fairbanks_DBH <- fb$DBH_cm[match(colnames(F_B), as.character(fb$ID))] * 10 # cm to mm
cat("Trees in Ring-Width File:", ncol(F_B), "\n")
cat("Trees with DBH:", sum(!is.na(Bluff_Fairbanks_DBH)), "\n")
stopifnot(!anyNA(Bluff_Fairbanks_DBH)) # stops if a tree has no DBH

# 1. raw BAI with the measured DBH
BAI_Bluff_Fairbanks <- bai.out(F_B, diam = data.frame(ID = colnames(F_B), D = Bluff_Fairbanks_DBH))

# 2. BAI with radius at least the summed ring widths (keeps BAI positive down to the pith; for detrending and RCS)
Bluff_Fairbanks_Dfix <- pmax(Bluff_Fairbanks_DBH, 2 * colSums(F_B, na.rm = TRUE))
BAIfix_Bluff_Fairbanks <- bai.out(F_B, diam = data.frame(ID = colnames(F_B), D = Bluff_Fairbanks_Dfix))

# 3. detrended BAI: 30-year spline (50% cut-off) fitted to log(BAI), index = exp(log BAI - spline)
#    (BAI in µm2 = mm2 x 1e6 only keeps the log values positive; it does not change the index)
FBDET <- exp(detrend(rwl = log(BAIfix_Bluff_Fairbanks * 1e6), method = "Spline", nyrs = 30, f = 0.5, difference = TRUE))

# 4. older cores: negative values are missing-value codes, keep series with at least 10 rings
Bluff_Fairbanks[!is.na(Bluff_Fairbanks) & Bluff_Fairbanks < 0] <- NA
Bluff_Fairbanks_Old <- Bluff_Fairbanks[, colSums(!is.na(Bluff_Fairbanks)) >= 10]
Bluff_Fairbanks_Meta <- filtered_data_FB[match(colnames(Bluff_Fairbanks_Old), filtered_data_FB$ID), ]
cat("Older series:", ncol(Bluff_Fairbanks_Old), " with metadata:", sum(!is.na(Bluff_Fairbanks_Meta$ID)), "\n")

# diameter at coring height; DBH if not measured; radius at least the summed ring widths
Dc <- Bluff_Fairbanks_Meta$diameter_at_coring_height_A_cm
DBH_old <- Bluff_Fairbanks_Meta$DBH_cm
D_old <- Dc
no_Dc <- is.na(Dc) | Dc <= 0
D_old[no_Dc] <- DBH_old[no_Dc]
D_old[is.na(D_old)] <- 0
D_old <- pmax(D_old * 10, 2 * colSums(Bluff_Fairbanks_Old, na.rm = TRUE)) # cm to mm
BAI_Bluff_Fairbanks_Old <- bai.out(Bluff_Fairbanks_Old, diam = data.frame(ID = colnames(Bluff_Fairbanks_Old), D = D_old))

# coring height -> breast height with each tree's own diameters: BAI x (DBH / diameter at coring height)^2
conv <- (DBH_old / Dc)^2
conv[!is.finite(conv) | conv <= 0] <- 1 # no conversion if a diameter is missing
BAI_Bluff_Fairbanks_Old <- BAI_Bluff_Fairbanks_Old * rep(conv, each = nrow(BAI_Bluff_Fairbanks_Old))
colnames(BAI_Bluff_Fairbanks_Old) <- paste0("old_", colnames(BAI_Bluff_Fairbanks_Old)) # no matching with the 2022 trees

# pith offsets: older cores years to pith + 1 (1 if not recorded); 2022 cores 1 (no pith estimates)
Bluff_Fairbanks_po_old <- Bluff_Fairbanks_Meta$years_to_pith + 1
Bluff_Fairbanks_po_old[is.na(Bluff_Fairbanks_po_old)] <- 1

# 5. RCS with the older cores of the same site (forest and treeline pooled: the 2015 data do not separate them)
All_FB <- combine.rwl(BAIfix_Bluff_Fairbanks, BAI_Bluff_Fairbanks_Old)
FB_po <- data.frame(series = c(colnames(BAIfix_Bluff_Fairbanks), colnames(BAI_Bluff_Fairbanks_Old)),
                    pith.offset = as.integer(c(rep(1, ncol(BAIfix_Bluff_Fairbanks)), Bluff_Fairbanks_po_old)))
FB_po <- FB_po[match(colnames(All_FB), FB_po$series), ]
cat("Trees in regional curve:", ncol(All_FB), "\n")
FBRCS <- rcs(rwl = All_FB, po = FB_po, biweight = TRUE, ratios = TRUE, rc.out = TRUE, make.plot = FALSE)
```

## BAI calculation: Alaska Range (AKA Denali) forest

```{r}
# DBH in mm, in the same order as the tree columns
Denali_Forest_DBH <- df$DBH_cm[match(colnames(D_F), as.character(df$ID))] * 10 # cm to mm
cat("Trees in Ring-Width File:", ncol(D_F), "\n")
cat("Trees with DBH:", sum(!is.na(Denali_Forest_DBH)), "\n")
stopifnot(!anyNA(Denali_Forest_DBH)) # stops if a tree has no DBH

# 1. raw BAI with the measured DBH
BAI_Denali_Forest <- bai.out(D_F, diam = data.frame(ID = colnames(D_F), D = Denali_Forest_DBH))

# 2. BAI with radius at least the summed ring widths (keeps BAI positive down to the pith; for detrending and RCS)
Denali_Forest_Dfix <- pmax(Denali_Forest_DBH, 2 * colSums(D_F, na.rm = TRUE))
BAIfix_Denali_Forest <- bai.out(D_F, diam = data.frame(ID = colnames(D_F), D = Denali_Forest_Dfix))

# 3. detrended BAI: 30-year spline (50% cut-off) fitted to log(BAI), index = exp(log BAI - spline)
#    (BAI in µm2 = mm2 x 1e6 only keeps the log values positive; it does not change the index)
DFDET <- exp(detrend(rwl = log(BAIfix_Denali_Forest * 1e6), method = "Spline", nyrs = 30, f = 0.5, difference = TRUE))

# 4. older cores: negative values are missing-value codes, keep series with at least 10 rings
Denali_Forest[!is.na(Denali_Forest) & Denali_Forest < 0] <- NA
Denali_Forest_Old <- Denali_Forest[, colSums(!is.na(Denali_Forest)) >= 10]
Denali_Forest_Meta <- filtered_data_DF[match(colnames(Denali_Forest_Old), filtered_data_DF$ID), ]
cat("Older series:", ncol(Denali_Forest_Old), " with metadata:", sum(!is.na(Denali_Forest_Meta$ID)), "\n")

# diameter at coring height; DBH if not measured; radius at least the summed ring widths
Dc <- Denali_Forest_Meta$diameter_at_coring_height_A_cm
DBH_old <- Denali_Forest_Meta$DBH_cm
D_old <- Dc
no_Dc <- is.na(Dc) | Dc <= 0
D_old[no_Dc] <- DBH_old[no_Dc]
D_old[is.na(D_old)] <- 0
D_old <- pmax(D_old * 10, 2 * colSums(Denali_Forest_Old, na.rm = TRUE)) # cm to mm
BAI_Denali_Forest_Old <- bai.out(Denali_Forest_Old, diam = data.frame(ID = colnames(Denali_Forest_Old), D = D_old))

# coring height -> breast height with each tree's own diameters: BAI x (DBH / diameter at coring height)^2
conv <- (DBH_old / Dc)^2
conv[!is.finite(conv) | conv <= 0] <- 1 # no conversion if a diameter is missing
BAI_Denali_Forest_Old <- BAI_Denali_Forest_Old * rep(conv, each = nrow(BAI_Denali_Forest_Old))
colnames(BAI_Denali_Forest_Old) <- paste0("old_", colnames(BAI_Denali_Forest_Old)) # no matching with the 2022 trees

# pith offsets: older cores years to pith + 1 (1 if not recorded); 2022 cores 1 (no pith estimates)
Denali_Forest_po_old <- Denali_Forest_Meta$years_to_pith + 1
Denali_Forest_po_old[is.na(Denali_Forest_po_old)] <- 1

# 5. RCS with the older cores of the same plot
All_DF <- combine.rwl(BAIfix_Denali_Forest, BAI_Denali_Forest_Old)
DF_po <- data.frame(series = c(colnames(BAIfix_Denali_Forest), colnames(BAI_Denali_Forest_Old)),
                    pith.offset = as.integer(c(rep(1, ncol(BAIfix_Denali_Forest)), Denali_Forest_po_old)))
DF_po <- DF_po[match(colnames(All_DF), DF_po$series), ]
cat("Trees in regional curve:", ncol(All_DF), "\n")
DFRCS <- rcs(rwl = All_DF, po = DF_po, biweight = TRUE, ratios = TRUE, rc.out = TRUE, make.plot = FALSE)
```

## BAI calculation: Alaska Range (AKA Denali) treeline

```{r}
# DBH in mm, in the same order as the tree columns
Denali_Treeline_DBH <- dt$DBH_cm[match(colnames(D_T), as.character(dt$ID))] * 10 # cm to mm
cat("Trees in Ring-Width File:", ncol(D_T), "\n")
cat("Trees with DBH:", sum(!is.na(Denali_Treeline_DBH)), "\n")
stopifnot(!anyNA(Denali_Treeline_DBH)) # stops if a tree has no DBH

# 1. raw BAI with the measured DBH
BAI_Denali_Treeline <- bai.out(D_T, diam = data.frame(ID = colnames(D_T), D = Denali_Treeline_DBH))

# 2. BAI with radius at least the summed ring widths (keeps BAI positive down to the pith; for detrending and RCS)
Denali_Treeline_Dfix <- pmax(Denali_Treeline_DBH, 2 * colSums(D_T, na.rm = TRUE))
BAIfix_Denali_Treeline <- bai.out(D_T, diam = data.frame(ID = colnames(D_T), D = Denali_Treeline_Dfix))

# 3. detrended BAI: 30-year spline (50% cut-off) fitted to log(BAI), index = exp(log BAI - spline)
#    (BAI in µm2 = mm2 x 1e6 only keeps the log values positive; it does not change the index)
DTDET <- exp(detrend(rwl = log(BAIfix_Denali_Treeline * 1e6), method = "Spline", nyrs = 30, f = 0.5, difference = TRUE))

# 4. older cores: negative values are missing-value codes, keep series with at least 10 rings
Denali_Treeline[!is.na(Denali_Treeline) & Denali_Treeline < 0] <- NA
Denali_Treeline_Old <- Denali_Treeline[, colSums(!is.na(Denali_Treeline)) >= 10]
Denali_Treeline_Meta <- filtered_data_DT[match(colnames(Denali_Treeline_Old), filtered_data_DT$ID), ]
cat("Older series:", ncol(Denali_Treeline_Old), " with metadata:", sum(!is.na(Denali_Treeline_Meta$ID)), "\n")

# diameter at coring height; DBH if not measured; radius at least the summed ring widths
Dc <- Denali_Treeline_Meta$diameter_at_coring_height_A_cm
DBH_old <- Denali_Treeline_Meta$DBH_cm
D_old <- Dc
no_Dc <- is.na(Dc) | Dc <= 0
D_old[no_Dc] <- DBH_old[no_Dc]
D_old[is.na(D_old)] <- 0
D_old <- pmax(D_old * 10, 2 * colSums(Denali_Treeline_Old, na.rm = TRUE)) # cm to mm
BAI_Denali_Treeline_Old <- bai.out(Denali_Treeline_Old, diam = data.frame(ID = colnames(Denali_Treeline_Old), D = D_old))

# coring height -> breast height with each tree's own diameters: BAI x (DBH / diameter at coring height)^2
conv <- (DBH_old / Dc)^2
conv[!is.finite(conv) | conv <= 0] <- 1 # no conversion if a diameter is missing
BAI_Denali_Treeline_Old <- BAI_Denali_Treeline_Old * rep(conv, each = nrow(BAI_Denali_Treeline_Old))
colnames(BAI_Denali_Treeline_Old) <- paste0("old_", colnames(BAI_Denali_Treeline_Old)) # no matching with the 2022 trees

# pith offsets: older cores years to pith + 1 (1 if not recorded); 2022 cores 1 (no pith estimates)
Denali_Treeline_po_old <- Denali_Treeline_Meta$years_to_pith + 1
Denali_Treeline_po_old[is.na(Denali_Treeline_po_old)] <- 1

# 5. RCS with the older cores of the same plot
All_DT <- combine.rwl(BAIfix_Denali_Treeline, BAI_Denali_Treeline_Old)
DT_po <- data.frame(series = c(colnames(BAIfix_Denali_Treeline), colnames(BAI_Denali_Treeline_Old)),
                    pith.offset = as.integer(c(rep(1, ncol(BAIfix_Denali_Treeline)), Denali_Treeline_po_old)))
DT_po <- DT_po[match(colnames(All_DT), DT_po$series), ]
cat("Trees in regional curve:", ncol(All_DT), "\n")
DTRCS <- rcs(rwl = All_DT, po = DT_po, biweight = TRUE, ratios = TRUE, rc.out = TRUE, make.plot = FALSE)
```

## BAI calculation: Brooks Range forest

```{r}
# DBH in mm, in the same order as the tree columns
Brooks_Range_Forest_DBH <- bf$DBH_cm[match(colnames(B_F), as.character(bf$ID))] * 10 # cm to mm
cat("Trees in Ring-Width File:", ncol(B_F), "\n")
cat("Trees with DBH:", sum(!is.na(Brooks_Range_Forest_DBH)), "\n")
stopifnot(!anyNA(Brooks_Range_Forest_DBH)) # stops if a tree has no DBH

# 1. raw BAI with the measured DBH
BAI_Brooks_Range_Forest <- bai.out(B_F, diam = data.frame(ID = colnames(B_F), D = Brooks_Range_Forest_DBH))

# 2. BAI with radius at least the summed ring widths (keeps BAI positive down to the pith; for detrending and RCS)
Brooks_Range_Forest_Dfix <- pmax(Brooks_Range_Forest_DBH, 2 * colSums(B_F, na.rm = TRUE))
BAIfix_Brooks_Range_Forest <- bai.out(B_F, diam = data.frame(ID = colnames(B_F), D = Brooks_Range_Forest_Dfix))

# 3. detrended BAI: 30-year spline (50% cut-off) fitted to log(BAI), index = exp(log BAI - spline)
#    (BAI in µm2 = mm2 x 1e6 only keeps the log values positive; it does not change the index)
BFDET <- exp(detrend(rwl = log(BAIfix_Brooks_Range_Forest * 1e6), method = "Spline", nyrs = 30, f = 0.5, difference = TRUE))

# 4. older cores: negative values are missing-value codes, keep series with at least 10 rings
Brooks_Range_Forest[!is.na(Brooks_Range_Forest) & Brooks_Range_Forest < 0] <- NA
Brooks_Range_Forest_Old <- Brooks_Range_Forest[, colSums(!is.na(Brooks_Range_Forest)) >= 10]
Brooks_Range_Forest_Meta <- filtered_data_BF[match(colnames(Brooks_Range_Forest_Old), filtered_data_BF$ID), ]
cat("Older series:", ncol(Brooks_Range_Forest_Old), " with metadata:", sum(!is.na(Brooks_Range_Forest_Meta$ID)), "\n")

# diameter at coring height; DBH if not measured; radius at least the summed ring widths
Dc <- Brooks_Range_Forest_Meta$diameter_at_coring_height_A_cm
DBH_old <- Brooks_Range_Forest_Meta$DBH_cm
D_old <- Dc
no_Dc <- is.na(Dc) | Dc <= 0
D_old[no_Dc] <- DBH_old[no_Dc]
D_old[is.na(D_old)] <- 0
D_old <- pmax(D_old * 10, 2 * colSums(Brooks_Range_Forest_Old, na.rm = TRUE)) # cm to mm
BAI_Brooks_Range_Forest_Old <- bai.out(Brooks_Range_Forest_Old, diam = data.frame(ID = colnames(Brooks_Range_Forest_Old), D = D_old))

# coring height -> breast height with each tree's own diameters: BAI x (DBH / diameter at coring height)^2
conv <- (DBH_old / Dc)^2
conv[!is.finite(conv) | conv <= 0] <- 1 # no conversion if a diameter is missing
BAI_Brooks_Range_Forest_Old <- BAI_Brooks_Range_Forest_Old * rep(conv, each = nrow(BAI_Brooks_Range_Forest_Old))
colnames(BAI_Brooks_Range_Forest_Old) <- paste0("old_", colnames(BAI_Brooks_Range_Forest_Old)) # no matching with the 2022 trees

# pith offsets: older cores years to pith + 1 (1 if not recorded); 2022 cores 1 (no pith estimates)
Brooks_Range_Forest_po_old <- Brooks_Range_Forest_Meta$years_to_pith + 1
Brooks_Range_Forest_po_old[is.na(Brooks_Range_Forest_po_old)] <- 1

# 5. RCS with the older cores of the same plot
All_BF <- combine.rwl(BAIfix_Brooks_Range_Forest, BAI_Brooks_Range_Forest_Old)
BF_po <- data.frame(series = c(colnames(BAIfix_Brooks_Range_Forest), colnames(BAI_Brooks_Range_Forest_Old)),
                    pith.offset = as.integer(c(rep(1, ncol(BAIfix_Brooks_Range_Forest)), Brooks_Range_Forest_po_old)))
BF_po <- BF_po[match(colnames(All_BF), BF_po$series), ]
cat("Trees in regional curve:", ncol(All_BF), "\n")
BFRCS <- rcs(rwl = All_BF, po = BF_po, biweight = TRUE, ratios = TRUE, rc.out = TRUE, make.plot = FALSE)
```

## BAI calculation: Brooks Range treeline

```{r}
# DBH in mm, in the same order as the tree columns
Brooks_Range_Treeline_DBH <- bt$DBH_cm[match(colnames(B_T), as.character(bt$ID))] * 10 # cm to mm
cat("Trees in Ring-Width File:", ncol(B_T), "\n")
cat("Trees with DBH:", sum(!is.na(Brooks_Range_Treeline_DBH)), "\n")
stopifnot(!anyNA(Brooks_Range_Treeline_DBH)) # stops if a tree has no DBH

# 1. raw BAI with the measured DBH
BAI_Brooks_Range_Treeline <- bai.out(B_T, diam = data.frame(ID = colnames(B_T), D = Brooks_Range_Treeline_DBH))

# 2. BAI with radius at least the summed ring widths (keeps BAI positive down to the pith; for detrending and RCS)
Brooks_Range_Treeline_Dfix <- pmax(Brooks_Range_Treeline_DBH, 2 * colSums(B_T, na.rm = TRUE))
BAIfix_Brooks_Range_Treeline <- bai.out(B_T, diam = data.frame(ID = colnames(B_T), D = Brooks_Range_Treeline_Dfix))

# 3. detrended BAI: 30-year spline (50% cut-off) fitted to log(BAI), index = exp(log BAI - spline)
#    (BAI in µm2 = mm2 x 1e6 only keeps the log values positive; it does not change the index)
BTDET <- exp(detrend(rwl = log(BAIfix_Brooks_Range_Treeline * 1e6), method = "Spline", nyrs = 30, f = 0.5, difference = TRUE))

# 4. older cores: negative values are missing-value codes, keep series with at least 10 rings
Brooks_Range_Treeline[!is.na(Brooks_Range_Treeline) & Brooks_Range_Treeline < 0] <- NA
Brooks_Range_Treeline_Old <- Brooks_Range_Treeline[, colSums(!is.na(Brooks_Range_Treeline)) >= 10]
Brooks_Range_Treeline_Meta <- filtered_data_BT[match(colnames(Brooks_Range_Treeline_Old), filtered_data_BT$ID), ]
cat("Older series:", ncol(Brooks_Range_Treeline_Old), " with metadata:", sum(!is.na(Brooks_Range_Treeline_Meta$ID)), "\n")

# diameter at coring height; DBH if not measured; radius at least the summed ring widths
Dc <- Brooks_Range_Treeline_Meta$diameter_at_coring_height_A_cm
DBH_old <- Brooks_Range_Treeline_Meta$DBH_cm
D_old <- Dc
no_Dc <- is.na(Dc) | Dc <= 0
D_old[no_Dc] <- DBH_old[no_Dc]
D_old[is.na(D_old)] <- 0
D_old <- pmax(D_old * 10, 2 * colSums(Brooks_Range_Treeline_Old, na.rm = TRUE)) # cm to mm
BAI_Brooks_Range_Treeline_Old <- bai.out(Brooks_Range_Treeline_Old, diam = data.frame(ID = colnames(Brooks_Range_Treeline_Old), D = D_old))

# coring height -> breast height with each tree's own diameters: BAI x (DBH / diameter at coring height)^2
conv <- (DBH_old / Dc)^2
conv[!is.finite(conv) | conv <= 0] <- 1 # no conversion if a diameter is missing
BAI_Brooks_Range_Treeline_Old <- BAI_Brooks_Range_Treeline_Old * rep(conv, each = nrow(BAI_Brooks_Range_Treeline_Old))
colnames(BAI_Brooks_Range_Treeline_Old) <- paste0("old_", colnames(BAI_Brooks_Range_Treeline_Old)) # no matching with the 2022 trees

# pith offsets: older cores years to pith + 1 (1 if not recorded); 2022 cores 1 (no pith estimates)
Brooks_Range_Treeline_po_old <- Brooks_Range_Treeline_Meta$years_to_pith + 1
Brooks_Range_Treeline_po_old[is.na(Brooks_Range_Treeline_po_old)] <- 1

# 5. RCS with the older cores of the same plot
All_BT <- combine.rwl(BAIfix_Brooks_Range_Treeline, BAI_Brooks_Range_Treeline_Old)
BT_po <- data.frame(series = c(colnames(BAIfix_Brooks_Range_Treeline), colnames(BAI_Brooks_Range_Treeline_Old)),
                    pith.offset = as.integer(c(rep(1, ncol(BAIfix_Brooks_Range_Treeline)), Brooks_Range_Treeline_po_old)))
BT_po <- BT_po[match(colnames(All_BT), BT_po$series), ]
cat("Trees in regional curve:", ncol(All_BT), "\n")
BTRCS <- rcs(rwl = All_BT, po = BT_po, biweight = TRUE, ratios = TRUE, rc.out = TRUE, make.plot = FALSE)
```

# Regional curves

```{r}
par(mfrow = c(2, 3))
plot(FBRCS$rc, type = "l", main = "Interior Alaska", xlab = "Cambial age (years)", ylab = "BAI (mm2)")
plot(DFRCS$rc, type = "l", main = "Alaska Range forest", xlab = "Cambial age (years)", ylab = "BAI (mm2)")
plot(DTRCS$rc, type = "l", main = "Alaska Range treeline", xlab = "Cambial age (years)", ylab = "BAI (mm2)")
plot(BFRCS$rc, type = "l", main = "Brooks Range forest", xlab = "Cambial age (years)", ylab = "BAI (mm2)")
plot(BTRCS$rc, type = "l", main = "Brooks Range treeline", xlab = "Cambial age (years)", ylab = "BAI (mm2)")
par(mfrow = c(1, 1))
```

### Combine all trees and calculate the means for the last 5, 10, 15 and 30 years

```{r}
BAI_all <- combine.rwl(list(BAI_Bluff_Fairbanks, BAI_Denali_Forest, BAI_Denali_Treeline,
                            BAI_Brooks_Range_Forest, BAI_Brooks_Range_Treeline))
DET_all <- combine.rwl(list(FBDET, DFDET, DTDET, BFDET, BTDET))
# RCS: keep only the 2022 trees (the older cores were used for the curves only)
RCS_all <- combine.rwl(list(FBRCS$rwi[, colnames(F_B)], DFRCS$rwi[, colnames(D_F)], DTRCS$rwi[, colnames(D_T)],
                            BFRCS$rwi[, colnames(B_F)], BTRCS$rwi[, colnames(B_T)]))

info <- read.table("Alaska_info_no_R.txt", row.names=1, header=T, check.names=F)

# sample names -> tree names in the ring-width files
tree <- sub("-.*$", "", rownames(info))   # BF278-1 -> BF278
tree <- sub("^DF", "", tree)              # DF45567 -> 45567
pref <- grepl("^(PE|BF|BT)", tree)
tree[pref] <- paste0(substr(tree[pref], 1, 2), as.integer(substr(tree[pref], 3, nchar(tree[pref])))) # PE035 -> PE35

# one column per sample (in the order of Alaska_info_no_R.txt); NA for samples without cores
BAI_m <- as.matrix(BAI_all)[, match(tree, colnames(BAI_all))] / 100 # mm2 -> cm2
DET_m <- as.matrix(DET_all)[, match(tree, colnames(DET_all))]
RCS_m <- as.matrix(RCS_all)[, match(tree, colnames(RCS_all))]
colnames(BAI_m) <- rownames(info)
colnames(DET_m) <- rownames(info)
colnames(RCS_m) <- rownames(info)

res <- data.frame(Sample = rownames(info), Plot = info$Plot,
                  BAI_5y  = colMeans(BAI_m[as.character(2018:2022), ], na.rm = TRUE),
                  BAI_10y = colMeans(BAI_m[as.character(2013:2022), ], na.rm = TRUE),
                  BAI_15y = colMeans(BAI_m[as.character(2008:2022), ], na.rm = TRUE),
                  BAI_30y = colMeans(BAI_m[as.character(1993:2022), ], na.rm = TRUE),
                  detBAI_5y  = colMeans(DET_m[as.character(2018:2022), ], na.rm = TRUE),
                  detBAI_10y = colMeans(DET_m[as.character(2013:2022), ], na.rm = TRUE),
                  detBAI_15y = colMeans(DET_m[as.character(2008:2022), ], na.rm = TRUE),
                  detBAI_30y = colMeans(DET_m[as.character(1993:2022), ], na.rm = TRUE),
                  rcsBAI_5y  = colMeans(RCS_m[as.character(2018:2022), ], na.rm = TRUE),
                  rcsBAI_10y = colMeans(RCS_m[as.character(2013:2022), ], na.rm = TRUE),
                  rcsBAI_15y = colMeans(RCS_m[as.character(2008:2022), ], na.rm = TRUE),
                  rcsBAI_30y = colMeans(RCS_m[as.character(1993:2022), ], na.rm = TRUE))
res[is.na(res)] <- NA # NaN (samples without cores) -> NA

cat("Samples with BAI:", sum(!is.na(res$BAI_5y)), "of", nrow(res), "\n")

# check against the values in Alaska_info_no_R.txt (differences should be ~0)
summary(res$BAI_5y - info$BAI_5y)
summary(res$detBAI_5y - info$detBAI_5y)
```

## Write the output

```{r}
# annual values 1990-2022 per sample (for ST1)
yrs <- as.character(1990:2022)
write.csv(t(BAI_m[yrs, ]), "BAI_annual.csv") # cm2
write.csv(t(DET_m[yrs, ]), "detBAI_annual.csv")
write.csv(t(RCS_m[yrs, ]), "rcsBAI_annual.csv")
write.csv(res, "BAI_results.csv", row.names = FALSE)

# Alaska_info_no_R with all three growth measures
info_new <- info
info_new[, colnames(res)[-(1:2)]] <- res[, -(1:2)]
write.table(cbind(Sample = rownames(info_new), info_new), "Alaska_info_no_R_BAI.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
```
