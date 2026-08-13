
#rm(list=ls())#if to delete files from environment####
rm(list=ls())
#####################################
install.packages("dplR")
install.packages("detrendeR")
install.packages("treeclim")
install.packages("ggplot2")
install.packages("dendRolAB")
library(ggplot2)
library(dplR)
library(detrendeR)
library(treeclim)
library(dendRolAB)
library(DendroEco)
library(tidyr)
install.packages("readxl")
library(readxl)
install.packages("dplyr")
library(dplyr)
#install.packages("rlang", type = "binary")
install.packages("devtools")
devtools::install_github("AllanBuras/dendRolAB")

#-------------###READ THE DATA### -------------------------
setwd("C:/Users/Samsung/Desktop/R-Alaska")

getwd()
list.files()

######### Load the data## 
FB <- read.rwl("FB_2022.rwl", format="tucson")
IDs1 <- read.ids(FB, stc=c(0,5,1))
F_B <-  treeMean(rwl=FB, ids = IDs1 , na.rm=TRUE)

DF <- read.rwl("DF_2022.rwl", format="tucson")
IDs2 <- read.ids(DF, stc=c(0,5,1))
D_F<-  treeMean(rwl=DF, ids = IDs2 , na.rm=TRUE)

DT<- read.rwl("DT_2022.rwl", format="tucson")
IDs3 <- read.ids(DT, stc=c(2,3,1))
D_T<-  treeMean(rwl=DT, ids = IDs3 , na.rm=TRUE)
colnames(D_T) <- paste0("PE", colnames(D_T))

BF <- read.rwl("BF_2022.rwl", format="tucson")
IDs4 <- read.ids(BF, stc=c(2,3,1))
B_F<-  treeMean(rwl=BF, ids = IDs4 , na.rm=TRUE)
colnames(B_F) <- paste0("BF", colnames(B_F))

BT<- read.rwl("BT_2022.rwl", format="tucson")
IDs5 <- read.ids(BT, stc=c(2,3,1))
B_T<-  treeMean(rwl=BT, ids = IDs5 , na.rm=TRUE)
colnames(B_T) <- paste0("BT", colnames(B_T))

################################################################
Bluff_Fairbanks <- read.rwl("IF_IT_2012.rwl", format="tucson")

Denali_Forest <- read.rwl("AF_2012.rwl", format="tucson")
Denali_Treeline <- read.rwl("AT_2012.rwl", format="tucson")

Brooks_Range_Forest <- read.rwl("BF_2012.rwl", format="tucson")
Brooks_Range_Treeline <- read.rwl("BT_2012.rwl", format="tucson")

########################Combine datasets#########################
All_FB <- combine.rwl(F_B, Bluff_Fairbanks)
All_DF <- combine.rwl(D_F, Denali_Forest)
All_DT <- combine.rwl(D_T, Denali_Treeline)
All_BF <- combine.rwl(B_F, Brooks_Range_Forest)
All_BT<- combine.rwl(B_T, Brooks_Range_Treeline)

#####################tree DBH datatree DBH data########################################
fb<- read.csv("FB.csv",sep=",", dec=",")
df <- read.csv("DF.csv",sep=",", dec=",")
dt<- read.csv("DT.csv",sep=",", dec=",")
bf<- read.csv("BF.csv",sep=",", dec=",")
bt<- read.csv("BT.csv",sep=",", dec=",")

###################################################################
## Load tree DBH data
metadata <- read_excel("Metadata.xlsx")
metadata
selected_data <- metadata[, c("ID","Site", "ID_Tag", "DBH_cm")]

#################data for each site#######################################
library(magrittr)
filtered_data_DF <- selected_data %>%
  filter(grepl("Denali NP, forest",Site) | grepl("Denali NP, forest, belt",Site))
filtered_data_DT<- selected_data %>%
  filter(grepl("Denali NP, treeline", Site))
filtered_data_FB <- selected_data %>%
  filter(grepl("Fairbanks bluff", Site))
filtered_data_BF <- selected_data %>%
  filter(grepl("Nutirwik Creek, S-facing slope, forest", Site))
filtered_data_BT <- selected_data %>%
  filter(grepl("Nutirwik Creek, S-facing slope, treeline", Site))

##########Data which have greater than 0,01 distance to pith value
filtered_data_DF1 <- filtered_data_DF %>%
 filter(DBH_cm > 0.01)
filtered_data_DT1 <- filtered_data_DT %>%
filter(DBH_cm > 0.01)
filtered_data_FB1 <- filtered_data_FB %>%
  filter(DBH_cm > 0.01)
filtered_data_BF1 <- filtered_data_BF %>%
  filter(DBH_cm > 0.01)
filtered_data_BT1 <- filtered_data_BT %>%
  filter(DBH_cm > 0.01)

#####################################################
Denali_Forest_DTP <- filtered_data_DF1 %>% select(ID,DBH_cm)
Denali_Treeline_DTP <- filtered_data_DT1 %>% select(ID,DBH_cm)
Brooks_Range_Forest_DTP <- filtered_data_BF1%>% select(ID,DBH_cm)
Brooks_Range_Treeline_DTP <- filtered_data_BT1 %>% select(ID,DBH_cm)
Bluff_Fairbanks_DTP <-filtered_data_FB1 %>% select(ID,DBH_cm)

##################################################################
# Ensure 'fb' has numeric values (if read.csv didn't handle dec="," correctly)
fb_clean <- fb %>%
  mutate(DBH_cm = as.numeric(gsub(",", ".", DBH_cm))) %>%
  select(ID, DBH_cm) # Keep only the columns that match Bluff_Fairbanks_DTP

# 2. Ensure IDs are the same type (Character) in both
fb_clean$ID <- as.character(fb_clean$ID)
Bluff_Fairbanks_DTP$ID <- as.character(Bluff_Fairbanks_DTP$ID)
#2. to bring both datasets together/ full_join 
#This creates DBH_cm.x (from fb_clean) and DBH_cm.y (from DTP)
Combined_Data <- full_join(fb_clean, Bluff_Fairbanks_DTP, by = "ID")

# 3. Create the final DBH column: 
# "If fb_clean has a value, use it. Otherwise, use the DTP value."
Bluff_Fairbanks_DBH <- Combined_Data %>%
  mutate(DBH_cm = coalesce(DBH_cm.x, DBH_cm.y)) %>%
  select(ID, DBH_cm)
rownames(Bluff_Fairbanks_DBH) <- Bluff_Fairbanks_DBH$ID

##################################################################
# 1. Ensure 'df' has numeric values (if read.csv didn't handle dec="," correctly)
df_clean <- df %>%
  mutate(DBH_cm = as.numeric(gsub(",", ".", DBH_cm))) %>%
  select(ID, DBH_cm) # Keep only the columns that match Bluff_Fairbanks_DTP

# Ensure IDs are the same type (Character) in both
df_clean$ID <- as.character(df_clean$ID)
Denali_Forest_DTP$ID <- as.character(Denali_Forest_DTP$ID)
# to bring both datasets together
# This creates DBH_cm.x (from fb_clean) and DBH_cm.y (from DTP)
Combined_Data <- full_join(df_clean, Denali_Forest_DTP, by = "ID")

# Create the final DBH column: 
# "If fb_clean has a value, use it. Otherwise, use the DTP value."
Denali_Forest_DBH <- Combined_Data %>%
  mutate(DBH_cm = coalesce(DBH_cm.x, DBH_cm.y)) %>%
  select(ID, DBH_cm)
rownames(Denali_Forest_DBH) <- Denali_Forest_DBH$ID


##################################################################
#  Ensure 'dT' has numeric values (if read.csv didn't handle dec="," correctly)
dt_clean <- dt %>%
  mutate(DBH_cm = as.numeric(gsub(",", ".", DBH_cm))) %>%
  select(ID, DBH_cm) # Keep only the columns that match Bluff_Fairbanks_DTP

# Ensure IDs are the same type (Character) in both
dt_clean$ID <- as.character(dt_clean$ID)
Denali_Treeline_DTP$ID <- as.character(Denali_Treeline_DTP$ID)
# Use a full_join to bring both datasets together
# This creates DBH_cm.x (from fb_clean) and DBH_cm.y (from DTP)
Combined_Data <- full_join(dt_clean, Denali_Treeline_DTP, by = "ID")

#  Create the final DBH column: 
# "If fb_clean has a value, use it. Otherwise, use the DTP value."
Denali_Treeline_DBH <- Combined_Data %>%
  mutate(DBH_cm = coalesce(DBH_cm.x, DBH_cm.y)) %>%
  select(ID, DBH_cm)
rownames(Denali_Treeline_DBH) <- Denali_Treeline_DBH$ID

##################################################################
#  Ensure 'bf' has numeric values (if read.csv didn't handle dec="," correctly)
bf_clean <- bf %>%
  mutate(DBH_cm = as.numeric(gsub(",", ".", DBH_cm))) %>%
  select(ID, DBH_cm) # Keep only the columns that match Bluff_Fairbanks_DTP

# 2. Ensure IDs are the same type (Character) in both
bf_clean$ID <- as.character(bf_clean$ID)
Brooks_Range_Forest_DTP$ID <- as.character(Brooks_Range_Forest_DTP$ID)
#2. Use a full_join to bring both datasets together
# This creates DBH_cm.x (from fb_clean) and DBH_cm.y (from DTP)
Combined_Data <- full_join(bf_clean, Brooks_Range_Forest_DTP, by = "ID")

# Create the final DBH column: 
# "If fb_clean has a value, use it. Otherwise, use the DTP value."
Brooks_Range_Forest_DBH <- Combined_Data %>%
  mutate(DBH_cm = coalesce(DBH_cm.x, DBH_cm.y)) %>%
  select(ID, DBH_cm)
rownames(Brooks_Range_Forest_DBH) <- Brooks_Range_Forest_DBH$ID


##################################################################
#  Ensure 'bt' has numeric values (if read.csv didn't handle dec="," correctly)
bt_clean <- bt %>%
  mutate(DBH_cm = as.numeric(gsub(",", ".", DBH_cm))) %>%
  select(ID, DBH_cm) # Keep only the columns that match Bluff_Fairbanks_DTP

# Ensure IDs are the same type (Character) in both
bt_clean$ID <- as.character(bt_clean$ID)
Brooks_Range_Treeline_DTP$ID <- as.character(Brooks_Range_Treeline_DTP$ID)
#Use a full_join to bring both datasets together
# This creates DBH_cm.x (from fb_clean) and DBH_cm.y (from DTP)
Combined_Data <- full_join(bt_clean, Brooks_Range_Treeline_DTP, by = "ID")

# Create the final DBH column: 
# "If fb_clean has a value, use it. Otherwise, use the DTP value."
Brooks_Range_Treeline_DBH <- Combined_Data %>%
  mutate(DBH_cm = coalesce(DBH_cm.x, DBH_cm.y)) %>%
  select(ID, DBH_cm)
rownames(Brooks_Range_Treeline_DBH) <- Brooks_Range_Treeline_DBH$ID


##############################################
# create data frame with tree IDs and diameters at coring height
is.numeric(Bluff_Fairbanks_DBH$DBH_cm)
sapply(Bluff_Fairbanks_DBH, class)
B <- transform(Bluff_Fairbanks_DBH, DBH = as.numeric(DBH_cm))
sapply(B, class)
B[,2]<-B[,2]*10# # conversion in mm# # # cm to mm conversion
rownames(B)<-B[,1]# # assure the correct assignment# # 
head(B)
B
B<-B[colnames(All_FB ),]#using the column of AV.RWJF to order A
Bluff_Fairbanks_DBH <- B[ ,c("ID","DBH_cm","DBH")]
years <- as.numeric(rownames(All_FB ))  # years in the tree-ring data

# Ensure your diameter data frame has NO missing values
Bluff_Fairbanks_DTP <- Bluff_Fairbanks_DBH %>% filter(!is.na(ID) & !is.na(DBH_cm))

#Force both sets of IDs to be character strings to avoid matching errors
tree_names <- as.character(colnames(All_FB ))
Bluff_Fairbanks_DBH$ID <- as.character(Bluff_Fairbanks_DBH$ID)

# Create the clean diameter table by specifically looking for those tree_names
# This ensures the order is IDENTICAL to the columns in Denali_Forest
Bluff_Fairbanks_DBH_Clean <- Bluff_Fairbanks_DTP[Bluff_Fairbanks_DBH$ID %in% tree_names, ]

# Sort it to match the column order exactly
Bluff_Fairbanks_DBH_Clean <- Bluff_Fairbanks_DBH_Clean[match(tree_names, Bluff_Fairbanks_DBH_Clean$ID), ]

# Check the numbers before running
cat("Trees in Ring-Width File:", ncol(All_FB), "\n")
cat("Trees in Diameter Table:", nrow(Bluff_Fairbanks_DBH_Clean), "\n")

#Find IDs that exist in BOTH the columns of Denali_Forest AND the ID column of DTP
common_trees <- intersect(colnames(All_FB), Bluff_Fairbanks_DBH$ID)

# Subset the ring-width data to keep only those common trees
Bluff_Fairbanks_Sub <- All_FB[, common_trees]

# Subset and reorder the diameter table to match exactly
Bluff_Fairbanks_DBH_Clean <- Bluff_Fairbanks_DBH[match(common_trees, Bluff_Fairbanks_DBH$ID), ]

# Final check: verify there are no NAs in the IDs
if(any(is.na(Bluff_Fairbanks_DBH_Clean$ID))) {
  stop("There are still NAs in your ID column!")
}

# Run the calculation using the new 'Sub' object
BAI_Bluff_Fairbanks <- bai.out(Bluff_Fairbanks_Sub, diam = Bluff_Fairbanks_DBH_Clean)

#################################################################
##############################################
##############################################
# create data frame with tree IDs and diameters at coring height
is.numeric(Denali_Forest_DBH$DBH_cm)
sapply(Denali_Forest_DBH, class)
B <- transform(Denali_Forest_DBH, DBH = as.numeric(DBH_cm))
sapply(B, class)
B[,2]<-B[,2]*10# # conversion in mm# # # cm to mm conversion
rownames(B)<-B[,1]# # assure the correct assignment# # 
head(B)
B
B<-B[colnames(All_DF ),]#using the column of AV.RWJF to order A
Denali_Forest_DBH <- B[ ,c("ID","DBH_cm","DBH")]
years <- as.numeric(rownames(All_DF ))  # years in the tree-ring data

# Ensure your diameter data frame has NO missing values
Denali_Forest_DTP <- Denali_Forest_DBH %>% filter(!is.na(ID) & !is.na(DBH_cm))

#Force both sets of IDs to be character strings to avoid matching errors
tree_names <- as.character(colnames(All_DF ))
Denali_Forest_DBH$ID <- as.character(Denali_Forest_DBH$ID)

#  Create the clean diameter table by specifically looking for those tree_names
# This ensures the order is IDENTICAL to the columns in Denali_Forest
Denali_Forest_DBH_Clean <- Denali_Forest_DTP[Denali_Forest_DBH$ID %in% tree_names, ]

#  Sort it to match the column order exactly
Denali_Forest_DBH_Clean <- Denali_Forest_DBH_Clean[match(tree_names, Denali_Forest_DBH_Clean$ID), ]

# 4. Check the numbers before running
cat("Trees in Ring-Width File:", ncol(All_DF), "\n")
cat("Trees in Diameter Table:", nrow(Denali_Forest_DBH_Clean), "\n")

#. Find IDs that exist in BOTH the columns of Denali_Forest AND the ID column of DTP
common_trees <- intersect(colnames(All_DF), Denali_Forest_DBH$ID)

# Subset the ring-width data to keep only those common trees
Denali_Forest_Sub <- All_DF[, common_trees]

#Subset and reorder the diameter table to match exactly
Denali_Forest_DBH_Clean <- Denali_Forest_DBH[match(common_trees, Denali_Forest_DBH$ID), ]

#Final check: verify there are no NAs in the IDs
if(any(is.na(Denali_Forest_DBH_Clean$ID))) {
  stop("There are still NAs in your ID column!")
}

# Run the calculation using the new 'Sub' object
BAI_Denali_Forest <- bai.out(Denali_Forest_Sub, diam = Denali_Forest_DBH_Clean)


#####################################################################
##############################################
# create data frame with tree IDs and diameters at coring height
is.numeric(Denali_Treeline_DBH$DBH_cm)
sapply(Denali_Treeline_DBH, class)
B <- transform(Denali_Treeline_DBH, DBH = as.numeric(DBH_cm))
sapply(B, class)
B[,2]<-B[,2]*10# # conversion in mm# # # cm to mm conversion
rownames(B)<-B[,1]# # assure the correct assignment# # 
head(B)
B
B<-B[colnames(All_DT ),]#using the column of AV.RWJF to order A
Denali_Treeline_DBH <- B[ ,c("ID","DBH_cm","DBH")]
years <- as.numeric(rownames(All_DT ))  # years in the tree-ring data

# Ensure your diameter data frame has NO missing values
Denali_Treeline_DTP <- Denali_Treeline_DBH %>% filter(!is.na(ID) & !is.na(DBH_cm))

#Force both sets of IDs to be character strings to avoid matching errors
tree_names <- as.character(colnames(All_DT ))
Denali_Treeline_DBH$ID <- as.character(Denali_Treeline_DBH$ID)

#  Create the clean diameter table by specifically looking for those tree_names
# This ensures the order is IDENTICAL to the columns in Denali_Treeline
Denali_Treeline_DBH_Clean <- Denali_Treeline_DTP[Denali_Treeline_DBH$ID %in% tree_names, ]

#  Sort it to match the column order exactly
Denali_Treeline_DBH_Clean <- Denali_Treeline_DBH_Clean[match(tree_names, Denali_Treeline_DBH_Clean$ID), ]

# Check the numbers before running
cat("Trees in Ring-Width File:", ncol(All_DT), "\n")
cat("Trees in Diameter Table:", nrow(Denali_Treeline_DBH_Clean), "\n")

# Find IDs that exist in BOTH the columns of Denali_Treeline AND the ID column of DTP
common_trees <- intersect(colnames(All_DT), Denali_Treeline_DBH$ID)

#  Subset the ring-width data to keep only those common trees
Denali_Treeline_Sub <- All_DT[, common_trees]

# Subset and reorder the diameter table to match exactly
Denali_Treeline_DBH_Clean <- Denali_Treeline_DBH[match(common_trees, Denali_Treeline_DBH$ID), ]

#  Final check: verify there are no NAs in the IDs
if(any(is.na(Denali_Treeline_DBH_Clean$ID))) {
  stop("There are still NAs in your ID column!")
}

#  Run the calculation using the new 'Sub' object
BAI_Denali_Treeline <- bai.out(Denali_Treeline_Sub, diam = Denali_Treeline_DBH_Clean)


##############################################
# create data frame with tree IDs and diameters at coring height
is.numeric(Brooks_Range_Forest_DBH$DBH_cm)
sapply(Brooks_Range_Forest_DBH, class)
B <- transform(Brooks_Range_Forest_DBH, DBH = as.numeric(DBH_cm))
sapply(B, class)
B[,2]<-B[,2]*10# # conversion in mm# # # cm to mm conversion
rownames(B)<-B[,1]# # assure the correct assignment# # 
head(B)
B
B<-B[colnames(All_BF ),]#using the column of AV.RWJF to order A
Brooks_Range_Forest_DBH <- B[ ,c("ID","DBH_cm","DBH")]
years <- as.numeric(rownames(All_BF ))  # years in the tree-ring data

# Ensure your diameter data frame has NO missing values
Brooks_Range_Forest_DTP <- Brooks_Range_Forest_DBH %>% filter(!is.na(ID) & !is.na(DBH_cm))

#1 Force both sets of IDs to be character strings to avoid matching errors
tree_names <- as.character(colnames(All_BF ))
Brooks_Range_Forest_DBH$ID <- as.character(Brooks_Range_Forest_DBH$ID)

#  Create the clean diameter table by specifically looking for those tree_names
# This ensures the order is IDENTICAL to the columns in Brooks_Range_Forest
Brooks_Range_Forest_DBH_Clean <- Brooks_Range_Forest_DTP[Brooks_Range_Forest_DBH$ID %in% tree_names, ]

# Sort it to match the column order exactly
Brooks_Range_Forest_DBH_Clean <- Brooks_Range_Forest_DBH_Clean[match(tree_names, Brooks_Range_Forest_DBH_Clean$ID), ]

#  Check the numbers before running
cat("Trees in Ring-Width File:", ncol(All_BF), "\n")
cat("Trees in Diameter Table:", nrow(Brooks_Range_Forest_DBH_Clean), "\n")

# Find IDs that exist in BOTH the columns of Brooks_Range_Forest AND the ID column of DTP
common_trees <- intersect(colnames(All_BF), Brooks_Range_Forest_DBH$ID)

#  Subset the ring-width data to keep only those common trees
Brooks_Range_Forest_Sub <- All_BF[, common_trees]

#  Subset and reorder the diameter table to match exactly
Brooks_Range_Forest_DBH_Clean <- Brooks_Range_Forest_DBH[match(common_trees, Brooks_Range_Forest_DBH$ID), ]

#  Final check: verify there are no NAs in the IDs
if(any(is.na(Brooks_Range_Forest_DBH_Clean$ID))) {
  stop("There are still NAs in your ID column!")
}

# Run the calculation using the new 'Sub' object
BAI_Brooks_Range_Forest <- bai.out(Brooks_Range_Forest_Sub, diam = Brooks_Range_Forest_DBH_Clean)

#####################################################################
##############################################
# create data frame with tree IDs and diameters at coring height
is.numeric(Brooks_Range_Treeline_DBH$DBH_cm)
sapply(Brooks_Range_Treeline_DBH, class)
B <- transform(Brooks_Range_Treeline_DBH, DBH = as.numeric(DBH_cm))
sapply(B, class)
B[,2]<-B[,2]*10# # conversion in mm# # # cm to mm conversion
rownames(B)<-B[,1]# # assure the correct assignment# # 
head(B)
B
B<-B[colnames(All_BT ),]#using the column of AV.RWJF to order A
Brooks_Range_Treeline_DBH <- B[ ,c("ID","DBH_cm","DBH")]
years <- as.numeric(rownames(All_BT ))  # years in the tree-ring data

# Ensure your diameter data frame has NO missing values
Brooks_Range_Treeline_DTP <- Brooks_Range_Treeline_DBH %>% filter(!is.na(ID) & !is.na(DBH_cm))

# Force both sets of IDs to be character strings to avoid matching errors
tree_names <- as.character(colnames(All_BT ))
Brooks_Range_Treeline_DBH$ID <- as.character(Brooks_Range_Treeline_DBH$ID)

# Create the clean diameter table by specifically looking for those tree_names
# This ensures the order is IDENTICAL to the columns in Brooks_Range_Treeline
Brooks_Range_Treeline_DBH_Clean <- Brooks_Range_Treeline_DTP[Brooks_Range_Treeline_DBH$ID %in% tree_names, ]

# Sort it to match the column order exactly
Brooks_Range_Treeline_DBH_Clean <- Brooks_Range_Treeline_DBH_Clean[match(tree_names, Brooks_Range_Treeline_DBH_Clean$ID), ]

#  Check the numbers before running
cat("Trees in Ring-Width File:", ncol(All_BT), "\n")
cat("Trees in Diameter Table:", nrow(Brooks_Range_Treeline_DBH_Clean), "\n")

#Find IDs that exist in BOTH the columns of Brooks_Range_Treeline AND the ID column of DTP
common_trees <- intersect(colnames(All_BT), Brooks_Range_Treeline_DBH$ID)

# Subset the ring-width data to keep only those common trees
Brooks_Range_Treeline_Sub <- All_BT[, common_trees]

# Subset and reorder the diameter table to match exactly
Brooks_Range_Treeline_DBH_Clean <- Brooks_Range_Treeline_DBH[match(common_trees, Brooks_Range_Treeline_DBH$ID), ]

# Final check: verify there are no NAs in the IDs
if(any(is.na(Brooks_Range_Treeline_DBH_Clean$ID))) {
  stop("There are still NAs in your ID column!")
}

# Run the calculation using the new 'Sub' object
BAI_Brooks_Range_Treeline <- bai.out(Brooks_Range_Treeline_Sub, diam = Brooks_Range_Treeline_DBH_Clean)

##########################################

write.csv(as.data.frame(BAI_Bluff_Fairbanks), "Bluff_Fairbanks_output.csv", row.names = TRUE)
write.csv(as.data.frame(BAI_Denali_Forest), "Denali_Forest_output.csv", row.names = TRUE)
write.csv(as.data.frame(BAI_Denali_Treeline), "Denali_Treeline_output.csv", row.names = TRUE)
write.csv(as.data.frame(BAI_Brooks_Range_Forest), "Brooks_Range_Forest_output.csv", row.names = TRUE)
write.csv(as.data.frame(BAI_Brooks_Range_Treeline), "Brooks_Range_Treeline_output.csv", row.names = TRUE)

#----------###SPLINE Detrend######--------------------
FBDET <- detrend(rwl =BAI_Bluff_Fairbanks, method = c("Spline"), nyrs = 30, f = 0.5, pos.slope = FALSE)
DFDET <- detrend(rwl =BAI_Denali_Forest, method = c("Spline"), nyrs = 30, f = 0.5, pos.slope = FALSE)
DTDET <- detrend(rwl =BAI_Denali_Treeline, method = c("Spline"), nyrs = 30, f = 0.5, pos.slope = FALSE)
BFDET <- detrend(rwl =BAI_Brooks_Range_Forest, method = c("Spline"), nyrs = 30, f = 0.5, pos.slope = FALSE)
BTDET <- detrend(rwl =BAI_Brooks_Range_Treeline, method = c("Spline"), nyrs = 30, f = 0.5, pos.slope = FALSE)


write.csv(as.data.frame(FBDET), "Bluff_Fairbanks_output1.csv", row.names = TRUE)
write.csv(as.data.frame(DFDET), "Denali_Forest_output1.csv", row.names = TRUE)
write.csv(as.data.frame(DTDET), "Denali_Treeline_output1.csv", row.names = TRUE)
write.csv(as.data.frame(BFDET), "Brooks_Range_Forest_output1.csv", row.names = TRUE)
write.csv(as.data.frame(BTDET), "Brooks_Range_Treeline_output1.csv", row.names = TRUE)
