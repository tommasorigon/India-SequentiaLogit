## ----loading the data----------------------------------------------------
rm(list = ls())
# The DS3 dataset is located in the following folder
load("ICPSR_36151/DS0003/36151-0003-Data.rda")

# The dataset have the following dimensions
dim(da36151.0003)

## ----trasforming data,message=FALSE--------------------------------------
library(dplyr)

# Conversion of the dataframe into a more manageable class
IHDS_II <- tbl_df(da36151.0003); rm(da36151.0003)

# Selection of relevant questions
dataset <- IHDS_II %>% transmute(state = STATEID, 
                                 age = EW6, 
                                 education = EW8, 
                                 child = EW9,  
                                 religion = ID11, 
                                 area = URBAN2011, 
                                 pregnant = FP1, 
                                 contraceptive = FP2A, 
                                 method = FP2B, 
                                 eligible= EWELIGIBLE)

## ----filtering eligible women--------------------------------------------
dataset <- dataset %>% filter(eligible=="(1) Yes 1")
dataset <- dataset %>% filter(age <= 49)

## ----table pregnant contraceptive----------------------------------------
table(dataset$pregnant, dataset$contraceptive, useNA="always")

## ----no pregnant---------------------------------------------------------
# Notice that the following command will exclude also missing values.
dataset <- filter(dataset, pregnant == "(0) No 0")

## ----method cleaning-----------------------------------------------------
dataset$method <- as.character(dataset$method)
dataset$method[dataset$contraceptive == "(0) No 0"] <- "1. No contraceptive method"

## ----filtering methods---------------------------------------------------
# Filtering missing values and "Others"
dataset <- filter(dataset, method!="(11) Others 11")

## ----filtering other NA--------------------------------------------------
# Summary of missing values
summary(subset(dataset,select=c(education, child)))
dataset <- filter(dataset, !is.na(contraceptive), !is.na(education),  !is.na(child))

## ----cleaning levels-----------------------------------------------------
library(stringr)
levels(dataset$state)     <- str_sub(levels(dataset$state), 6, -4)
levels(dataset$education) <- str_sub(levels(dataset$education), 6, -2)
levels(dataset$area)      <- str_sub(levels(dataset$area), 5, -3)
levels(dataset$contraceptive) <- str_sub(levels(dataset$contraceptive), 5, -3)
levels(dataset$method)    <- str_sub(levels(dataset$method), 6, -3)
levels(dataset$religion)  <- str_sub(levels(dataset$religion), 5, -3)

# Reorder alphabetically the levels
dataset$state <- factor(as.character(dataset$state))
# State levels
levels(dataset$state)[levels(dataset$state)=="Delhi"] <- "NCT of Delhi"

# Uttar pradesh is now the baseline level.
new_levels    <- c("Uttar Pradesh",levels(dataset$state)[-which(levels(dataset$state)=="Uttar Pradesh")])
dataset$state <- factor(dataset$state,levels=new_levels); rm(new_levels)
levels(dataset$state)

## ----grouping levels-----------------------------------------------------
# Grouping levels of education
levels(dataset$education) <- c("None", rep("Low", 5), rep("Intermediate", 6), rep("High", 5))

# Grouping levels of religion
levels(dataset$religion) <- c("Hindu", "Muslim", "Christian", rep("Other", 6))

# Relabeling numbers of child
dataset$child[dataset$child > 1] <- "More than 1"
dataset$child <- as.factor(dataset$child)
dataset$child <- factor(dataset$child,levels=levels(dataset$child)[c(3,1,2)])

## ----relabeling method---------------------------------------------------
# Converting variable into a factor
dataset$method <- factor(as.character(dataset$method))

# Relabeling method levels
levels(dataset$method)[c(6,7,10)] <- "2. Sterilization"   
levels(dataset$method)[c(7,8)] <- "3. Natural methods"
levels(dataset$method)[c(1,2,3,4,5)] <- "4. Modern methods"

# Reorder the levels
dataset$method <- factor(as.character(dataset$method))

knitr::kable(t(as.matrix(table(dataset$method))))

## ----saving output-------------------------------------------------------
# Save the dataset into a workspace.
save.image("dataset.RData")
# A short description
str(dataset)
