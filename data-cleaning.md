# Preliminary operations for the IHDS-II dataset

## The Indian Human Development Survey-II (IHDS-II)

This short document explains in detail all the preliminary operations performed to the [IDHS-II](http://ihds.info/IHDS-II) dataset in our paper. From the [official documentation](http://www.icpsr.umich.edu/icpsrweb/content/DSDR/idhs-II-data-guide.html):

> The India Human Development Survey-II (IHDS-II), 2011-12 is a nationally representative, multi-topic survey of 42,152 households in 1,420 villages and 1,042 urban neighborhoods across India. These data are mostly re-interviews of households interviewed for IHDS (ICPSR 22626) in 2004-05.

We cannot redistribute the data, but they can be downloaded from the [Data Sharing for Demographic Research Archive](http://www.icpsr.umich.edu/icpsrweb/ICPSR/studies/36151) at ICPSR. Download will require a registration but is completely free. We focus on the `R` version of the dataset called **DS3: Eligible Women**. The download provides a zip directory `ICPSR_36151.zip` which contains data and additional documentation. The `R` dataset `36151-0003-Data.rda` of interest is in the sub-directory `DS0003`.

Eligible women are ever-married women aged 15 - 49. Non-eligible women will be excluded from the dataset.


```r
rm(list = ls())
# The DS3 dataset is located in the following folder
load("ICPSR_36151/DS0003/36151-0003-Data.rda")

# The dataset have the following dimensions
dim(da36151.0003)
```

```
## [1] 39523   580
```

## Selection of the variables of interest

Among all the 580 variables we select only those of interest for our analysis (or relevant for data cleaning operations). Refer to Section 1.1 of the paper for a discussion on the response, and the covariates of interest. The variables of interest are coded according a [codebook](http://www.icpsr.umich.edu/cgi-bin/file?comp=none&study=36151&ds=3&file_id=1207405&path=ICPSR). The exact questions posed to each eligible woman are instead described in the [questionnaire](http://www.icpsr.umich.edu/cgi-bin/file?comp=none&study=36151&ds=3&file_id=1212084&path=ICPSR). We considered the following variables:

- `STATEID`. Relabeled as `state`, represents the name of the State where each woman lives.
- `EW6`. Relabeled as `age`, represents the age of each eligible woman in 2011.
- `EW8`. Relabeled as `education`, represents the year of education completed.
- `EW9`. Relabeled as `child`, represents the number of children alive.
- `ID11`. Relabeled as `religion`, represents the religion of the head of the household.
- `URBAN2011`. Relabeled as `area`, represents the urban residence from census 2011.
- `FP1`.  Relabeled as `pregnant`, defines whether the woman is currently pregnant or not.
- `FP2A`. Relabeled as `contraceptive` indicates whether the woman is using a contraceptive method, provided that she is not pregnant.
- `FP2B`. Relabeled as `method`, represents the main contraceptive method adopted, provided that the woman is using contraceptives.
- `EWELIGIBLE`. Relabeled as `eligible`, defines whether the woman is eligible or not (women 15-49, ever married).



```r
library(dplyr)

# Conversion of the dataframe into a more manageable class
IHDS_II <- tbl_df(da36151.0003); rm(da36151.0003)

# Selection of relevant variables
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
```

## Selection of eligible women

We select only eligible women (i.e. ever married women and aged 15 - 49). Notice that among the eligible women there is **one woman** who declared to be 81. We discarded this woman as well, since it seems to be either a transcription error or another kind of gross error.


```r
dataset <- dataset %>% filter(eligible=="(1) Yes 1")
dataset <- dataset %>% filter(age <= 49)
```

Eligible women comprise a total of `35281` observations.

## Cleaning missing or not coherent values 

#### 1. Pregnancy and contraceptive usage

Some of the selected variables present missing values. In certain circumstances these are induced by the structure of the questionnaire itself. For instance, a pregnant woman is not asked to declare whether she is using or not a contraceptive method. 


```r
table(dataset$pregnant, dataset$contraceptive, useNA="always")
```

```
##            
##             (0) No 0 (1) Yes 1  <NA>
##   (0) No 0      8066     23522   242
##   (1) Yes 1        0         0  1700
##   <NA>            28         5  1718
```

We decided to exclude from the analysis pregnant women (`1700`), and those women who do not declare their pregnancy status (`1718 + 28 + 5`), comprising a total of `1700 + 1718 + 28 + 5 = 3451` cases. The `28 + 5 = 33` missing cases reported in the table above, are women who declared to be `unsure` about their pregnancy status. We held them out from the dataset as well.


```r
dataset <- filter(dataset, pregnant == "(0) No 0")
```

#### 2. Contraceptive usage and method usage

Similarly, the `method` variable has some "structural" missingness, due to the fact that only women using contraceptives are asked about the contraceptive `method`. Hence, women having missing values in the `method` variable, but non-missing values in the `contraceptive` variable, should be kept as cases of women using no contraceptive methods. Instead, women having missing values both in the `method` and `contraceptive` variables, should be removed, since we do not have information on contraceptive behavior for these women. To perform this operation, we first create an additional category in the `method` variable (called `1. No contraceptive method`) to distinguish between these two different cases of missingness.


```r
dataset$method <- as.character(dataset$method)
dataset$method[dataset$contraceptive == "(0) No 0"] <- "1. No contraceptive method"
```

Then, consistent with the above discussion, we hold out only the real cases of missingness (879 observations), as well as those women who declare to use "others" contraceptive methods (417 observations).


```r
# Filtering missing values and "Others"
dataset <- filter(dataset, method!="(11) Others 11")
```

#### 3. Other missingness

There are still some missing values in the covariates, but their number is negligible (3 and 7), and therefore we can safely exclude them from the analysis. 


```r
# Summary of missing values
summary(subset(dataset,select=c(education, child)))
```

```
##              education         child       
##  (00) none 0      :11084   Min.   : 0.000  
##  (10) Secondary 10: 2961   1st Qu.: 2.000  
##  (05) 5th class 5 : 2771   Median : 2.000  
##  (08) 8th class 8 : 2389   Mean   : 2.517  
##  (09) 9th class 9 : 2235   3rd Qu.: 3.000  
##  (Other)          : 9091   Max.   :13.000  
##  NA's             :    3   NA's   :7
```

```r
dataset <- filter(dataset, !is.na(contraceptive), !is.na(education),  !is.na(child))
```

The final dataset comprises a total of `30524` observations out of the original `35281` elegible women.

## Relabeling variables and categories

In this section, we slightly modify the name of some categories for aesthetical reasons. This is accomplished via the function `str_sub()` which deletes unnecessary numbers and ids in the name of some categories. We also fix, for each categorical variable, its baseline consistent with those considered in the paper. Note also that, for the `state` variable the category `Delhi` is relabelled `NCT of Delhi` and `Uttar Pradesh` is fixed as baseline category.


```r
# This part "modifies" the name of some variables - just for aesthetical reasons.
library(stringr)
levels(dataset$state)     <- str_sub(levels(dataset$state), 6, -4)
levels(dataset$education) <- str_sub(levels(dataset$education), 6, -2)
levels(dataset$area)      <- str_sub(levels(dataset$area), 5, -3)
levels(dataset$contraceptive) <- str_sub(levels(dataset$contraceptive), 5, -3)
levels(dataset$method)    <- str_sub(levels(dataset$method), 6, -3)
levels(dataset$religion)  <- str_sub(levels(dataset$religion), 5, -3)

# Re-order alphabetically the categories
dataset$state <- factor(as.character(dataset$state))
# State categories
levels(dataset$state)[levels(dataset$state)=="Delhi"] <- "NCT of Delhi"

# Uttar pradesh is now the baseline category.
new_levels    <- c("Uttar Pradesh",levels(dataset$state)[-which(levels(dataset$state)=="Uttar Pradesh")])
dataset$state <- factor(dataset$state,levels=new_levels); rm(new_levels)
levels(dataset$state)
```

```
##  [1] "Uttar Pradesh"      "Andhra Pradesh"     "Arunachal Pradesh" 
##  [4] "Assam"              "Bihar"              "Chandigarh"        
##  [7] "Chhattisgarh"       "Dadra+Nagar Haveli" "Daman & Diu"       
## [10] "NCT of Delhi"       "Goa"                "Gujarat"           
## [13] "Haryana"            "Himachal Pradesh"   "Jammu & Kashmir"   
## [16] "Jharkhand"          "Karnataka"          "Kerala"            
## [19] "Madhya Pradesh"     "Maharashtra"        "Manipur"           
## [22] "Meghalaya"          "Mizoram"            "Nagaland"          
## [25] "Orissa"             "Pondicherry"        "Punjab"            
## [28] "Rajasthan"          "Sikkim"             "Tamil Nadu"        
## [31] "Tripura"            "Uttarakhand"        "West Bengal"
```

The qualitative variables `education`, `religion`, and `child`, are measured on several categories. To maintain simple interpretation, and consistent with other studies on contraceptive preferences in India, we aggregate some of these categories. In particular:

- `education`. We create a four level factor for education (No education, Low, Intermediate, and High).
- `religion`. We consider the most frequent religions in India (Hindu, Muslim, Christian), and we grouped the other religions in a single category: Others.
- `child`. We consider a three level factor for the number of children: none, one, or more than one. The `More than 1` category is fixed as baseline.

The `R` code to perform these operations is:

```r
# Grouping categories of education
levels(dataset$education) <- c("None", rep("Low", 5), rep("Intermediate", 6), rep("High", 5))

# Grouping categories of religion
levels(dataset$religion) <- c("Hindu", "Muslim", "Christian", rep("Other", 6))

# Relabeling numbers of children
dataset$child[dataset$child > 1] <- "More than 1"
dataset$child <- as.factor(dataset$child)
dataset$child <- factor(dataset$child,levels=levels(dataset$child)[c(3,1,2)])
```

The variable `method` is the response variable, and is composed by the following categories:

- `1. No contraceptive method`
- `2. Sterilization`. Includes female sterilization, hysterectomy and male sterilization.
- `3. Natural methods`. Includes periodic abstinence and withdrawal.
- `4. Modern methods`. Includes condom, copper T/IUD, diaphgram/jelly, injectible contraception and oral pill.

See Section 1.1 of the paper for a justification on the definition of these four categories. The `R` code to create the final response variable is:

```r
# Converting variable into a factor
dataset$method <- factor(as.character(dataset$method))

# Relabeling method categories
levels(dataset$method)[c(6,7,10)] <- "2. Sterilization"   
levels(dataset$method)[c(7,8)] <- "3. Natural methods"
levels(dataset$method)[c(1,2,3,4,5)] <- "4. Modern methods"

# Reorder the categories
dataset$method <- factor(as.character(dataset$method))
```

The final response variable is composed as follows


```r
knitr::kable(t(as.matrix(table(dataset$method))),format='markdown')
```



| 1. No contraceptive method| 2. Sterilization| 3. Natural methods| 4. Modern methods|
|--------------------------:|----------------:|------------------:|-----------------:|
|                       8059|            15323|               2914|              4228|

## Final dataset

The final dataset is saved in the workspace and summarized below.


```r
# Save the dataset into a workspace.
save.image("dataset.RData")
# A short description
str(dataset)
```

```
## Classes 'tbl_df', 'tbl' and 'data.frame':	30524 obs. of  10 variables:
##  $ state        : Factor w/ 33 levels "Uttar Pradesh",..: 15 15 15 15 15 15 15 15 15 15 ...
##  $ age          : num  49 26 33 43 47 38 25 42 37 25 ...
##  $ education    : Factor w/ 4 levels "None","Low","Intermediate",..: 1 3 4 1 1 1 1 1 1 1 ...
##  $ child        : Factor w/ 3 levels "More than 1",..: 1 1 1 1 1 1 3 1 1 1 ...
##  $ religion     : Factor w/ 4 levels "Hindu","Muslim",..: 2 2 2 2 2 2 2 2 2 2 ...
##  $ area         : Factor w/ 2 levels "rural","urban": 1 1 1 1 1 1 1 1 1 1 ...
##  $ pregnant     : Factor w/ 2 levels "(0) No 0","(1) Yes 1": 1 1 1 1 1 1 1 1 1 1 ...
##  $ contraceptive: Factor w/ 2 levels "No","Yes": 2 2 2 2 2 2 2 2 2 1 ...
##  $ method       : Factor w/ 4 levels "1. No contraceptive method",..: 3 4 2 2 3 2 4 3 3 1 ...
##  $ eligible     : Factor w/ 2 levels "(0) No 0","(1) Yes 1": 2 2 2 2 2 2 2 2 2 2 ...
```


## Descriptive analysis

The following commands reproduce Figure 2 of our paper,  describing an empirical relation between the variable `age` and the probability of the different contraceptive choices.


```r
library(dummies)
library(reshape2)
library(ggplot2)

data.plot <- aggregate(dummy(dataset$method)  ~ age, mean, data = dataset)[,-1]
colnames(data.plot) <- levels(dataset$method)
data.plot <- cbind(Age=15:49,melt(data.plot))
p0 <- ggplot(data = data.plot, aes(x = Age, y = value)) + geom_point() + ylab("") + xlab("") + theme_bw() + facet_grid(~variable)

ggsave("img/frequencies.pdf",p0,device="pdf",width=8,height=2.666)
ggsave("img/frequencies.jpg",p0,device="jpg",width=8,height=2.666)
```

![](https://raw.githubusercontent.com/tommasorigon/India-SequentiaLogit/master/img/frequencies.jpg)


