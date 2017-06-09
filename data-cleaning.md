# Preliminary operations for the IHDS-II dataset
Tommaso Rigon  



## The Indian Human Development Survey-II (IHDS-II)

This short tutorial explain in detail all the preliminary operations performed to the [IDHS-II](http://ihds.info/IHDS-II) dataset. The code used is also available [here](https://github.com/tommasorigon/India-SequentiaLogit/blob/master/data-cleaning.R). From the [official documentation](http://www.icpsr.umich.edu/icpsrweb/content/DSDR/idhs-II-data-guide.html):

> The India Human Development Survey-II (IHDS-II), 2011-12 is a nationally representative, multi-topic survey of 42,152 households in 1,420 villages and 1,042 urban neighborhoods across India. These data are mostly re-interviews of households interviewed for IHDS (ICPSR 22626) in 2004-05.

We cannot re-distribute the data they can be downloaded from the [Data Sharing for Demographic Research Archive](http://www.icpsr.umich.edu/icpsrweb/ICPSR/studies/36151) at ICPSR. Download will require a registration but is free. 

We will use the dataset called **DS3: Eligible Women**. Eligible women are ever-married women aged 15 - 49. Those ever-married women older than 49 years that were interviewed in the initial IHDS wave, **are included in the dataset even if they are not eligible anymore**.


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

## Selection of the variable of interest

Among all the $580$ variables we select only those which we will include in our analysis. The variables of interest are coded according a [codebook](http://www.icpsr.umich.edu/cgi-bin/file?comp=none&study=36151&ds=3&file_id=1207405&path=ICPSR). The exact question posed to each eligible woman is instead described in the [questionnaire](http://www.icpsr.umich.edu/cgi-bin/file?comp=none&study=36151&ds=3&file_id=1212084&path=ICPSR). We considered the following variables

- `STATEID`. Relabeled as `state`, contains the name of the State.
- `EW6`. Relabeled as `age`, represents the age of each eligible woman in 2011.
- `EW8`. Relabeled as `education`, represents the year of education completed (illiterate=00,5th class=05, bachelors=15, above bachelors=16).
- `EW9`. Relabeled as `child`, represents the number of children alive.
- `ID11`. Relabeled as `religion`, represents the religion of the head of the household.
- `URBAN2011`. Relabeled as `area`, represents the urban residence from census 2011.
- `FP1`.  Relabeled as `pregnant`, represents whether the woman is currently pregnant or not.
- `FP2A`. Relabeled as `contraceptive` indicates whether a woman is using a contraceptive method, provided that she is not pregnant.
- `FP2B`. Relabeled as `method`, represents the main contraceptive method adopted, provided that she is using one contraception.
- `EWELIGIBLE`. Relabeled as `eligible`, represents the eligible women in the dataset (women 15-49, ever married).



```r
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
```

## Selecting eligible women

First of all, we select only eligible women, that is, ever married women and aged 15 - 49. Notice that among eligible women there is **one woman** who declared to be 81; we discarded that women as well, since it seems to be either a transcription error or another kind of gross error.


```r
dataset <- dataset %>% filter(eligible=="(1) Yes 1")
dataset <- dataset %>% filter(age <= 49)
```

Eligible units comprises a total of 35281 observations.

## Missing values 

#### 1. Pregnancy and contraceptive usage
Some of the selected variables present missing values. In certain circumstances these are structural: e.g. a pregnant woman is not asked to declare whether she is using or not a contraceptive method. 


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

We decided to exclude from the analysis pregnant women (`1700`) and those women who do not declared their pregnancy status (`1718 + 28 + 5`), comprising a total of `1700 + 1718 + 28 + 5 = 4185` cases. The `28 + 5 = 33` missing cases ---reported in the table above--- are women who declared to be `unsure` about their pregnancy status; we holded them out from the dataset as well.


```r
# Notice that the following command will exclude also missing values.
dataset <- filter(dataset, pregnant == "(0) No 0")
```

#### 2. Contraceptive usage and method usage

Similarly, the `method` variable present some "structural" missingness. Women not using any kind of contraception were not aked about the contraception method. Therefore, we create an additional level called `1. No contraceptive method`.


```r
dataset$method <- as.character(dataset$method)
dataset$method[dataset$contraceptive == "(0) No 0"] <- "1. No contraceptive method"
```

Then, we excluded from the analysis the remaining missing values (879 observations), as well as those women who declare to use "others" contraceptive methods (417 observations).


```r
# Filtering missing values and "Others"
dataset <- filter(dataset, method!="(11) Others 11")
```


#### 3. Other missingness

Unfortunately, there are still some missing values among the some variables but their number is negligible (3 and 7), and therefore we can safely exclude them from the analysis. 


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

The dataset now comprises a total of 30524 observations out of `35281`, which is the number of eligible units.

## Relabeling variables and categories

For aesthetical reasons, we changed the name of some categories. The category `Delhi` is relabeled `NCT of Delhi` and `Uttar Pradesh` is setted as baseline level.


```r
library(stringr)
levels(dataset$state)     <- str_sub(levels(dataset$state), 6, -4)
levels(dataset$education) <- str_sub(levels(dataset$education), 6, -2)
levels(dataset$area)      <- str_sub(levels(dataset$area), 5, -3)
levels(dataset$contraceptive) <- str_sub(levels(dataset$contraceptive), 5, -3)
levels(dataset$method)    <- str_sub(levels(dataset$method), 6, -3)
levels(dataset$religion)  <- str_sub(levels(dataset$religion), 5, -3)

# Re-order alphabetically the levels
dataset$state <- factor(as.character(dataset$state))
# State levels
levels(dataset$state)[levels(dataset$state)=="Delhi"] <- "NCT of Delhi"

# Uttar pradesh is now the baseline level.
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

Then, we grouped some of the covariates' labels as follow:

- `education`. We create a four level factor: `No education`, `Low` , `Intermediate` and `High`.
- `religion`. We considered the most frequent religions in India (`Hindu`, `Muslim`, `Christian`), and we gropued the other religions together.
- `child`. We considered a three level factor for the number of child: `None`, `One` or `More than one`. The \texttt{More than 1} category is setted as baseline.



```r
# Grouping levels of education
levels(dataset$education) <- c("None", rep("Low", 5), rep("Intermediate", 6), rep("High", 5))

# Grouping levels of religion
levels(dataset$religion) <- c("Hindu", "Muslim", "Christian", rep("Other", 6))

# Relabeling numbers of child
dataset$child[dataset$child > 1] <- "More than 1"
dataset$child <- as.factor(dataset$child)
dataset$child <- factor(dataset$child,levels=levels(dataset$child)[c(3,1,2)])
```

The variable `method` is the response variable, which is composed by the following levels:

- `1. No contraceptive method`
- `2. Sterilization`. Includes female sterilization, hysteroctomy and male sterilization.
- `3. Natural methods`. Includes periodic abstinence and withdrawal
- `4. Modern methods`. Includes condom, copper T/IUD, diaphgram / jelly, injectible contraception and oral pill


```r
# Converting variable into a factor
dataset$method <- factor(as.character(dataset$method))

# Relabeling method levels
levels(dataset$method)[c(6,7,10)] <- "2. Sterilization"   
levels(dataset$method)[c(7,8)] <- "3. Natural methods"
levels(dataset$method)[c(1,2,3,4,5)] <- "4. Modern methods"

# Reorder the levels
dataset$method <- factor(as.character(dataset$method))

knitr::kable(t(as.matrix(table(dataset$method))))
```

 1. No contraceptive method   2. Sterilization   3. Natural methods   4. Modern methods
---------------------------  -----------------  -------------------  ------------------
                       8059              15323                 2914                4228

## Final dataset


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

