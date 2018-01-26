# Results: DIC/WAIC indexes, figures and tables

## Results

This documentation reproduces the main results provided in Section 4 of the paper, including the computation of the DIC and WAIC indexes, the figures and the tables.

The results of posterior computation are available in the [`workspaces`](https://github.com/tommasorigon/India-SequentiaLogit/blob/master/workspaces) folder, and have been previously obtained in the [`estimation.md`](https://github.com/tommasorigon/India-SequentiaLogit/blob/master/estimation.md) document. We load all these files in the memory, along with the cleaned dataset obtained in the [`data-cleaning.md`](https://github.com/tommasorigon/India-SequentiaLogit/blob/master/data-cleaning.md) document, and the [`R core functions`](https://github.com/tommasorigon/India-SequentiaLogit/blob/master/core_functions.R). Again note that the `dataset` **is not made available** in this GitHub repository.



```r
library(dplyr)     # Required to manipulate the dataset
library(splines)   # Required for computing the B-spline basis
library(reshape2)  # Manipulating data
library(ggplot2)   # Graphical library
library(gridExtra)

# Load the results of the MCMC chain
rm(list=ls())

# Load the clean dataset
load("dataset.RData")

# Load the results of the MCMC chain
load("workspaces/ranef.RData")
load("workspaces/ranef_s.RData")
load("workspaces/dp_ranef.RData")
load("workspaces/dp_ranef_s.RData")

# Core functions. We will only need the function "IC"
source("core_functions.R")
```

## Information criteria (DIC and WAIC)

Following [Gelman et al. (2014)](https://link.springer.com/article/10.1007/s11222-013-9416-2) we compute the DIC (Deviance Information Criterion) and WAIC (Watanabe-Akaike information criterion) indexes. This can be done by using the `IC` function, included in the [`core_functions.R`](https://github.com/tommasorigon/India-SequentiaLogit/blob/master/core_functions.R) file. For each model, the IC function returns both the DIC, and the WAIC.

Notice that the factorization of our likelihood allows to decompose the DIC and WAIC indexes, that is

> DIC = DIC_Usage + DIC_Reversibility + DIC_Method

and similarly

> WAIC = WAIC_Usage + WAIC_Reversibility + WAIC_Method

The code to obtain Table 1 in the paper is provided below.


```r
tab <- rbind(IC(fit1_ranef) + IC(fit2_ranef) + IC(fit3_ranef),
IC(fit1_ranef_s) + IC(fit2_ranef_s) + IC(fit3_ranef_s),
IC(fit1_dp_ranef) + IC(fit2_dp_ranef) + IC(fit3_dp_ranef),
IC(fit1_dp_ranef_s) + IC(fit2_dp_ranef_s) + IC(fit3_dp_ranef_s))

tab[,2] <- -2*tab[,2]
tab <- tab[,c(1,2)]

rownames(tab) <- c("baseline","splines","mixture", "mixture + splines")
colnames(tab) <- c("DIC","-2*WAIC")
knitr::kable(round(tab,digits=2),format="markdown")
```



|                  |      DIC|  -2*WAIC|
|:-----------------|--------:|--------:|
|baseline          | 53507.70| 53503.41|
|splines           | 53092.90| 53088.62|
|mixture           | 53505.00| 53498.06|
|mixture + splines | 53091.25| 53083.61|


## Effects of the variables area, education, religion and child

In the following, we compute the posterior mean of the effects of the variables `area`, `education`, `religion` and `child`, together with the 0.95 credible intervals. These quantities are shown in Table 3 of the paper.


```r
alpha <- 0.05 # alpha = 1 - Coverage level

tab <-cbind(
  # Usage choice
  cbind(apply(fit1_dp_ranef_s$beta_Fix,2,mean),apply(fit1_dp_ranef_s$beta_Fix,2,function(x) quantile(x,alpha/2)),apply(fit1_dp_ranef_s$beta_Fix,2,function(x) quantile(x,1-alpha/2))),
  # Reversibility choice
  cbind(apply(fit2_dp_ranef_s$beta_Fix,2,mean),apply(fit2_dp_ranef_s$beta_Fix,2,function(x) quantile(x,alpha/2)),apply(fit2_dp_ranef_s$beta_Fix,2,function(x) quantile(x,1-alpha/2))),
  # Method choice
  cbind(apply(fit3_dp_ranef_s$beta_Fix,2,mean),apply(fit3_dp_ranef_s$beta_Fix,2,function(x) quantile(x,alpha/2)),apply(fit3_dp_ranef_s$beta_Fix,2,function(x) quantile(x,1-alpha/2)))
)

colnames(tab) <- c("Usage (mean)","Usage (lower)","Usage (upper)",
                   "Reversibility (mean)","Reversibility (lower)","Reversibility (upper)", 
                   "Method (mean)","Method (lower)", "Method (upper)")
rownames(tab) <- c("no child","one child","urban","muslim","christian","other","low","intermediate","high")

knitr::kable(round(tab,digits=2),format="markdown")
```



|             | Usage (mean)| Usage (lower)| Usage (upper)| Reversibility (mean)| Reversibility (lower)| Reversibility (upper)| Method (mean)| Method (lower)| Method (upper)|
|:------------|------------:|-------------:|-------------:|--------------------:|---------------------:|---------------------:|-------------:|--------------:|--------------:|
|no child     |        -3.71|         -3.88|         -3.54|                 2.20|                  1.69|                  2.75|         -0.19|          -0.61|           0.22|
|one child    |        -1.37|         -1.45|         -1.28|                 2.19|                  2.06|                  2.32|         -0.19|          -0.34|          -0.04|
|urban        |         0.24|          0.17|          0.31|                 0.29|                  0.21|                  0.38|          0.45|           0.32|           0.58|
|muslim       |        -0.43|         -0.52|         -0.34|                 1.24|                  1.12|                  1.36|          0.13|          -0.03|           0.29|
|christian    |        -0.26|         -0.48|         -0.03|                 0.00|                 -0.29|                  0.29|          0.36|          -0.13|           0.88|
|other        |         0.08|         -0.11|          0.29|                 0.46|                  0.25|                  0.66|          0.30|           0.02|           0.59|
|low          |         0.14|          0.05|          0.23|                 0.08|                 -0.03|                  0.19|          0.43|           0.25|           0.60|
|intermediate |         0.20|          0.12|          0.27|                 0.50|                  0.40|                  0.59|          0.71|           0.56|           0.86|
|high         |         0.27|          0.16|          0.37|                 1.28|                  1.14|                  1.41|          1.16|           0.97|           1.34|


## Age effects

The `age` effect is obtained by computing the posterior mean of each function f_k(), along with the 0.95 point-wise credible intervals. The first step consists in constructing the B-spline basis function evaluated on the grid of values 15,...,49. This is done by using the `spline.des` function from the `splines` R package. In the [`core_functions.R`](https://github.com/tommasorigon/India-SequentiaLogit/blob/master/core_functions.R) file, the construction of the `B` matrix proceeds exactly in the same way. 


```r
# Knots placement
inner_knots <- 40; degree <- 3
xl    <- min(dataset$age); xr <- max(dataset$age); dx <- (xr - xl) / (inner_knots-1)
knots <- seq(xl - degree * dx, xr + degree * dx, by = dx)

# Fixed quantities
B        <- spline.des(knots, 15:49, degree + 1, 0 * 15:49, outer.ok=TRUE)$design
```

The final number of columns of `B`, coinciding with the total number of splines coefficients, is


```r
ncol(B)
```

```
## [1] 42
```

Then, the Figure 3 can be obtained as follow


```r
# Posterior sample for each f_k()
eta1_spline <- t(B%*%t(fit1_dp_ranef_s$beta_spline))
eta2_spline <- t(B%*%t(fit2_dp_ranef_s$beta_spline))
eta3_spline <- t(B%*%t(fit3_dp_ranef_s$beta_spline))

# Plots
data.plot <- data.frame(x=15:49, y=colMeans(eta1_spline), 
                        ymax=apply(eta1_spline,2,function(x) quantile(x,0.975)),
                        ymin=apply(eta1_spline,2,function(x) quantile(x,1-0.975)),
                        k="Usage choice")
data.plot <- rbind(data.plot,data.frame(x=15:49,y=apply(eta2_spline,2,mean),
                                        ymax=apply(eta2_spline,2,function(x) quantile(x,0.975)),
                                        ymin=apply(eta2_spline,2,function(x) quantile(x,1-0.975)),
                                        k="Reversibility choice"))
data.plot <- rbind(data.plot,data.frame(x=15:49,y=apply(eta3_spline,2,mean),
                                        ymax=apply(eta3_spline,2,function(x) quantile(x,0.975)),
                                        ymin=apply(eta3_spline,2,function(x) quantile(x,1-0.975)),
                                        k="Method choice"))

p.spline <- ggplot(data = data.plot, aes(x = x, y = y,ymin=ymin,ymax=ymax)) + geom_line()  +theme_bw()+ xlab("") + ylab("") + geom_ribbon(alpha=0.25) + facet_grid(.~k)

ggsave("img/Age_effect.pdf",p.spline,device="pdf",width=8,height=2.6666)
ggsave("img/Age_effect.jpg",p.spline,device="jpg",width=8,height=2.6666)
```


![](https://github.com/tommasorigon/India-SequentiaLogit/blob/master/img/Age_effect.jpg)


## State-specific effects

In the following, we reproduce the boxplots for the state-specific effects (see Figure 4).


```r
# Relabel the names of the random effects appropriately
colnames(fit1_dp_ranef_s$beta_RF) <- colnames(fit2_dp_ranef_s$beta_RF) <- colnames(fit3_dp_ranef_s$beta_RF) <- levels(dataset$state)[-1]

# Usage choice
data.plot <- melt(as.matrix(fit1_dp_ranef_s$beta_RF))
data.plot$k <- "Usage choice"
# Ordering the levels according to the median
data.plot$Var2    <- factor(data.plot$Var2, levels = levels(data.plot$Var2)[order(apply(fit1_dp_ranef_s$beta_RF,2, median))])
p1 <- ggplot(data = data.plot, aes(x = Var2, y = value)) + geom_boxplot(outlier.alpha = 0) + theme_bw() + theme(axis.text.x = element_text(angle = 50,hjust = 1)) + theme(legend.position = "none")+ xlab("") + ylab("State effect") + facet_grid(~k) +ylim(c(-5,2.5))

# Reversibility choice
data.plot <- melt(as.matrix(fit2_dp_ranef_s$beta_RF))
data.plot$k <- "Reversibility choice"
# Ordering the levels according to the median
data.plot$Var2    <- factor(data.plot$Var2, levels = levels(data.plot$Var2)[order(apply(fit2_dp_ranef_s$beta_RF,2, median))])
p2 <- ggplot(data = data.plot, aes(x = Var2, y = value)) + geom_boxplot(outlier.alpha = 0) + theme_bw() + theme(axis.text.x = element_text(angle = 50,hjust = 1)) + theme(legend.position = "none")+ xlab("") + ylab("State effect") + facet_grid(~k) + ylim(c(-10,2.5))

# Method choice
data.plot <- melt(as.matrix(fit3_dp_ranef_s$beta_RF))
data.plot$k <- "Method choice"
# Ordering the levels according to the median
data.plot$Var2    <- factor(data.plot$Var2, levels = levels(data.plot$Var2)[order(apply(fit3_dp_ranef_s$beta_RF,2, median))])
p3 <- ggplot(data = data.plot, aes(x = Var2, y = value)) + geom_boxplot(outlier.alpha = 0) + theme_bw() + theme(axis.text.x = element_text(angle = 50,hjust = 1)) + theme(legend.position = "none")+ xlab("") + ylab("State effect") + facet_grid(~k) +ylim(-2,10)

ggsave("img/State.pdf",grid.arrange(p1,p2,p3,ncol=1),device="pdf",width=7.8,height=9)
ggsave("img/State.jpg",grid.arrange(p1,p2,p3,ncol=1),device="jpg",width=7.8,height=9)
```

![](https://github.com/tommasorigon/India-SequentiaLogit/blob/master/img/State.jpg)

## State-specific effects: clustering and maps

Clustering States will require additional libraries. We will make use of the ['mcclust'](https://cran.r-project.org/web/packages/mcclust/index.html) R package, which is based on the paper of [Fritsch and Ickstadt (2009)](https://projecteuclid.org/download/pdf_1/euclid.ba/1340370282). In particular, the `medv` function allows to reproduce the work of [Medvedovic et al. (2002)](https://www.cs.princeton.edu/~bee/courses/read/medvedovic-bioinformatics-2002.pdf).

The following code allows to obtain the cluster labels for the `Usage choice`, `Reversibility choice` and the `Method choice` models. Interesting clusters (mostly associated with the `Usage choice` model) are also reported.


```r
library(mcclust)

# First of all, we compute the dissimilarity matrix between States based on the MCMC labels
distance1 <- comp.psm(fit1_dp_ranef_s$S)
distance2 <- comp.psm(fit2_dp_ranef_s$S)
distance3 <- comp.psm(fit3_dp_ranef_s$S)

# Clustering states as in Medvedovic et al. (2002)
clust1_med  <- medv(distance1)
clust2_med  <- medv(distance2)
clust3_med  <- medv(distance3)

## Usage choice
levels(dataset$state)[-1][clust1_med==1]
```

```
##  [1] "Andhra Pradesh"     "Arunachal Pradesh"  "Bihar"             
##  [4] "Chandigarh"         "Chhattisgarh"       "Dadra+Nagar Haveli"
##  [7] "Daman & Diu"        "NCT of Delhi"       "Goa"               
## [10] "Gujarat"            "Haryana"            "Himachal Pradesh"  
## [13] "Jammu & Kashmir"    "Jharkhand"          "Karnataka"         
## [16] "Kerala"             "Madhya Pradesh"     "Maharashtra"       
## [19] "Orissa"             "Pondicherry"        "Punjab"            
## [22] "Rajasthan"          "Sikkim"             "Tamil Nadu"        
## [25] "Tripura"            "Uttarakhand"        "West Bengal"
```

```r
levels(dataset$state)[-1][clust1_med==2]
```

```
## [1] "Assam"     "Manipur"   "Meghalaya" "Mizoram"   "Nagaland"
```

The graphical representation of India requires additional external files, which **are not provided** in the GitHub repository, but can be dowloaded for free from the [Global Administrative Areas](http://www.gadm.org/) website. We downloaded the data for India in the `shapefile` format. Refer to the documentation of the [Global Administrative Areas](http://www.gadm.org/) for further information about the data source.

In order to make the map data compatible with ours, there are some extra steps to do. Having loaded the data in memory, we compare which States are present in the map but not included in the survey. In many cases, it is just a matter of relabeling (e.g. from `Dadra and Nagar Haveli` to `Dadra+Nagar Haveli`). However, some cases require special attention:

1. As mentioned, `Dadra and Nagar Haveli`, `Daman and Diu`, `Jammu and Kashmir`,`Odisha` and `Puducherry` can be relabelled to be coherent with the names of the survey.
2. For `Andaman and Nicobar` and `Lakshadweep`, both small islands, we do not have any observation in our dataset, and they will be ignored.
3. `Telengana` is grouped with `Andhra Pradesh`. In fact, `Telengana` became a separate State recently, but at the time of the interview it was part of `Andhra Pradesh`.


```r
library(rgeos)
library(maptools) # Required to load the shp data
library(mapproj)

# Load the data 
states.shp <- readShapeSpatial("IND_adm_shp/IND_adm1.shp")

# State names are converted to character
states.shp$NAME_1 <- as.character(states.shp$NAME_1)

# States that are present in the map file but are not present in the survey
states.shp$NAME_1[!(states.shp$NAME_1 %in% levels(dataset$state))]
```

```
## [1] "Andaman and Nicobar"    "Dadra and Nagar Haveli"
## [3] "Daman and Diu"          "Jammu and Kashmir"     
## [5] "Lakshadweep"            "Odisha"                
## [7] "Puducherry"             "Telangana"
```

```r
# Some states have just a different name
states.shp$NAME_1[1]
```

```
## [1] "Andaman and Nicobar"
```

```r
states.shp$NAME_1[8]  <- "Dadra+Nagar Haveli"
states.shp$NAME_1[9]  <- "Daman & Diu"
states.shp$NAME_1[14] <- "Jammu & Kashmir"
states.shp$NAME_1[18]
```

```
## [1] "Lakshadweep"
```

```r
states.shp$NAME_1[26] <- "Orissa"
states.shp$NAME_1[27] <- "Pondicherry"

# Telengana has to be merged with Andhra Pradesh
states.shp$NAME_1[33] <- "Andhra Pradesh"; states.shp$ID_1[33] <- 2
states.shp$NAME_1 <- as.factor(states.shp$NAME_1)

# fortify shape file to get into dataframe
states.shp.f <- fortify(states.shp, region = "ID_1")
```

The maps in Figure 5 are obtained as follow, leveraging the posterior median of each state-specific effect in the three different models. The baseline State `Uttar Pradesh` was fixed equal to 0.


```r
# Usage choice
state_effects <- c(0,apply(fit1_dp_ranef_s$beta_RF,2,median)) 
id         <- states.shp$NAME_1 %in% levels(dataset$state)
data.plot  <- data.frame(NAME_1 = states.shp$NAME_1[id], id = states.shp$ID_1[id])
data.plot <- merge(data.plot,data.frame(NAME_1=levels(dataset$state),Effect=state_effects),by="NAME_1")

# Merge with coefficients and Reorder
merge.shp.coef <- merge(states.shp.f, data.plot, by = "id", all.x = TRUE)
final.plot     <- merge.shp.coef[order(merge.shp.coef$order), ]

# Final plot - Attention, It could be very slow!
final.plot$k <- "Usage Choice"
p4 <- ggplot() + geom_polygon(data = final.plot, aes(x = long, y = lat, group = group, fill=Effect),col = "black") + coord_map()   +xlab("Longitude") + ylab("Latitude") + facet_grid(~k) + theme_bw() + theme(legend.position="none") + scale_fill_gradient(low = "grey28",high = "white", breaks=-6:6)

# Reversibility choice
state_effects <- c(0,apply(fit2_dp_ranef_s$beta_RF,2,median)) 
id         <- states.shp$NAME_1 %in% levels(dataset$state)
data.plot  <- data.frame(NAME_1 = states.shp$NAME_1[id], id = states.shp$ID_1[id])
data.plot <- merge(data.plot,data.frame(NAME_1=levels(dataset$state),Effect=state_effects),by="NAME_1")

# Merge with coefficients and Reorder
merge.shp.coef <- merge(states.shp.f, data.plot, by = "id", all.x = TRUE)
final.plot     <- merge.shp.coef[order(merge.shp.coef$order), ]

# Final plot
final.plot$k <- "Reversibility Choice"
p5 <- ggplot() + geom_polygon(data = final.plot, aes(x = long, y = lat, group = group, fill=Effect),col = "black") + coord_map()   +xlab("Longitude") + ylab("Latitude") + facet_grid(~k) + theme_bw() + theme(legend.position="none") + scale_fill_gradient(low = "grey28",high = "white", breaks=-6:6)

# Method choice 
state_effects <- c(0,apply(fit3_dp_ranef_s$beta_RF,2,median))
id         <- states.shp$NAME_1 %in% levels(dataset$state)
data.plot  <- data.frame(NAME_1 = states.shp$NAME_1[id], id = states.shp$ID_1[id])
data.plot <- merge(data.plot,data.frame(NAME_1=levels(dataset$state),Effect=state_effects),by="NAME_1")

# Merge with coefficients and Reorder
merge.shp.coef <- merge(states.shp.f, data.plot, by = "id", all.x = TRUE)
final.plot     <- merge.shp.coef[order(merge.shp.coef$order), ]

# Final plot
final.plot$k <- "Method Choice"
p6 <- ggplot() + geom_polygon(data = final.plot, aes(x = long, y = lat, group = group, fill=Effect),col = "black") + coord_map()   +xlab("Longitude") + ylab("Latitude") + facet_grid(~k) + theme_bw() + theme(legend.position="none") + scale_fill_gradient(low = "grey28",high = "white", breaks=-6:6)

# Commented since it could be very slow!
ggsave("img/map.pdf",grid.arrange(p4,p5,p6,ncol=3),device="pdf",width=7.5,height=3.375)
ggsave("img/map.jpg",grid.arrange(p4,p5,p6,ncol=3),device="jpg",width=7.5,height=3.375)
```

![](https://github.com/tommasorigon/India-SequentiaLogit/blob/master/img/map.jpg)

