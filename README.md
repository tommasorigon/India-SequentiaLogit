## Bayesian semiparametric modelling of contraceptive behaviour in India via sequential logistic regressions

This repository explain in detail all operations described in the paper [Rigon, Durante and Torelli (2016). *Bayesian semiparametric modelling of contraceptive behavior in India via sequential logistic regressions* [https://arxiv.org/abs/1405.7555]](https://arxiv.org/abs/1405.7555). The aim of this document is to fully reproduce the results of our paper and to provide additional tools (e.g. a Shiny application), which could help in understanding the results of our model.

In doing so, we organized this additional documentation in four sections, which should be executed in order. The files are:

- A [`data-cleaning.md`](https://github.com/tommasorigon/India-SequentiaLogit/blob/master/data-cleaning.md) document, in which all the preliminary operations are discussed, including the original source of the data.
- A [`estimation.md`](https://github.com/tommasorigon/India-SequentiaLogit/blob/master/estimation.md) document, in which the model discussed in Section 3 of our paper is fitted using MCMC.
- A [`results.md`](https://github.com/tommasorigon/India-SequentiaLogit/blob/master/results.md) document, in which graph and tables of Section 4 of our paper are reproduced.
- A [`predictive_pref.md`](https://github.com/tommasorigon/India-SequentiaLogit/blob/master/predictive_perf.md) document, in which out-of-sample performances are computed for our model and for competing models. 

## Shiny web-application

We builded a [Shiny app](https://github.com/tommasorigon/India-SequentiaLogit/tree/master/SequentiaLogisticApp) in which are displayed the posterior probabilities for the contraceptive usage of a given woman. The characteristic of this woman (e.g. `age`, `child`) can be defined interactively by the user.

The application relies on three `R` packages: `shiny`, `shinydashboard` and `ggplot2`, which should be installed in the local environment, if not already present. Then, the Shiny application can be launched by executing the following `R` code.

```r
# install.packages(shiny)
# install.packages(shinydashboard)
# install.packages(ggplot2)

runGitHub("India-SequentiaLogit","tommasorigon", subdir = "SequentiaLogisticApp")
```
