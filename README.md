## Bayesian semiparametric modelling of contraceptive behaviour in India via sequential logistic regressions

This repository is associated with the paper [Rigon, Durante and Torelli (2016). *Bayesian semiparametric modelling of contraceptive behavior in India via sequential logistic regressions*](https://arxiv.org/abs/1405.7555), and aims at providing detailed materials to fully reproduce the analyses in the paper. An interactive Shiny application is also made available to interactively visualize relevant quantities of direct interest for current family planning policies in India.

The documentation is organized in four sections described below.  

- [`data-cleaning.md`](https://github.com/tommasorigon/India-SequentiaLogit/blob/master/data-cleaning.md): contains detailed guidelines and code to download the raw data, and perform preliminary operations to clean the dataset.
- A [`estimation.md`](https://github.com/tommasorigon/India-SequentiaLogit/blob/master/estimation.md): contains detailed guidelines and code to perform posterior computation for the model presented in Section 2 of the paper, and the associated sub-models described in Section 4. Some convergence diagnostics for the MCMC chains are also discussed.
- A [`results.md`](https://github.com/tommasorigon/India-SequentiaLogit/blob/master/results.md): contains detailed guidelines and code to reproduce the graphs and tables in Section 4, where posterior inference under the proposed model is performed.
- A [`predictive_pref.md`](https://github.com/tommasorigon/India-SequentiaLogit/blob/master/predictive_perf.md): contains detailed guidelines and code to study out-of-sample predictive performance, and compare results with relevant competitors.

These sections should be executed in order.

In the repository we also made available the file [`core_functions.R`](https://github.com/tommasorigon/India-SequentiaLogit/blob/master/core_functions.R) which contains the core functions---such as the Gibbs sampling algorithm described in detail in Section 3 of the paper.

## Shiny application

We also provide a [Shiny application](https://github.com/tommasorigon/India-SequentiaLogit/tree/master/SequentiaLogisticApp) in which, for every combination of covariates characterizing the different women profiles, the posterior probabilities of the corresponding contraceptive preferences, and the sequential probabilities in Figure 1, are computed. The characteristics of each woman (e.g. `AGE`, `CHILD`, etc.) can be defined interactively by the user.

The application relies on three `R` packages: `shiny`, `shinydashboard` and `ggplot2`. If not already available, these packages **should be installed** using the `R` function `install.packages()`. Once these packages are installed, the Shiny application can be launched by executing the following `R` code.

```r
# # Command required to install the necessary packages
# install.packages("shiny")
# install.packages("shinydashboard")
# install.packages("ggplot2")

shiny::runGitHub("India-SequentiaLogit","tommasorigon", subdir = "SequentiaLogisticApp")
```
