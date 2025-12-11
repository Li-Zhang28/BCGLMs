# BCGLMs: Bayesian Compostional GLMs for Analyzing Microbiome Data

# Overview

This R package provide Functions for fitting Bayesian compositional models for analyzing microbiome data, including generalized linear models (GLMs), cumulative logistic model and survival outcome.
We employ a structured regularized horseshoe prior for the compositional coefficients and a soft sum-to-zero restriction on coefficients through the prior distribution. 
The models are fitted using efficient Hamiltonian Monte Carlo algorithm.

Author: Li Zhang [lzhang28@uab.edu](mailto:lzhang28@uab.edu)

# Installation 
For users running R version 4.4.1, the following code can be used to install the required packages and BCGLMs:

##
<pre> ## Install Bioconductor dependency
if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
## Install phyloseq
##Documentation: Check out [Phyloseq R package](https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html) for more versions
BiocManager::install("phyloseq")
## Install remotes for GitHub installation
install.packages("remotes")
## Install brms (development version)
remotes::install_github("paul-buerkner/brms")
## Install BhGLM
## NOTE: BhGLM is currently incompatible with R > 4.4
remotes::install_github("nyiuab/BhGLM", force = TRUE)
## Install BCGLMs
remotes::install_github("Li-Zhang28/BCGLMs", force = TRUE, build_vignettes = TRUE)
 </pre>

# Example
Fitting Bayesian compositional GLMs methods for microbiome data with continuous outcome.

##
<pre> dat = sim_c(n = 400, p = 100) 
          sim = similarity(dat$x) 
          fit = bcglm(x = dat$x, y = dat$y, family = gaussian, df_local = 1, df_global = 1, similarity = sim) 
          summary(fit) 
          fixef(fit) 
          mcmc_plot(fit, variable = "^b_X", regex = TRUE) 
          plot(fit, variable = c("b_X18", "b_X20", "b_X22", "b_X24", "b_X26", "b_X28"), nvariables = 6) </pre>
