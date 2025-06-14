# BCGLMs: Bayesian Compostional GLMs for Analyzing Microbiome Data

# Overview

This R package provide Functions for fitting Bayesian compositional models for analyzing microbiome data, including generalized linear models (GLMs), cumulative logistic model and survival outcome.
We employ a structured regularized horseshoe prior for the compositional coefficients and a soft sum-to-zero restriction on coefficients through the prior distribution. 
The models are fitted using efficient Hamiltonian Monte Carlo algorithm.

Author: Li Zhang [lzhang28@uab.edu](mailto:lzhang28@uab.edu)

# Installation 

Install the released version of remotes from CRAN:
##
    library(remotes)
    install_github("Li-Zhang28/BCGLMs", force=T, build_vignettes=T)
