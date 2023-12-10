# SLAM

The `slam` package is used for implementing the Monte Carlo EM algorithm proposed in the paper ``Semiparametric Latent ANOVA Model for Event-Related Potentials". The algorithm estimates the latency of ERP signals with the full posterior distribution at both subject level and group level.


## Installation and Usage

You can install `slam` from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("chenghanyustats/slam")
```

Some important functions are

- `mcem_slam()` implement the proposed MCEM algorithm for SLAM.

- `Rsolnp::solnp()` is used to optimize the hyperparameters of the kernel function used in the GPR.

- `get_pi_t_sig()` does the regression curve fitting and computes the associated credible interval.

- `plot_pred_gp_f_y()` plots the fitted results from `get_pi_t_sig()`.

- `get_amp_*()` family estimates the amplitude of the ERP components.


## Package Dependency

The required R packages include 

    packages <- c("Matrix", "matrixcalc", "emulator", "doParallel", "Rsolnp", 
              "mvnfast", "extraDistr", "truncnorm", "invgamma", "ggridges",
              "ggplot2", "tidyverse")

The R version >= 3.5 is required.



### Vignette

The package provide a vignette `slam.Rmd`(`slam.pdf`) for implementing the proposed method including simulation, and ERP case studies. The code reproduces some key figures shown in the paper.

