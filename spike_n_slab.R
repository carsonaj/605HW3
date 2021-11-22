library(MASS)
source("var_sel_data_gen.R")
source("var_sel_MCMC.R")
source("eval.R")

spike_n_slab <- function(data) {
    # get data:
    x <- data$x
    y <- data$y

    n <- nrow(x)
    p <- ncol(x)

    # Run MCMC:
    # params
    g <- p^2
    kap <- 2
    burnin <- 1000
    n_samp <- 50

    # init gamma0
    f <- floor(n / 2)
    gam_init <- rep(0, p)
    gam_init[1:f] <- rep(1, f)

    # mcmc
    gams <- mcmc(x, y, g, kap, gam_init, burnin, n_samp)
    gam_probs <- colSums(gams) / n_samp

    gam_hat <- gam_init
    ind <- which(gam_probs >= .5)
    gam_hat[ind] <- 1

    return(gam_hat)
}
