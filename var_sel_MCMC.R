library(MASS)

single_flip <- function(gam1) {
    p <- length(gam1)
    ind <- sample(1:p, 1)
    gam2 <- gam1
    gam2[ind] <- as.numeric(!gam2[ind])

    return(gam2)
}

double_flip <- function(gam1) {
    inds <- which(gam1 == 1)
    inds_c <- which(gam1 == 0)
    ind <- sample(inds, 1)
    ind_c <- sample(inds_c, 1)

    gam2 <- gam1
    gam2[ind] <- as.numeric(!gam2[ind])
    gam2[ind_c] <- as.numeric(!gam2[ind_c])

    return(gam2)
}

get_candidate_gam2 <- function(gam1) {
    p <- length(gam1)
    p_gam1 <- sum(gam1)
    if (p_gam1 == 0 || p_gam1 == p) {
        gam2 <- single_flip(gam1)
    } else {
        nbhd <- sample(1:2, 1)
        if (nbhd == 1) {
            gam2 <- single_flip(gam1)
        } else {
            gam2 <- double_flip(gam1)
        }
    }

    return(gam2)
}

# "uniform" on gam2 when currently at gam1:
s <- function(gam1) {
    p <- length(gam1)
    inds <- which(gam1 == 1)
    inds_c <- which(gam1 == 0)
    l <- length(inds)
    l_c <- length(inds_c)

    if(l == 0 || l_c == 0) {
        prob <- 1 / p
    } else {
        prob <- .5 / p + .5 / (l * l_c)
    }

    return(prob)
}
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

get_x_gam <- function(x, gam) {
    inds <- which(gam == 1)
    x_gam <- x[ , inds]

    return(x_gam)
}

get_phi_gam <- function(x, gam) {
    n <- nrow(x)
    p_gam <- sum(gam)
    if (p_gam == 0) {
        phi_gam <- matrix(0, n, n)
    } else {
        x_gam <- get_x_gam(x, gam)
        prod_inv <- ginv(crossprod(x_gam, x_gam))
        phi_gam <- tcrossprod(x_gam %*% prod_inv, x_gam)
    }

    return(phi_gam)
}

get_r_gam <- function(y, phi_gam) {
    norm_y_sq <- crossprod(y, y)
    r_gam <- crossprod(y, phi_gam %*% y) / norm_y_sq

    return(r_gam)
}

# likelihood p(y|gam):
likelihood <- function(y, x, gam, g) {
    n <- nrow(x)
    p <- ncol(x)
    p_gam <- sum(gam)
    if (p_gam > n) {
        like <- 0
    } else {
        phi_gam <- get_phi_gam(x, gam)
        r_gam <- get_r_gam(y, phi_gam)

        left <- log(gamma(n / 2)) + (n / 2) * log(1 + g) -
            (n / 2) * log(pi) - (n / 2) * log(crossprod(y, y))
        right <- (p_gam / 2) * log(1 + g) - (n / 2) * log (1 + g * (1 - r_gam^2))
        loglike <- left + right

        like <- exp(loglike)
    }

    return(like)
}

# prior p(gam)
prior <- function(gam, kap, n) {
    p <- length(gam)
    p_gam <- sum(gam)
    indic <- as.numeric(p_gam <= n)
    prob <- 1 / p^(kap * p_gam) * indic

    return(prob)
}

# proportional posterior pi(gam|y):
post <- function(y, x, gam, g, kap) {
    n <- nrow(x)
    p <- ncol(x)
    like_prob <- likelihood(y, x, gam, g)
    prior_prob <- prior(gam, kap, n)

    if (prior_prob == 0) {
        prob <- 0
    } else {
        prob <- like_prob * prior_prob
    }

    return(prob)
}
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# MH accept prob:
get_accept_prob <- function(gam1, gam2, y, x, g, kap) {
    top <- post(y, x, gam2, g, kap) * s(gam2)
    bottom <- post(y, x, gam1, g, kap) * s(gam1)

    if (top == 0) {
        r <- 0
    } else {
        r <- min(1, top / bottom)
    }

    return(r)
}

# get next gam state:
get_gam2 <- function(gam1, y, x, g, kap) {
    gam2 <- get_candidate_gam2(gam1)
    accept_prob <- get_accept_prob(gam1, gam2, y, x, g, kap)
    accept <- sample(1:2, 1, prob = c(1 - accept_prob, accept_prob))

    if (accept == 1) {
        gam2 <- gam1
    }

    return(gam2)
}
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# MCMC sampler for gamma
mcmc <- function(x, y, g, kap, gam_init, burnin, n_samp) {
    p <- ncol(x)
    gam1 <- gam_init
    for (i in 1:burnin) {
        gam1 <- get_gam2(gam1, y, x, g, kap)
    }

    gams <- matrix(0, n_samp, p)

    for (i in 1:n_samp) {
        gam1 <- get_gam2(gam1, y, x, g, kap)
        gams[i, ] <- gam1
    }

    return(gams)
}

