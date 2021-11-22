zero_one <- function(gam0, gam_hat) {
    p <- length(gam0)
    s <- 0
    for (i in 1:p) {
        g01 <- (gam_hat[i] == 0) && (gam0[i] == 1)
        g10 <- (gam_hat[i] == 1) && (gam0[i] == 0)

        if (g01 || g10) {
            s <- s + 1
        }
    }

    return(s)
}

supp_size <- function(gam_hat) {
    s <- sum(gam_hat == 1)

    return(s)
}
