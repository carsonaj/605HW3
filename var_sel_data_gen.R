get_mse <- function(z) {
    zbar <- mean(z)
    se <- sum(z-zbar)^2
    mse <- mean(se)

    return(mse)
}

get_sig <- function(p, rho) {
    sig <- matrix(rho, p, p)
    diag(sig) <- rep(1, p)

    return(sig)
}

get_x <- function(n, p, sig) {
    x <- mvnfast::rmvn(n, rep(0, p), sig)
    x <- scale(x)

    return(x)
}

get_b0 <- function(p, k0) {
    b0 <- rep(0, p)
    b0[1:k0] <- 1

    return(b0)
}

get_f <- function(x, b0) {
    f <- as.vector(x %*% b0)

    return(f)
}

get_ep <- function(n) {
    ep <- as.vector(mvnfast::rmvn(1, rep(0, n), diag(n)))

    return(ep)
}

get_sig0 <- function(snr, f, ep) {
    mse_f <- get_mse(f)
    mse_ep <- get_mse(ep)
    sig0 <- sqrt(mse_f / (snr * mse_ep))

    return(sig0)
}

get_y <- function(f, sig0, ep) {
    y <- f + sig0 * ep

    return(y)
}

get_data <- function(n, p, snr) {
    rho <- .5
    k0 <- 10

    sig <- get_sig(p, rho)
    x <- get_x(n, p, sig)
    b0 <- get_b0(p, k0)
    f <- get_f(x, b0)
    ep <- get_ep(n)
    sig0 <- get_sig0(snr, f, ep)
    y <- get_y(f, sig0, ep)

    data <- list("x" = x, "y" = y, "b0" = b0, "sig" = sig, "f" = f,
                      "ep" = ep, "sig0" = sig0)

    return(data)


}
