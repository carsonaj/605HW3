library(MASS)
source("var_sel_data_gen.R")
source("var_sel_MCMC.R")

# Data gen for tests:
set.seed(1)
rho <- .5
k0 <- 10
n <- 100
p <- 200
snr <- 2

sig <- get_sig(p, rho)
x <- get_x(n, p, sig)
b0 <- get_b0(p, k0)
f <- get_f(x, b0)
ep <- get_ep(n)
sig0 <- get_sig0(snr, f, ep)
y <- get_y(f, sig0, ep)

# test for get_candidate_gam2:
gam1 <- rep(0, 5)
gam2_cand <- get_candidate_gam2(gam1)
gam2_cand

gam1 <- rep(1, 5)
gam2_cand <- get_candidate_gam2(gam1)
gam2_cand

gam1 <- c(1, 0, 0, 1, 1, 0)
gam2_cand <- get_candidate_gam2(gam1)
gam2_cand

# test for get_x_gam:
gam0 <- rep(0, p)
gam1 <- rep(1, p)
x_gam1 <- get_x_gam(x, gam0)
x_gam0 <- get_x_gam(x, gam1)

# test for get_phi_gam:
phi_gam0 <- get_phi_gam(x, gam0)
phi_gam1 <- get_phi_gam(x, gam1)

# test for get_r_gam:
r_gam0 <- get_r_gam(y, phi_gam0)
r_gam1 <- get_r_gam(y, phi_gam1)

# test for likelihood:
l0 <- likelihood(y, x, gam0, g)
l1 <- likelihood(y, x, gam1, g)

# test for

# test for get_accept_prob:
gam1 <- rep(0, p)
gam1[1:floor(n / 2)] <- rep(1, floor(n / 2))
gam2 <- rep(0, p)
gam2[1:(n+1)] <- rep(1, n+1)

accept_prob <- get_accept_prob(gam1, gam2, y, x, g, kap)
accept_prob

