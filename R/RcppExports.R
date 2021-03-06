# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

logsumC <- function(x) {
    .Call(`_binMixtC_logsumC`, x)
}

logdiffC <- function(x) {
    .Call(`_binMixtC_logdiffC`, x)
}

ran <- function(n, p, mu, v) {
    .Call(`_binMixtC_ran`, n, p, mu, v)
}

mvrnormArma <- function(n, mu, sigma) {
    .Call(`_binMixtC_mvrnormArma`, n, mu, sigma)
}

genmixtC <- function(n, p, mu, sigma) {
    .Call(`_binMixtC_genmixtC`, n, p, mu, sigma)
}

zerobin1 <- function(bin) {
    .Call(`_binMixtC_zerobin1`, bin)
}

buildbinC <- function(x, R) {
    .Call(`_binMixtC_buildbinC`, x, R)
}

buildbinmargC <- function(X, R) {
    .Call(`_binMixtC_buildbinmargC`, X, R)
}

dmixtC <- function(x, pi, mu, v) {
    .Call(`_binMixtC_dmixtC`, x, pi, mu, v)
}

loglimixtC <- function(x, pi, mu, v) {
    .Call(`_binMixtC_loglimixtC`, x, pi, mu, v)
}

loglimultC <- function(bin, pi, mu, v) {
    .Call(`_binMixtC_loglimultC`, bin, pi, mu, v)
}

emuniv_nomixtC <- function(bin, mu0, sigma0, eps, it) {
    .Call(`_binMixtC_emuniv_nomixtC`, bin, mu0, sigma0, eps, it)
}

emunivC <- function(bin, pi0, mu0, sigma0, eps, it) {
    .Call(`_binMixtC_emunivC`, bin, pi0, mu0, sigma0, eps, it)
}

margmu <- function(mu, dim) {
    .Call(`_binMixtC_margmu`, mu, dim)
}

margv <- function(v, dim) {
    .Call(`_binMixtC_margv`, v, dim)
}

loglimargC <- function(bin, pi, mu, v) {
    .Call(`_binMixtC_loglimargC`, bin, pi, mu, v)
}

emdimC <- function(bin, pi0, mu0, sigma0) {
    .Call(`_binMixtC_emdimC`, bin, pi0, mu0, sigma0)
}

embingauscompC <- function(bin, pi0, mu0, sigma0, eps, it) {
    .Call(`_binMixtC_embingauscompC`, bin, pi0, mu0, sigma0, eps, it)
}

binmixtC <- function(data, cl, R, it, eps, nrep) {
    .Call(`_binMixtC_binmixtC`, data, cl, R, it, eps, nrep)
}

binmixtclassicC <- function(data, cl, R, it, eps, eps1, it1, nrep) {
    .Call(`_binMixtC_binmixtclassicC`, data, cl, R, it, eps, eps1, it1, nrep)
}

tabulate2 <- function(x, max) {
    .Call(`_binMixtC_tabulate2`, x, max)
}

zerobin <- function(bin) {
    .Call(`_binMixtC_zerobin`, bin)
}

buildbinmultC <- function(x, R) {
    .Call(`_binMixtC_buildbinmultC`, x, R)
}

choosecpp <- function(n, k) {
    .Call(`_binMixtC_choosecpp`, n, k)
}

combncpp <- function(N, K) {
    .Call(`_binMixtC_combncpp`, N, K)
}

buildbinbivmargC <- function(X, R) {
    .Call(`_binMixtC_buildbinbivmargC`, X, R)
}

loglimultmbivC <- function(bin, pi, mu, v) {
    .Call(`_binMixtC_loglimultmbivC`, bin, pi, mu, v)
}

margmubiv <- function(mu, dim) {
    .Call(`_binMixtC_margmubiv`, mu, dim)
}

margvbiv <- function(v, dim) {
    .Call(`_binMixtC_margvbiv`, v, dim)
}

loglimargbivC <- function(bin, pi, mu, v) {
    .Call(`_binMixtC_loglimargbivC`, bin, pi, mu, v)
}

emdimbivpimuC <- function(bin, pi0, mu0, sigma0) {
    .Call(`_binMixtC_emdimbivpimuC`, bin, pi0, mu0, sigma0)
}

sigmadimfunC <- function(bin, pcol, p11, p12, a11, a12, b11, b12, pi0, mu0, mu1, sigma0) {
    .Call(`_binMixtC_sigmadimfunC`, bin, pcol, p11, p12, a11, a12, b11, b12, pi0, mu0, mu1, sigma0)
}

embinbivcompC <- function(bin, pi0, mu0, sigma0, eps, it) {
    .Call(`_binMixtC_embinbivcompC`, bin, pi0, mu0, sigma0, eps, it)
}

binmixtbivC <- function(data, cl, R, it, eps, nrep) {
    .Call(`_binMixtC_binmixtbivC`, data, cl, R, it, eps, nrep)
}

binmixtbivclassicC <- function(data, cl, R, it, eps, eps1, it1, nrep) {
    .Call(`_binMixtC_binmixtbivclassicC`, data, cl, R, it, eps, eps1, it1, nrep)
}

