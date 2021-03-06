// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// logsumC
double logsumC(DoubleVector x);
RcppExport SEXP _binMixtC_logsumC(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DoubleVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(logsumC(x));
    return rcpp_result_gen;
END_RCPP
}
// logdiffC
double logdiffC(DoubleVector x);
RcppExport SEXP _binMixtC_logdiffC(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DoubleVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(logdiffC(x));
    return rcpp_result_gen;
END_RCPP
}
// ran
DoubleVector ran(int n, DoubleVector p, DoubleVector mu, DoubleVector v);
RcppExport SEXP _binMixtC_ran(SEXP nSEXP, SEXP pSEXP, SEXP muSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< DoubleVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< DoubleVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< DoubleVector >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(ran(n, p, mu, v));
    return rcpp_result_gen;
END_RCPP
}
// mvrnormArma
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma);
RcppExport SEXP _binMixtC_mvrnormArma(SEXP nSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(mvrnormArma(n, mu, sigma));
    return rcpp_result_gen;
END_RCPP
}
// genmixtC
List genmixtC(int n, DoubleVector p, List mu, List sigma);
RcppExport SEXP _binMixtC_genmixtC(SEXP nSEXP, SEXP pSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< DoubleVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< List >::type mu(muSEXP);
    Rcpp::traits::input_parameter< List >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(genmixtC(n, p, mu, sigma));
    return rcpp_result_gen;
END_RCPP
}
// zerobin1
List zerobin1(List bin);
RcppExport SEXP _binMixtC_zerobin1(SEXP binSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type bin(binSEXP);
    rcpp_result_gen = Rcpp::wrap(zerobin1(bin));
    return rcpp_result_gen;
END_RCPP
}
// buildbinC
List buildbinC(arma::vec x, int R);
RcppExport SEXP _binMixtC_buildbinC(SEXP xSEXP, SEXP RSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type R(RSEXP);
    rcpp_result_gen = Rcpp::wrap(buildbinC(x, R));
    return rcpp_result_gen;
END_RCPP
}
// buildbinmargC
List buildbinmargC(arma::mat X, arma::vec R);
RcppExport SEXP _binMixtC_buildbinmargC(SEXP XSEXP, SEXP RSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type R(RSEXP);
    rcpp_result_gen = Rcpp::wrap(buildbinmargC(X, R));
    return rcpp_result_gen;
END_RCPP
}
// dmixtC
arma::vec dmixtC(arma::vec x, arma::vec pi, arma::vec mu, arma::vec v);
RcppExport SEXP _binMixtC_dmixtC(SEXP xSEXP, SEXP piSEXP, SEXP muSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pi(piSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(dmixtC(x, pi, mu, v));
    return rcpp_result_gen;
END_RCPP
}
// loglimixtC
double loglimixtC(arma::vec x, arma::vec pi, arma::vec mu, arma::vec v);
RcppExport SEXP _binMixtC_loglimixtC(SEXP xSEXP, SEXP piSEXP, SEXP muSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pi(piSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(loglimixtC(x, pi, mu, v));
    return rcpp_result_gen;
END_RCPP
}
// loglimultC
double loglimultC(Rcpp::List bin, arma::vec pi, arma::vec mu, arma::vec v);
RcppExport SEXP _binMixtC_loglimultC(SEXP binSEXP, SEXP piSEXP, SEXP muSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type bin(binSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pi(piSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(loglimultC(bin, pi, mu, v));
    return rcpp_result_gen;
END_RCPP
}
// emuniv_nomixtC
List emuniv_nomixtC(List bin, arma::vec mu0, arma::vec sigma0, double eps, int it);
RcppExport SEXP _binMixtC_emuniv_nomixtC(SEXP binSEXP, SEXP mu0SEXP, SEXP sigma0SEXP, SEXP epsSEXP, SEXP itSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type bin(binSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sigma0(sigma0SEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type it(itSEXP);
    rcpp_result_gen = Rcpp::wrap(emuniv_nomixtC(bin, mu0, sigma0, eps, it));
    return rcpp_result_gen;
END_RCPP
}
// emunivC
List emunivC(List bin, arma::vec pi0, arma::vec mu0, arma::vec sigma0, double eps, int it);
RcppExport SEXP _binMixtC_emunivC(SEXP binSEXP, SEXP pi0SEXP, SEXP mu0SEXP, SEXP sigma0SEXP, SEXP epsSEXP, SEXP itSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type bin(binSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pi0(pi0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sigma0(sigma0SEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type it(itSEXP);
    rcpp_result_gen = Rcpp::wrap(emunivC(bin, pi0, mu0, sigma0, eps, it));
    return rcpp_result_gen;
END_RCPP
}
// margmu
arma::vec margmu(List mu, int dim);
RcppExport SEXP _binMixtC_margmu(SEXP muSEXP, SEXP dimSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type mu(muSEXP);
    Rcpp::traits::input_parameter< int >::type dim(dimSEXP);
    rcpp_result_gen = Rcpp::wrap(margmu(mu, dim));
    return rcpp_result_gen;
END_RCPP
}
// margv
arma::vec margv(List v, int dim);
RcppExport SEXP _binMixtC_margv(SEXP vSEXP, SEXP dimSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type dim(dimSEXP);
    rcpp_result_gen = Rcpp::wrap(margv(v, dim));
    return rcpp_result_gen;
END_RCPP
}
// loglimargC
double loglimargC(Rcpp::List bin, arma::vec pi, List mu, List v);
RcppExport SEXP _binMixtC_loglimargC(SEXP binSEXP, SEXP piSEXP, SEXP muSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type bin(binSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pi(piSEXP);
    Rcpp::traits::input_parameter< List >::type mu(muSEXP);
    Rcpp::traits::input_parameter< List >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(loglimargC(bin, pi, mu, v));
    return rcpp_result_gen;
END_RCPP
}
// emdimC
List emdimC(List bin, arma::vec pi0, arma::vec mu0, arma::vec sigma0);
RcppExport SEXP _binMixtC_emdimC(SEXP binSEXP, SEXP pi0SEXP, SEXP mu0SEXP, SEXP sigma0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type bin(binSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pi0(pi0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sigma0(sigma0SEXP);
    rcpp_result_gen = Rcpp::wrap(emdimC(bin, pi0, mu0, sigma0));
    return rcpp_result_gen;
END_RCPP
}
// embingauscompC
List embingauscompC(List bin, arma::vec pi0, List mu0, List sigma0, double eps, int it);
RcppExport SEXP _binMixtC_embingauscompC(SEXP binSEXP, SEXP pi0SEXP, SEXP mu0SEXP, SEXP sigma0SEXP, SEXP epsSEXP, SEXP itSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type bin(binSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pi0(pi0SEXP);
    Rcpp::traits::input_parameter< List >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< List >::type sigma0(sigma0SEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type it(itSEXP);
    rcpp_result_gen = Rcpp::wrap(embingauscompC(bin, pi0, mu0, sigma0, eps, it));
    return rcpp_result_gen;
END_RCPP
}
// binmixtC
List binmixtC(arma::mat data, int cl, arma::vec R, int it, double eps, int nrep);
RcppExport SEXP _binMixtC_binmixtC(SEXP dataSEXP, SEXP clSEXP, SEXP RSEXP, SEXP itSEXP, SEXP epsSEXP, SEXP nrepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type cl(clSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type R(RSEXP);
    Rcpp::traits::input_parameter< int >::type it(itSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type nrep(nrepSEXP);
    rcpp_result_gen = Rcpp::wrap(binmixtC(data, cl, R, it, eps, nrep));
    return rcpp_result_gen;
END_RCPP
}
// binmixtclassicC
List binmixtclassicC(arma::mat data, int cl, arma::vec R, int it, double eps, double eps1, int it1, int nrep);
RcppExport SEXP _binMixtC_binmixtclassicC(SEXP dataSEXP, SEXP clSEXP, SEXP RSEXP, SEXP itSEXP, SEXP epsSEXP, SEXP eps1SEXP, SEXP it1SEXP, SEXP nrepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type cl(clSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type R(RSEXP);
    Rcpp::traits::input_parameter< int >::type it(itSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< double >::type eps1(eps1SEXP);
    Rcpp::traits::input_parameter< int >::type it1(it1SEXP);
    Rcpp::traits::input_parameter< int >::type nrep(nrepSEXP);
    rcpp_result_gen = Rcpp::wrap(binmixtclassicC(data, cl, R, it, eps, eps1, it1, nrep));
    return rcpp_result_gen;
END_RCPP
}
// tabulate2
NumericVector tabulate2(NumericVector x, int max);
RcppExport SEXP _binMixtC_tabulate2(SEXP xSEXP, SEXP maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type max(maxSEXP);
    rcpp_result_gen = Rcpp::wrap(tabulate2(x, max));
    return rcpp_result_gen;
END_RCPP
}
// zerobin
List zerobin(List bin);
RcppExport SEXP _binMixtC_zerobin(SEXP binSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type bin(binSEXP);
    rcpp_result_gen = Rcpp::wrap(zerobin(bin));
    return rcpp_result_gen;
END_RCPP
}
// buildbinmultC
List buildbinmultC(arma::mat x, int R);
RcppExport SEXP _binMixtC_buildbinmultC(SEXP xSEXP, SEXP RSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type R(RSEXP);
    rcpp_result_gen = Rcpp::wrap(buildbinmultC(x, R));
    return rcpp_result_gen;
END_RCPP
}
// choosecpp
int choosecpp(int n, int k);
RcppExport SEXP _binMixtC_choosecpp(SEXP nSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(choosecpp(n, k));
    return rcpp_result_gen;
END_RCPP
}
// combncpp
arma::mat combncpp(int N, int K);
RcppExport SEXP _binMixtC_combncpp(SEXP NSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(combncpp(N, K));
    return rcpp_result_gen;
END_RCPP
}
// buildbinbivmargC
List buildbinbivmargC(arma::mat X, arma::vec R);
RcppExport SEXP _binMixtC_buildbinbivmargC(SEXP XSEXP, SEXP RSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type R(RSEXP);
    rcpp_result_gen = Rcpp::wrap(buildbinbivmargC(X, R));
    return rcpp_result_gen;
END_RCPP
}
// loglimultmbivC
double loglimultmbivC(Rcpp::List bin, arma::vec pi, Rcpp::List mu, Rcpp::List v);
RcppExport SEXP _binMixtC_loglimultmbivC(SEXP binSEXP, SEXP piSEXP, SEXP muSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type bin(binSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pi(piSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type mu(muSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(loglimultmbivC(bin, pi, mu, v));
    return rcpp_result_gen;
END_RCPP
}
// margmubiv
List margmubiv(List mu, arma::vec dim);
RcppExport SEXP _binMixtC_margmubiv(SEXP muSEXP, SEXP dimSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type dim(dimSEXP);
    rcpp_result_gen = Rcpp::wrap(margmubiv(mu, dim));
    return rcpp_result_gen;
END_RCPP
}
// margvbiv
List margvbiv(List v, arma::vec dim);
RcppExport SEXP _binMixtC_margvbiv(SEXP vSEXP, SEXP dimSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type v(vSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type dim(dimSEXP);
    rcpp_result_gen = Rcpp::wrap(margvbiv(v, dim));
    return rcpp_result_gen;
END_RCPP
}
// loglimargbivC
double loglimargbivC(Rcpp::List bin, arma::vec pi, List mu, List v);
RcppExport SEXP _binMixtC_loglimargbivC(SEXP binSEXP, SEXP piSEXP, SEXP muSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type bin(binSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pi(piSEXP);
    Rcpp::traits::input_parameter< List >::type mu(muSEXP);
    Rcpp::traits::input_parameter< List >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(loglimargbivC(bin, pi, mu, v));
    return rcpp_result_gen;
END_RCPP
}
// emdimbivpimuC
List emdimbivpimuC(Rcpp::List bin, arma::vec pi0, List mu0, List sigma0);
RcppExport SEXP _binMixtC_emdimbivpimuC(SEXP binSEXP, SEXP pi0SEXP, SEXP mu0SEXP, SEXP sigma0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type bin(binSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pi0(pi0SEXP);
    Rcpp::traits::input_parameter< List >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< List >::type sigma0(sigma0SEXP);
    rcpp_result_gen = Rcpp::wrap(emdimbivpimuC(bin, pi0, mu0, sigma0));
    return rcpp_result_gen;
END_RCPP
}
// sigmadimfunC
List sigmadimfunC(List bin, arma::vec pcol, arma::mat p11, arma::mat p12, arma::mat a11, arma::mat a12, arma::mat b11, arma::mat b12, arma::vec pi0, List mu0, List mu1, List sigma0);
RcppExport SEXP _binMixtC_sigmadimfunC(SEXP binSEXP, SEXP pcolSEXP, SEXP p11SEXP, SEXP p12SEXP, SEXP a11SEXP, SEXP a12SEXP, SEXP b11SEXP, SEXP b12SEXP, SEXP pi0SEXP, SEXP mu0SEXP, SEXP mu1SEXP, SEXP sigma0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type bin(binSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pcol(pcolSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type p11(p11SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type p12(p12SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type a11(a11SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type a12(a12SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type b11(b11SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type b12(b12SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pi0(pi0SEXP);
    Rcpp::traits::input_parameter< List >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< List >::type mu1(mu1SEXP);
    Rcpp::traits::input_parameter< List >::type sigma0(sigma0SEXP);
    rcpp_result_gen = Rcpp::wrap(sigmadimfunC(bin, pcol, p11, p12, a11, a12, b11, b12, pi0, mu0, mu1, sigma0));
    return rcpp_result_gen;
END_RCPP
}
// embinbivcompC
List embinbivcompC(List bin, arma::vec pi0, List mu0, List sigma0, double eps, int it);
RcppExport SEXP _binMixtC_embinbivcompC(SEXP binSEXP, SEXP pi0SEXP, SEXP mu0SEXP, SEXP sigma0SEXP, SEXP epsSEXP, SEXP itSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type bin(binSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pi0(pi0SEXP);
    Rcpp::traits::input_parameter< List >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< List >::type sigma0(sigma0SEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type it(itSEXP);
    rcpp_result_gen = Rcpp::wrap(embinbivcompC(bin, pi0, mu0, sigma0, eps, it));
    return rcpp_result_gen;
END_RCPP
}
// binmixtbivC
List binmixtbivC(arma::mat data, int cl, arma::vec R, int it, double eps, int nrep);
RcppExport SEXP _binMixtC_binmixtbivC(SEXP dataSEXP, SEXP clSEXP, SEXP RSEXP, SEXP itSEXP, SEXP epsSEXP, SEXP nrepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type cl(clSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type R(RSEXP);
    Rcpp::traits::input_parameter< int >::type it(itSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type nrep(nrepSEXP);
    rcpp_result_gen = Rcpp::wrap(binmixtbivC(data, cl, R, it, eps, nrep));
    return rcpp_result_gen;
END_RCPP
}
// binmixtbivclassicC
List binmixtbivclassicC(arma::mat data, int cl, arma::vec R, int it, double eps, double eps1, int it1, int nrep);
RcppExport SEXP _binMixtC_binmixtbivclassicC(SEXP dataSEXP, SEXP clSEXP, SEXP RSEXP, SEXP itSEXP, SEXP epsSEXP, SEXP eps1SEXP, SEXP it1SEXP, SEXP nrepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type cl(clSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type R(RSEXP);
    Rcpp::traits::input_parameter< int >::type it(itSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< double >::type eps1(eps1SEXP);
    Rcpp::traits::input_parameter< int >::type it1(it1SEXP);
    Rcpp::traits::input_parameter< int >::type nrep(nrepSEXP);
    rcpp_result_gen = Rcpp::wrap(binmixtbivclassicC(data, cl, R, it, eps, eps1, it1, nrep));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_binMixtC_logsumC", (DL_FUNC) &_binMixtC_logsumC, 1},
    {"_binMixtC_logdiffC", (DL_FUNC) &_binMixtC_logdiffC, 1},
    {"_binMixtC_ran", (DL_FUNC) &_binMixtC_ran, 4},
    {"_binMixtC_mvrnormArma", (DL_FUNC) &_binMixtC_mvrnormArma, 3},
    {"_binMixtC_genmixtC", (DL_FUNC) &_binMixtC_genmixtC, 4},
    {"_binMixtC_zerobin1", (DL_FUNC) &_binMixtC_zerobin1, 1},
    {"_binMixtC_buildbinC", (DL_FUNC) &_binMixtC_buildbinC, 2},
    {"_binMixtC_buildbinmargC", (DL_FUNC) &_binMixtC_buildbinmargC, 2},
    {"_binMixtC_dmixtC", (DL_FUNC) &_binMixtC_dmixtC, 4},
    {"_binMixtC_loglimixtC", (DL_FUNC) &_binMixtC_loglimixtC, 4},
    {"_binMixtC_loglimultC", (DL_FUNC) &_binMixtC_loglimultC, 4},
    {"_binMixtC_emuniv_nomixtC", (DL_FUNC) &_binMixtC_emuniv_nomixtC, 5},
    {"_binMixtC_emunivC", (DL_FUNC) &_binMixtC_emunivC, 6},
    {"_binMixtC_margmu", (DL_FUNC) &_binMixtC_margmu, 2},
    {"_binMixtC_margv", (DL_FUNC) &_binMixtC_margv, 2},
    {"_binMixtC_loglimargC", (DL_FUNC) &_binMixtC_loglimargC, 4},
    {"_binMixtC_emdimC", (DL_FUNC) &_binMixtC_emdimC, 4},
    {"_binMixtC_embingauscompC", (DL_FUNC) &_binMixtC_embingauscompC, 6},
    {"_binMixtC_binmixtC", (DL_FUNC) &_binMixtC_binmixtC, 6},
    {"_binMixtC_binmixtclassicC", (DL_FUNC) &_binMixtC_binmixtclassicC, 8},
    {"_binMixtC_tabulate2", (DL_FUNC) &_binMixtC_tabulate2, 2},
    {"_binMixtC_zerobin", (DL_FUNC) &_binMixtC_zerobin, 1},
    {"_binMixtC_buildbinmultC", (DL_FUNC) &_binMixtC_buildbinmultC, 2},
    {"_binMixtC_choosecpp", (DL_FUNC) &_binMixtC_choosecpp, 2},
    {"_binMixtC_combncpp", (DL_FUNC) &_binMixtC_combncpp, 2},
    {"_binMixtC_buildbinbivmargC", (DL_FUNC) &_binMixtC_buildbinbivmargC, 2},
    {"_binMixtC_loglimultmbivC", (DL_FUNC) &_binMixtC_loglimultmbivC, 4},
    {"_binMixtC_margmubiv", (DL_FUNC) &_binMixtC_margmubiv, 2},
    {"_binMixtC_margvbiv", (DL_FUNC) &_binMixtC_margvbiv, 2},
    {"_binMixtC_loglimargbivC", (DL_FUNC) &_binMixtC_loglimargbivC, 4},
    {"_binMixtC_emdimbivpimuC", (DL_FUNC) &_binMixtC_emdimbivpimuC, 4},
    {"_binMixtC_sigmadimfunC", (DL_FUNC) &_binMixtC_sigmadimfunC, 12},
    {"_binMixtC_embinbivcompC", (DL_FUNC) &_binMixtC_embinbivcompC, 6},
    {"_binMixtC_binmixtbivC", (DL_FUNC) &_binMixtC_binmixtbivC, 6},
    {"_binMixtC_binmixtbivclassicC", (DL_FUNC) &_binMixtC_binmixtbivclassicC, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_binMixtC(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
