#include "evomod.h"

namespace evomod {

arma::vec rnorm_vec(arma::uword n) {
  arma::vec out(n);
  for (arma::uword i = 0; i < n; ++i) out[i] = R::norm_rand();
  return out;
}

arma::mat locus_effects(const arma::mat& alleles) {
  arma::uword loci = alleles.n_cols / 2;
  arma::mat Y(alleles.n_rows, loci);
  for (arma::uword l = 0; l < loci; ++l)
    Y.col(l) = alleles.col(2 * l) + alleles.col(2 * l + 1);
  return Y;
}

// ---- Additive ---------------------------------------------------------------

arma::vec Additive::genotypic_value(const arma::mat& alleles) const {
  return arma::sum(alleles, 1);
}

arma::mat Additive::local_effects(const arma::mat& alleles) const {
  return locus_effects(alleles);
}

// ---- Multilinear ------------------------------------------------------------

arma::vec Multilinear::bilinear(const arma::vec& a, const arma::vec& b) const {
  arma::vec out(epsilon_.n_slices);
  for (arma::uword k = 0; k < epsilon_.n_slices; ++k)
    out[k] = arma::as_scalar(a.t() * epsilon_.slice(k) * b);
  return out;
}

arma::vec Multilinear::genotypic_value(const arma::mat& alleles) const {
  arma::mat Y = locus_effects(alleles);
  arma::vec S = arma::sum(Y, 1);
  arma::mat YY = Y * Y.t();
  arma::vec x = S;
  for (arma::uword k = 0; k < epsilon_.n_slices; ++k) {
    const arma::mat& E = epsilon_.slice(k);
    x[k] += 0.5 * (arma::as_scalar(S.t() * E * S) - arma::accu(E % YY));
  }
  return x;
}

arma::mat Multilinear::local_effects(const arma::mat& alleles) const {
  arma::mat Y = locus_effects(alleles);
  arma::vec S = arma::sum(Y, 1);
  arma::mat L(Y.n_rows, Y.n_cols);
  for (arma::uword i = 0; i < Y.n_cols; ++i)
    L.col(i) = Y.col(i) + bilinear(Y.col(i), S - Y.col(i));
  return L;
}

// ---- Mutation ---------------------------------------------------------------

GaussianMutation::GaussianMutation(double rate, const arma::mat& cov) : rate_(rate) {
  if (!arma::chol(chol_lower_, cov, "lower"))
    Rcpp::stop("mutation covariance matrix is not positive definite");
}

void GaussianMutation::mutate(arma::mat& alleles) const {
  arma::uword n_alleles = alleles.n_cols;
  int n_mut = static_cast<int>(R::rbinom(static_cast<double>(n_alleles), rate_));
  for (int i = 0; i < n_mut; ++i) {
    arma::uword a = static_cast<arma::uword>(R::unif_rand() * n_alleles);
    if (a == n_alleles) a = n_alleles - 1;
    alleles.col(a) += chol_lower_ * rnorm_vec(alleles.n_rows);
  }
}

// ---- Fitness ----------------------------------------------------------------

GaussianFitness::GaussianFitness(const arma::vec& theta, const arma::mat& omega)
  : theta_(theta) {
  if (!arma::inv_sympd(omega_inv_, omega))
    Rcpp::stop("omega is not positive definite");
}

arma::vec GaussianFitness::log_fitness(const arma::mat& z) const {
  arma::mat d = z.each_col() - theta_;
  return -0.5 * arma::sum(d % (omega_inv_ * d), 0).t();
}

// ---- Factories ----------------------------------------------------------------

std::unique_ptr<Architecture> make_architecture(const Rcpp::List& spec) {
  std::string type = Rcpp::as<std::string>(spec["type"]);
  if (type == "additive") return std::unique_ptr<Architecture>(new Additive());
  if (type == "multilinear") {
    arma::cube epsilon = Rcpp::as<arma::cube>(spec["epsilon_slices"]);
    return std::unique_ptr<Architecture>(new Multilinear(epsilon));
  }
  Rcpp::stop("unknown architecture type: " + type);
}

std::unique_ptr<Mutation> make_mutation(const Rcpp::List& spec) {
  std::string type = Rcpp::as<std::string>(spec["type"]);
  if (type == "gaussian")
    return std::unique_ptr<Mutation>(new GaussianMutation(
      Rcpp::as<double>(spec["rate"]), Rcpp::as<arma::mat>(spec["cov"])));
  Rcpp::stop("unknown mutation type: " + type);
}

}  // namespace evomod
