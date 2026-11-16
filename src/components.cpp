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

// ---- BinaryB -----------------------------------------------------------------

arma::vec BinaryB::genotypic_value(const arma::mat& alleles) const {
  arma::uword p = alleles.n_rows - 1;
  return alleles.rows(1, p) * alleles.row(0).t();
}

arma::mat BinaryB::allele_effects(const arma::mat& alleles) const {
  arma::uword p = alleles.n_rows - 1;
  arma::mat E = alleles.rows(1, p);
  E.each_row() %= alleles.row(0);
  return E;
}

arma::mat BinaryB::local_effects(const arma::mat& alleles) const {
  arma::uword p = alleles.n_rows - 1, loci = alleles.n_cols / 2;
  arma::mat L(p, loci);
  for (arma::uword l = 0; l < loci; ++l)
    L.col(l) = alleles(arma::span(1, p), 2 * l) * alleles(0, 2 * l) +
               alleles(arma::span(1, p), 2 * l + 1) * alleles(0, 2 * l + 1);
  return L;
}

// ---- ContributionMatrix ---------------------------------------------------------

void ContributionMatrix::unpack(const arma::mat& alleles, arma::vec& a, arma::mat& C) const {
  a.set_size(k_);
  C.set_size(p_, k_);
  for (arma::uword j = 0; j < k_; ++j) a[j] = alleles(0, 2 * j) + alleles(0, 2 * j + 1);
  for (arma::uword j = 0; j < k_; ++j)
    for (arma::uword i = 0; i < p_; ++i) {
      arma::uword l = k_ + j * p_ + i;
      C(i, j) = alleles(0, 2 * l) + alleles(0, 2 * l + 1);
    }
}

arma::vec ContributionMatrix::genotypic_value(const arma::mat& alleles) const {
  arma::vec a;
  arma::mat C;
  unpack(alleles, a, C);
  return C * a;
}

arma::mat ContributionMatrix::local_effects(const arma::mat& alleles) const {
  arma::vec a;
  arma::mat C;
  unpack(alleles, a, C);
  arma::mat L(p_, k_ + p_ * k_, arma::fill::zeros);
  for (arma::uword j = 0; j < k_; ++j) {
    L.col(j) = C.col(j) * a[j];                       // weight gene j set to 0
    for (arma::uword i = 0; i < p_; ++i)
      L(i, k_ + j * p_ + i) = C(i, j) * a[j];         // contribution c_ij set to 0
  }
  return L;
}

// ---- Mutation ---------------------------------------------------------------

GaussianMutation::GaussianMutation(double rate, const arma::mat& cov, double persistence)
  : rate_(rate), persistence_(persistence) {
  if (!arma::chol(chol_lower_, cov, "lower"))
    Rcpp::stop("mutation covariance matrix is not positive definite");
}

void GaussianMutation::mutate(arma::mat& alleles) const {
  arma::uword n_alleles = alleles.n_cols;
  int n_mut = static_cast<int>(R::rbinom(static_cast<double>(n_alleles), rate_));
  for (int i = 0; i < n_mut; ++i) {
    arma::uword a = static_cast<arma::uword>(R::unif_rand() * n_alleles);
    if (a == n_alleles) a = n_alleles - 1;
    if (persistence_ == 1.0) alleles.col(a) += chol_lower_ * rnorm_vec(alleles.n_rows);
    else alleles.col(a) = persistence_ * alleles.col(a) + chol_lower_ * rnorm_vec(alleles.n_rows);
  }
}

void GaussianMutation::mutate_one(arma::mat& alleles) const {
  arma::uword n_alleles = alleles.n_cols;
  arma::uword a = std::min(static_cast<arma::uword>(R::unif_rand() * n_alleles), n_alleles - 1);
  if (persistence_ == 1.0) alleles.col(a) += chol_lower_ * rnorm_vec(alleles.n_rows);
  else alleles.col(a) = persistence_ * alleles.col(a) + chol_lower_ * rnorm_vec(alleles.n_rows);
}

void BinaryBMutation::mutate_one(arma::mat& alleles) const {
  arma::uword n_alleles = alleles.n_cols, p = alleles.n_rows - 1;
  double y_events = rate_y_ * n_alleles, b_events = rate_b_ * n_alleles * p;
  if (R::unif_rand() * (y_events + b_events) < y_events) {
    arma::uword a = std::min(static_cast<arma::uword>(R::unif_rand() * n_alleles), n_alleles - 1);
    double old = persistence_ == 1.0 ? alleles(0, a) : persistence_ * alleles(0, a);
    alleles(0, a) = old + sd_y_ * R::norm_rand();
  } else {
    arma::uword n_entries = n_alleles * p;
    arma::uword pos = std::min(static_cast<arma::uword>(R::unif_rand() * n_entries), n_entries - 1);
    arma::uword a = pos / p, k = 1 + pos % p;
    if (!random_state_ || R::unif_rand() < 0.5) alleles(k, a) = 1.0 - alleles(k, a);
  }
}

void BinaryBMutation::mutate(arma::mat& alleles) const {
  arma::uword n_alleles = alleles.n_cols, p = alleles.n_rows - 1;
  int n_y = static_cast<int>(R::rbinom(static_cast<double>(n_alleles), rate_y_));
  for (int i = 0; i < n_y; ++i) {
    arma::uword a = std::min(static_cast<arma::uword>(R::unif_rand() * n_alleles), n_alleles - 1);
    double old = persistence_ == 1.0 ? alleles(0, a) : persistence_ * alleles(0, a);
    alleles(0, a) = old + sd_y_ * R::norm_rand();
  }
  arma::uword n_entries = n_alleles * p;
  int n_b = static_cast<int>(R::rbinom(static_cast<double>(n_entries), rate_b_));
  for (int i = 0; i < n_b; ++i) {
    arma::uword pos = std::min(static_cast<arma::uword>(R::unif_rand() * n_entries), n_entries - 1);
    arma::uword a = pos / p, k = 1 + pos % p;
    if (!random_state_ || R::unif_rand() < 0.5) alleles(k, a) = 1.0 - alleles(k, a);
  }
}

void ContributionMutation::mutate(arma::mat& alleles) const {
  arma::uword w = 2 * weight_loci_, c = alleles.n_cols - w;
  int n_w = static_cast<int>(R::rbinom(static_cast<double>(w), rate_w_));
  for (int i = 0; i < n_w; ++i) {
    arma::uword a = std::min(static_cast<arma::uword>(R::unif_rand() * w), w - 1);
    double old = persistence_ == 1.0 ? alleles(0, a) : persistence_ * alleles(0, a);
    alleles(0, a) = old + sd_ * R::norm_rand();
  }
  int n_c = static_cast<int>(R::rbinom(static_cast<double>(c), rate_c_));
  for (int i = 0; i < n_c; ++i) {
    arma::uword a = w + std::min(static_cast<arma::uword>(R::unif_rand() * c), c - 1);
    double old = persistence_ == 1.0 ? alleles(0, a) : persistence_ * alleles(0, a);
    alleles(0, a) = old + sd_ * R::norm_rand();
  }
}

void ContributionMutation::mutate_one(arma::mat& alleles) const {
  arma::uword w = 2 * weight_loci_, c = alleles.n_cols - w;
  double w_events = rate_w_ * w, c_events = rate_c_ * c;
  arma::uword a;
  if (R::unif_rand() * (w_events + c_events) < w_events)
    a = std::min(static_cast<arma::uword>(R::unif_rand() * w), w - 1);
  else
    a = w + std::min(static_cast<arma::uword>(R::unif_rand() * c), c - 1);
  double old = persistence_ == 1.0 ? alleles(0, a) : persistence_ * alleles(0, a);
  alleles(0, a) = old + sd_ * R::norm_rand();
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

// ---- Recombination ------------------------------------------------------------

void FreeRecombination::offspring(const arma::mat& sire, const arma::mat& dame,
                                  arma::mat& kid) const {
  arma::uword loci = kid.n_cols / 2;
  for (arma::uword l = 0; l < loci; ++l) {
    kid.col(2 * l)     = sire.col(2 * l + (R::unif_rand() < 0.5 ? 0 : 1));
    kid.col(2 * l + 1) = dame.col(2 * l + (R::unif_rand() < 0.5 ? 0 : 1));
  }
}

void LinkageMap::gamete(const arma::mat& parent, arma::mat& kid, arma::uword slot) const {
  arma::uword loci = kid.n_cols / 2;
  arma::uword h = R::unif_rand() < 0.5 ? 0 : 1;
  for (arma::uword l = 0; l < loci; ++l) {
    if (l > 0 && R::unif_rand() < rates_[l - 1]) h = 1 - h;
    kid.col(2 * l + slot) = parent.col(2 * l + h);
  }
}

void LinkageMap::offspring(const arma::mat& sire, const arma::mat& dame, arma::mat& kid) const {
  gamete(sire, kid, 0);
  gamete(dame, kid, 1);
}

// ---- Factories ----------------------------------------------------------------

std::unique_ptr<Architecture> make_architecture(const Rcpp::List& spec) {
  std::string type = Rcpp::as<std::string>(spec["type"]);
  if (type == "additive") return std::unique_ptr<Architecture>(new Additive());
  if (type == "binary_b") return std::unique_ptr<Architecture>(new BinaryB());
  if (type == "contribution_matrix")
    return std::unique_ptr<Architecture>(new ContributionMatrix(
      Rcpp::as<arma::uword>(spec["traits"]), Rcpp::as<arma::uword>(spec["genes"])));
  if (type == "multilinear") {
    arma::cube epsilon = Rcpp::as<arma::cube>(spec["epsilon_slices"]);
    return std::unique_ptr<Architecture>(new Multilinear(epsilon));
  }
  Rcpp::stop("unknown architecture type: " + type);
}

std::unique_ptr<Mutation> make_mutation(const Rcpp::List& spec) {
  std::string type = Rcpp::as<std::string>(spec["type"]);
  // models saved before persistence existed are random walks
  double persistence = spec.containsElementNamed("persistence") ?
    Rcpp::as<double>(spec["persistence"]) : 1.0;
  if (type == "gaussian")
    return std::unique_ptr<Mutation>(new GaussianMutation(
      Rcpp::as<double>(spec["rate"]), Rcpp::as<arma::mat>(spec["cov"]), persistence));
  if (type == "binary_b")
    return std::unique_ptr<Mutation>(new BinaryBMutation(
      Rcpp::as<double>(spec["rate"]), Rcpp::as<double>(spec["sd"]),
      Rcpp::as<double>(spec["rate_b"]), Rcpp::as<bool>(spec["random_state"]), persistence));
  if (type == "contribution")
    return std::unique_ptr<Mutation>(new ContributionMutation(
      Rcpp::as<double>(spec["rate_weight"]), Rcpp::as<double>(spec["rate_contribution"]),
      Rcpp::as<double>(spec["sd"]), Rcpp::as<arma::uword>(spec["weight_loci"]), persistence));
  Rcpp::stop("unknown mutation type: " + type);
}

std::unique_ptr<Recombination> make_recombination(const Rcpp::List& spec) {
  std::string type = Rcpp::as<std::string>(spec["type"]);
  if (type == "free") return std::unique_ptr<Recombination>(new FreeRecombination());
  if (type == "map")
    return std::unique_ptr<Recombination>(new LinkageMap(Rcpp::as<arma::vec>(spec["rates"])));
  Rcpp::stop("unknown recombination type: " + type);
}

}  // namespace evomod
