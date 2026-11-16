// Generation loop and the C++ functions called from R.
#include "evomod.h"

using namespace evomod;

namespace {

// Indices drawn with replacement, probability proportional to `w`.
arma::uvec sample_weighted(const arma::vec& w, arma::uword n) {
  arma::vec cum = arma::cumsum(w);
  double total = cum[cum.n_elem - 1];
  arma::uvec out(n);
  for (arma::uword i = 0; i < n; ++i) {
    double u = R::unif_rand() * total;
    const double* pos = std::upper_bound(cum.begin(), cum.end(), u);
    arma::uword idx = static_cast<arma::uword>(pos - cum.begin());
    out[i] = std::min(idx, cum.n_elem - 1);
  }
  return out;
}

// Genic part of G: 2 x the sum over loci of the covariance of allele effects
// (pooling both alleles of every individual), the G matrix expected without
// linkage disequilibrium or departures from Hardy-Weinberg. Population moments,
// scaled by n / (n - 1) like arma::cov. Empty if the architecture is not a
// sum of allele effects.
arma::mat genic_covariance(const Architecture& arch, const arma::cube& alleles) {
  arma::uword n = alleles.n_slices, loci = alleles.n_cols / 2;
  arma::mat E0 = arch.allele_effects(alleles.slice(0));
  if (E0.n_elem == 0) return arma::mat();
  arma::uword p = E0.n_rows;
  arma::mat sums(p, loci, arma::fill::zeros);
  arma::cube cross(p, p, loci, arma::fill::zeros);
  for (arma::uword i = 0; i < n; ++i) {
    arma::mat E = i == 0 ? E0 : arch.allele_effects(alleles.slice(i));
    for (arma::uword l = 0; l < loci; ++l) {
      arma::mat a = E.cols(2 * l, 2 * l + 1);
      sums.col(l) += arma::sum(a, 1);
      cross.slice(l) += a * a.t();
    }
  }
  arma::mat genic(p, p, arma::fill::zeros);
  double m = 2.0 * n;
  for (arma::uword l = 0; l < loci; ++l) {
    arma::vec mu = sums.col(l) / m;
    genic += 2.0 * (cross.slice(l) / m - mu * mu.t());
  }
  return genic * (static_cast<double>(n) / (n - 1));
}

struct Records {
  std::vector<int> generation;
  std::vector<arma::vec> theta, mean_x, mean_z;
  std::vector<arma::mat> G, P, G_genic;

  void add(int gen, const arma::vec& th, const arma::mat& X, const arma::mat& Z,
           const Architecture& arch, const arma::cube& alleles) {
    generation.push_back(gen);
    theta.push_back(th);
    mean_x.push_back(arma::mean(X, 1));
    mean_z.push_back(arma::mean(Z, 1));
    G.push_back(arma::cov(X.t()));
    P.push_back(arma::cov(Z.t()));
    G_genic.push_back(genic_covariance(arch, alleles));
  }

  Rcpp::List to_list(arma::uword p) const {
    arma::uword r = generation.size();
    arma::mat th(p, r), mx(p, r), mz(p, r);
    arma::cube g(p, p, r), pp(p, p, r);
    bool has_genic = r > 0 && G_genic[0].n_elem > 0;
    arma::cube gg(has_genic ? p : 0, has_genic ? p : 0, has_genic ? r : 0);
    for (arma::uword i = 0; i < r; ++i) {
      th.col(i) = theta[i]; mx.col(i) = mean_x[i]; mz.col(i) = mean_z[i];
      g.slice(i) = G[i]; pp.slice(i) = P[i];
      if (has_genic) gg.slice(i) = G_genic[i];
    }
    Rcpp::List out = Rcpp::List::create(
      Rcpp::Named("generation") = generation, Rcpp::Named("theta") = th,
      Rcpp::Named("mean_x") = mx, Rcpp::Named("mean_z") = mz,
      Rcpp::Named("G") = g, Rcpp::Named("P") = pp);
    if (has_genic) out["G_genic"] = gg;
    return out;
  }
};

void phenotypes(const Architecture& arch, double env_sd, const arma::cube& alleles,
                arma::mat& X, arma::mat& Z) {
  for (arma::uword i = 0; i < alleles.n_slices; ++i) {
    X.col(i) = arch.genotypic_value(alleles.slice(i));
    Z.col(i) = X.col(i) + env_sd * rnorm_vec(X.n_rows);
  }
}

}  // namespace

// Run a population through a sequence of regimes.
//
// alleles: genome rows x alleles x individuals. traits: number of traits.
// z: current phenotypes (traits x individuals), or a 0 x 0 matrix to draw them. Each regime is a list with
// generations, selective (bool), theta (optimum at the start), delta_theta
// (added to theta before selection in every generation) and omega.
// Moments are recorded after every generation that is a multiple of
// record_every, and for the starting state if record_start is true.
// [[Rcpp::export]]
Rcpp::List evolve_cpp(arma::cube alleles, int traits, arma::mat z, Rcpp::List architecture,
                      Rcpp::List mutation, Rcpp::List recombination, double env_sd,
                      Rcpp::List regimes,
                      int generation, int record_every, bool record_start) {
  std::unique_ptr<Architecture> arch = make_architecture(architecture);
  std::unique_ptr<Mutation> mut = make_mutation(mutation);
  std::unique_ptr<Recombination> recomb = make_recombination(recombination);

  arma::uword p = traits, n_alleles = alleles.n_cols, n = alleles.n_slices;
  arma::mat X(p, n), Z(p, n);
  for (arma::uword i = 0; i < n; ++i) X.col(i) = arch->genotypic_value(alleles.slice(i));
  if (z.n_elem == 0) {
    Z = X;
    for (arma::uword i = 0; i < n; ++i) Z.col(i) += env_sd * rnorm_vec(p);
  } else {
    Z = z;
  }

  Records rec;
  arma::vec theta = Rcpp::as<arma::vec>(Rcpp::as<Rcpp::List>(regimes[0])["theta"]);
  if (record_start) rec.add(generation, theta, X, Z, *arch, alleles);

  arma::cube offspring(alleles.n_rows, n_alleles, n);
  arma::vec w(n);

  for (R_xlen_t r = 0; r < regimes.size(); ++r) {
    Rcpp::List regime = regimes[r];
    int n_gen = Rcpp::as<int>(regime["generations"]);
    bool selective = Rcpp::as<bool>(regime["selective"]);
    theta = Rcpp::as<arma::vec>(regime["theta"]);
    arma::vec delta = Rcpp::as<arma::vec>(regime["delta_theta"]);
    arma::mat omega = Rcpp::as<arma::mat>(regime["omega"]);

    for (int g = 0; g < n_gen; ++g) {
      Rcpp::checkUserInterrupt();
      theta += delta;
      if (selective) {
        arma::vec logw = GaussianFitness(theta, omega).log_fitness(Z);
        w = arma::exp(logw - logw.max());
      } else {
        w.ones();
      }
      for (arma::uword i = 0; i < n; ++i) mut->mutate(alleles.slice(i));
      arma::uvec sires = sample_weighted(w, n), dames = sample_weighted(w, n);
      for (arma::uword i = 0; i < n; ++i)
        recomb->offspring(alleles.slice(sires[i]), alleles.slice(dames[i]), offspring.slice(i));
      alleles.swap(offspring);
      phenotypes(*arch, env_sd, alleles, X, Z);
      ++generation;
      if (generation % record_every == 0) rec.add(generation, theta, X, Z, *arch, alleles);
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("alleles") = alleles, Rcpp::Named("x") = X, Rcpp::Named("z") = Z,
    Rcpp::Named("theta") = theta, Rcpp::Named("generation") = generation,
    Rcpp::Named("records") = rec.to_list(p));
}

// [[Rcpp::export]]
arma::mat genotypic_values_cpp(arma::cube alleles, int traits, Rcpp::List architecture) {
  std::unique_ptr<Architecture> arch = make_architecture(architecture);
  arma::mat X(traits, alleles.n_slices);
  for (arma::uword i = 0; i < alleles.n_slices; ++i)
    X.col(i) = arch->genotypic_value(alleles.slice(i));
  return X;
}

// [[Rcpp::export]]
arma::mat genic_covariance_cpp(arma::cube alleles, Rcpp::List architecture) {
  return genic_covariance(*make_architecture(architecture), alleles);
}

// [[Rcpp::export]]
arma::mat local_effects_cpp(arma::mat alleles, Rcpp::List architecture) {
  return make_architecture(architecture)->local_effects(alleles);
}

// [[Rcpp::export]]
arma::mat mutate_cpp(arma::mat alleles, Rcpp::List mutation) {
  make_mutation(mutation)->mutate(alleles);
  return alleles;
}

// [[Rcpp::export]]
arma::vec log_fitness_cpp(arma::mat z, arma::vec theta, arma::mat omega) {
  return GaussianFitness(theta, omega).log_fitness(z);
}
