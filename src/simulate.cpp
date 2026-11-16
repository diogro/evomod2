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

struct Records {
  std::vector<int> generation;
  std::vector<arma::vec> theta, mean_x, mean_z;
  std::vector<arma::mat> G, P;

  void add(int gen, const arma::vec& th, const arma::mat& X, const arma::mat& Z) {
    generation.push_back(gen);
    theta.push_back(th);
    mean_x.push_back(arma::mean(X, 1));
    mean_z.push_back(arma::mean(Z, 1));
    G.push_back(arma::cov(X.t()));
    P.push_back(arma::cov(Z.t()));
  }

  Rcpp::List to_list(arma::uword p) const {
    arma::uword r = generation.size();
    arma::mat th(p, r), mx(p, r), mz(p, r);
    arma::cube g(p, p, r), pp(p, p, r);
    for (arma::uword i = 0; i < r; ++i) {
      th.col(i) = theta[i]; mx.col(i) = mean_x[i]; mz.col(i) = mean_z[i];
      g.slice(i) = G[i]; pp.slice(i) = P[i];
    }
    return Rcpp::List::create(
      Rcpp::Named("generation") = generation, Rcpp::Named("theta") = th,
      Rcpp::Named("mean_x") = mx, Rcpp::Named("mean_z") = mz,
      Rcpp::Named("G") = g, Rcpp::Named("P") = pp);
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
                      Rcpp::List mutation, double env_sd, Rcpp::List regimes,
                      int generation, int record_every, bool record_start) {
  std::unique_ptr<Architecture> arch = make_architecture(architecture);
  std::unique_ptr<Mutation> mut = make_mutation(mutation);

  arma::uword p = traits, n_alleles = alleles.n_cols, n = alleles.n_slices;
  arma::uword loci = n_alleles / 2;
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
  if (record_start) rec.add(generation, theta, X, Z);

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
      for (arma::uword i = 0; i < n; ++i) {
        const arma::mat& s = alleles.slice(sires[i]);
        const arma::mat& d = alleles.slice(dames[i]);
        arma::mat& kid = offspring.slice(i);
        for (arma::uword l = 0; l < loci; ++l) {
          kid.col(2 * l)     = s.col(2 * l + (R::unif_rand() < 0.5 ? 0 : 1));
          kid.col(2 * l + 1) = d.col(2 * l + (R::unif_rand() < 0.5 ? 0 : 1));
        }
      }
      alleles.swap(offspring);
      phenotypes(*arch, env_sd, alleles, X, Z);
      ++generation;
      if (generation % record_every == 0) rec.add(generation, theta, X, Z);
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
