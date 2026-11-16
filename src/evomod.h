// Core components of the evomod2 simulation engine.
//
// A model is a combination of components, each an abstract interface:
//
//   Architecture  genotype -> genotypic value map
//   Mutation      changes allelic effects
//   Fitness       phenotype -> fitness (log scale)
//
// To add a component, derive from the interface, implement its virtual
// functions, and add a case to the matching make_*() factory (components.cpp)
// so it can be built from the list that the R constructor returns.
//
// Genomes are stored as matrices of allelic effect vectors: traits x alleles,
// with alleles 2l and 2l + 1 (0-based) belonging to locus l.

#ifndef EVOMOD2_H
#define EVOMOD2_H

#include <RcppArmadillo.h>
#include <memory>

namespace evomod {

// ---- Genotype -> genotypic value ------------------------------------------

class Architecture {
public:
  virtual ~Architecture() {}
  // genotypic value of a genome (traits x alleles)
  virtual arma::vec genotypic_value(const arma::mat& alleles) const = 0;
  // effect of each locus in the background of the rest of the genome
  // (traits x loci): change in genotypic value when the locus is set to the
  // reference genotype (all effects zero)
  virtual arma::mat local_effects(const arma::mat& alleles) const = 0;
};

// Locus effects: sum of the two allele vectors of each locus (traits x loci).
arma::mat locus_effects(const arma::mat& alleles);

// x = sum_i Y_i
class Additive : public Architecture {
public:
  arma::vec genotypic_value(const arma::mat& alleles) const override;
  arma::mat local_effects(const arma::mat& alleles) const override;
};

// Pairwise multilinear epistasis (Hansen & Wagner 2001):
// x = S + sum_{i<j} eps(Y_i, Y_j), with eps(a, b)[k] = a' E_k b and
// E_k = epsilon.slice(k) symmetric. Computed as
// S + (eps(S, S) - sum_i eps(Y_i, Y_i)) / 2.
class Multilinear : public Architecture {
public:
  explicit Multilinear(const arma::cube& epsilon) : epsilon_(epsilon) {}
  arma::vec genotypic_value(const arma::mat& alleles) const override;
  arma::mat local_effects(const arma::mat& alleles) const override;
private:
  arma::cube epsilon_;  // slice k is E_k (traits x traits)
  arma::vec bilinear(const arma::vec& a, const arma::vec& b) const;
};

// ---- Mutation -------------------------------------------------------------

class Mutation {
public:
  virtual ~Mutation() {}
  // mutate one genome in place
  virtual void mutate(arma::mat& alleles) const = 0;
};

// Each allele mutates with probability `rate`; a mutation adds a draw from
// MvNormal(0, cov) to its effect vector.
class GaussianMutation : public Mutation {
public:
  GaussianMutation(double rate, const arma::mat& cov);
  void mutate(arma::mat& alleles) const override;
private:
  double rate_;
  arma::mat chol_lower_;  // cov = L L'
};

// ---- Fitness --------------------------------------------------------------

class Fitness {
public:
  virtual ~Fitness() {}
  // log fitness of each individual (columns of z)
  virtual arma::vec log_fitness(const arma::mat& z) const = 0;
};

// Gaussian surface: log w = -1/2 (z - theta)' omega^-1 (z - theta)
class GaussianFitness : public Fitness {
public:
  GaussianFitness(const arma::vec& theta, const arma::mat& omega);
  arma::vec log_fitness(const arma::mat& z) const override;
private:
  arma::vec theta_;
  arma::mat omega_inv_;
};

// ---- Factories: R lists -> components -------------------------------------

std::unique_ptr<Architecture> make_architecture(const Rcpp::List& spec);
std::unique_ptr<Mutation> make_mutation(const Rcpp::List& spec);

// ---- Random numbers (R's generator, so set.seed() applies) ----------------

arma::vec rnorm_vec(arma::uword n);

}  // namespace evomod

#endif
