// Core components of the evomod2 simulation engine.
//
// A model is a combination of components, each an abstract interface:
//
//   Architecture  genotype -> genotypic value map
//   Mutation      changes allelic effects
//   Fitness       phenotype -> fitness (log scale)
//   Recombination parents -> offspring genome (meiosis)
//
// To add a component, derive from the interface, implement its virtual
// functions, and add a case to the matching make_*() factory (components.cpp)
// so it can be built from the list that the R constructor returns.
//
// Genomes are matrices with one column per allele (alleles 2l and 2l + 1,
// 0-based, belong to locus l). The rows hold the allele's state, as defined
// by the architecture: an effect vector (one row per trait) for Additive and
// Multilinear; the scalar effect y followed by the column of the binary
// matrix B for BinaryB.

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
  // effect vector of each allele (traits x alleles) when the genotypic value
  // is their sum; an empty matrix for architectures where it is not
  virtual arma::mat allele_effects(const arma::mat& alleles) const { return arma::mat(); }
};

// Locus effects: sum of the two allele vectors of each locus (traits x loci).
arma::mat locus_effects(const arma::mat& alleles);

// x = sum_i Y_i
class Additive : public Architecture {
public:
  arma::vec genotypic_value(const arma::mat& alleles) const override;
  arma::mat local_effects(const arma::mat& alleles) const override;
  arma::mat allele_effects(const arma::mat& alleles) const override { return alleles; }
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

// Binary ontogenetic matrix (Melo & Marroig 2015): row 0 holds the allelic
// effects y, rows 1..p the binary matrix B (traits x alleles); x = B y.
class BinaryB : public Architecture {
public:
  arma::vec genotypic_value(const arma::mat& alleles) const override;
  arma::mat local_effects(const arma::mat& alleles) const override;
  arma::mat allele_effects(const arma::mat& alleles) const override;
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

// Mutation for BinaryB genomes. Each allele's y mutates with probability
// rate_y, adding Normal(0, sd_y); each entry of B mutates with probability
// rate_b. A mutated entry flips; with random_state it gets a random value
// instead (0 or 1, so it changes half of the time, as in the C code used for
// the published simulations).
class BinaryBMutation : public Mutation {
public:
  BinaryBMutation(double rate_y, double sd_y, double rate_b, bool random_state)
    : rate_y_(rate_y), sd_y_(sd_y), rate_b_(rate_b), random_state_(random_state) {}
  void mutate(arma::mat& alleles) const override;
private:
  double rate_y_, sd_y_, rate_b_;
  bool random_state_;
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

// ---- Recombination --------------------------------------------------------

// Builds an offspring genome from two parents: column 2l of the offspring gets
// one of the sire's two alleles at locus l, column 2l + 1 one of the dam's.
class Recombination {
public:
  virtual ~Recombination() {}
  virtual void offspring(const arma::mat& sire, const arma::mat& dame, arma::mat& kid) const = 0;
};

// Every locus independent (each allele drawn with probability 1/2).
class FreeRecombination : public Recombination {
public:
  void offspring(const arma::mat& sire, const arma::mat& dame, arma::mat& kid) const override;
};

// Loci in a line: rates[l] is the probability of a crossover between locus l
// and l + 1 (0.5 between chromosomes). A gamete starts on a random haplotype
// and switches haplotype at each crossover.
class LinkageMap : public Recombination {
public:
  explicit LinkageMap(const arma::vec& rates) : rates_(rates) {}
  void offspring(const arma::mat& sire, const arma::mat& dame, arma::mat& kid) const override;
private:
  arma::vec rates_;
  void gamete(const arma::mat& parent, arma::mat& kid, arma::uword slot) const;
};

// ---- Factories: R lists -> components -------------------------------------

std::unique_ptr<Architecture> make_architecture(const Rcpp::List& spec);
std::unique_ptr<Mutation> make_mutation(const Rcpp::List& spec);
std::unique_ptr<Recombination> make_recombination(const Rcpp::List& spec);

// ---- Random numbers (R's generator, so set.seed() applies) ----------------

arma::vec rnorm_vec(arma::uword n);

}  // namespace evomod

#endif
