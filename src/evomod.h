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

// Contribution matrix (do O & Whitlock 2023): k weight genes a_j and a
// traits x k matrix of contribution genes c_ij, each a locus holding one scalar
// per allele (the two alleles of a locus are summed). Loci 0..k-1 are the
// weight genes; locus k + j * p + i is c_ij. x_i = sum_j c_ij a_j.
class ContributionMatrix : public Architecture {
public:
  ContributionMatrix(arma::uword traits, arma::uword genes) : p_(traits), k_(genes) {}
  arma::vec genotypic_value(const arma::mat& alleles) const override;
  arma::mat local_effects(const arma::mat& alleles) const override;
private:
  arma::uword p_, k_;
  void unpack(const arma::mat& alleles, arma::vec& a, arma::mat& C) const;
};

// ---- Mutation -------------------------------------------------------------

class Mutation {
public:
  virtual ~Mutation() {}
  // mutate one genome in place (one generation of mutation)
  virtual void mutate(arma::mat& alleles) const = 0;
  // apply exactly one mutation event, at a site chosen in proportion to the
  // mutation rates
  virtual void mutate_one(arma::mat& alleles) const = 0;
  // expected number of mutation events per genome per generation
  virtual double expected_mutations(const arma::mat& alleles) const = 0;
};

// Each allele mutates with probability `rate`; a mutation replaces its effect
// vector e by persistence * e + MvNormal(0, cov): persistence = 1 is a random
// walk (effects accumulate), 0 the house of cards model (a fresh draw,
// unrelated to the old effect; Kingman 1978, Turelli 1984). The same holds for
// the other mutation models below.
class GaussianMutation : public Mutation {
public:
  GaussianMutation(double rate, const arma::mat& cov, double persistence = 1.0);
  void mutate(arma::mat& alleles) const override;
  void mutate_one(arma::mat& alleles) const override;
  double expected_mutations(const arma::mat& alleles) const override {
    return rate_ * alleles.n_cols;
  }
private:
  double rate_;
  arma::mat chol_lower_;  // cov = L L'
  double persistence_;
};

// Mutation for BinaryB genomes. Each allele's y mutates with probability
// rate_y, adding Normal(0, sd_y); each entry of B mutates with probability
// rate_b. A mutated entry flips; with random_state it gets a random value
// instead (0 or 1, so it changes half of the time, as in the C code used for
// the published simulations).
class BinaryBMutation : public Mutation {
public:
  BinaryBMutation(double rate_y, double sd_y, double rate_b, bool random_state,
                  double persistence = 1.0)
    : rate_y_(rate_y), sd_y_(sd_y), rate_b_(rate_b), random_state_(random_state),
      persistence_(persistence) {}
  void mutate(arma::mat& alleles) const override;
  void mutate_one(arma::mat& alleles) const override;
  double expected_mutations(const arma::mat& alleles) const override {
    return rate_y_ * alleles.n_cols + rate_b_ * alleles.n_cols * (alleles.n_rows - 1);
  }
private:
  double rate_y_, sd_y_, rate_b_;
  bool random_state_;
  double persistence_;
};

// Mutation for contribution matrix genomes: alleles of the first
// `weight_loci` loci mutate with probability rate_weight, the others with
// rate_contribution; a mutation adds Normal(0, sd) to the allele's value.
class ContributionMutation : public Mutation {
public:
  ContributionMutation(double rate_weight, double rate_contribution, double sd,
                       arma::uword weight_loci, double persistence = 1.0)
    : rate_w_(rate_weight), rate_c_(rate_contribution), sd_(sd), weight_loci_(weight_loci),
      persistence_(persistence) {}
  void mutate(arma::mat& alleles) const override;
  void mutate_one(arma::mat& alleles) const override;
  double expected_mutations(const arma::mat& alleles) const override {
    return rate_w_ * 2 * weight_loci_ + rate_c_ * (alleles.n_cols - 2 * weight_loci_);
  }
private:
  double rate_w_, rate_c_, sd_;
  arma::uword weight_loci_;
  double persistence_;
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
