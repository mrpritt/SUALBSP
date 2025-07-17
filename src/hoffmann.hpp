#pragma once

#include <algorithm>
#include <vector>

#include "instance.hpp"
#include "util.hpp"

using namespace sualbsp;

struct hoffmann_t {
  std::vector<short int> deg;
  std::vector<short int> eligible;
  std::vector<short int> tasks;
  std::vector<short int> best_tasks;
  std::vector<double> weights;
  std::vector<bool> pending;

  short int n_best_tasks;
  double best;
  unsigned n_loads;

  hoffmann_t(const Instance &I);

  void set_degrees(const std::vector<std::vector<Task>> &);
  void set_weights(const std::vector<Time> &);
  double total_weight() const { return std::accumulate(weights.begin(), weights.end(), 0.0); }
  void compute_eligible(unsigned &);
  unsigned short num_assigned() const { return std::count(pending.begin(), pending.end(), false); }

  void compute_standard_weights(const Instance &);
  void compute_sw_weights(const Instance &, double, double, double);
  void compute_s_weights(const Instance &, double, double, double, bool);
  const double s_epsilon = 0.000001;
};

unsigned short int Hoffmann(Instance &I);
unsigned short int SW_Hoffmann(Instance &I, unsigned short lb);
unsigned short int FH_Hoffmann(Instance &I, unsigned short lb);
unsigned short int S_Hoffmann(Instance &I, unsigned short lb);
