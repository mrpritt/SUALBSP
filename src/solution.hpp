#pragma once

#include <iostream>
#include <vector>

#include "dimensions.hpp"
#include "instance.hpp"
#include "logging.hpp"
#include "lowerbounds.hpp"
#include "prbp.hpp"
#include "util.hpp"

struct Solution {
  std::vector<Task> pi;
  std::vector<Task> a;
  unsigned m;

  Solution() : m(0) {}
  Solution(std::vector<Task> pi, unsigned m) : pi(pi), m(m) {}
  bool isValid() const { return pi.size() > 0; }

  void read(std::istream&);
  void reverse(const sualbsp::Instance&);
};

struct ExtendedSolution : public Solution {
  double time;
  unsigned iter;
  short unsigned lb;
  ExtendedSolution() : time(logging::elapsed()), iter(0), lb(0) {}

  bool set(sualbsp::Instance &I, const Solution &s) {
    if (s.m < lb)
      throw "Error: solution less than lower bound";

    if (s.m < m) {
      Solution::operator=(s);
      time = logging::elapsed();
      iter = S.iter;
      new_incumbent(I);
      assign(I, pi, a);
      return true;
    }
    return false;
  }

  bool set(sualbsp::Instance &I, const std::vector<unsigned short int> &assignment, const std::vector<unsigned short int> &solution, unsigned short ms) {
    if (ms < m) {
      iter = S.iter;
      time = logging::elapsed();
      m = ms;
      pi = assignment;
      a = solution;
      new_incumbent(I);
      return true;
    }
    return false;
  }
};

extern ExtendedSolution incumbent;

bool is_feasible(const sualbsp::Instance &, const std::vector<Task> &);
