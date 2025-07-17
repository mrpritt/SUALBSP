#include "solution.hpp"

#include <algorithm>
using namespace std;
using namespace sualbsp;

bool is_feasible(const Instance &I, const std::vector<Task> &π) {
  if (π.size() != I.n)
    return false;
  for (unsigned j = 1; j != I.n; ++j)
    for (unsigned i = 0; i != j; ++i)
      if (I.D[π[j]][π[i]])
        return false;
  return true;
}

void Solution::read(std::istream &in) {
  pi.clear();

  string line;
  getline(in, line);
  istringstream iss(line);
  Task t;
  while (iss >> t)
    pi.push_back(t - 1);
}

void Solution::reverse(const Instance &I) {
  for (unsigned i = 0; i != I.n; ++i)
    pi[i] = I.n - 1 - pi[i];
  ::reverse(pi.begin(), pi.end());
  assign(I, pi, a);
}

ExtendedSolution incumbent;
