#pragma once

#include <functional>
#include <vector>

#include <boost/multi_array.hpp>

#include "bitvector.hpp"

namespace sualbsp {
using namespace std;
using namespace boost;

struct Instance {
  unsigned n;
  unsigned c;
  unsigned optm;

  vector<Time> t0, t;
  vector<BitRow> d;
  multi_array<Time, 2> sf0, sf;
  multi_array<Time, 2> sb0, sb;

  vector<vector<Task>> P, F;
  vector<vector<Task>> Ps, Fs, Is;
  unsigned maxPs, maxFs;
  vector<vector<vector<unsigned short>>> O;
  vector<BitRow> D;
  vector<Time> sfi;
  vector<Time> sbi;

  bool setupDirected;
  bool fwdSymmetric;
  bool bwdSymmetric;
  bool isTriangular;
  vector<int> insertionlb;
  bool isStrictTriangular;

  vector<Time> ta, tn;
  vector<unsigned short> E, T;

  Instance() : n(0) {}
  Instance(Instance &&other);
  Instance(const Instance &other);
  Instance &operator=(Instance other);
  void swap(Instance &other);

  void allocate();
  void readMP(istream &,bool=true);
  void readSBF(istream &, bool=true);
  void read(istream &, std::string,bool=true);
  void writeSBF(ostream &);

  void reverse();

  Time get_pw(Task i) const {
    Time pw = t[i];
    for (auto j : Fs[i])
      pw += t[j];
    return pw;
  }

  Time get_ω(Task i) const {
    Time ω = inf_time;
    for (unsigned j = 0; j != n; ++j) {
      if (canForwardExt(i, j))
        ω = min(ω, sf[i][j]);
      if (canBackwardExt(i, j))
        ω = min(ω, sb[i][j]);
    }
    return ω;
  }

  Time get_α(Task i) const {
    Time α = inf_time;
    for (unsigned h = 0; h != n; ++h) {
      if (canForwardExt(h, i))
        α = min(α, sf[h][i]);
      if (canBackwardExt(h, i))
        α = min(α, sb[h][i]);
    }
    return α;
  }

  Time get_ωτ(Task i, const BitVector &) const;
  Time get_αμ(Task i, const BitVector &) const;
  Time get_ωμ(Task i, const BitVector &) const;
  Time get_ατ(Task i, const BitVector &) const;

  pair<unsigned, unsigned> FS(Task i, unsigned um) const { return {E[i], um + 2 - T[i]}; }
  bool canAssign(Task i, unsigned k, unsigned um) const { return E[i] <= k && k <= um + 1 - T[i]; }
  bool canAssignFwd(Task i, Task j, unsigned k, unsigned um) const;
  bool canAssignBwd(Task i, Task j, unsigned k, unsigned um) const;
  bool canShare(Task i, Task j, unsigned um) const;

  bool canForward(Task i, Task j) const { return !((D[i][j] && !d[i][j]) || D[j][i] || i == j); }
  bool canBackward(Task i, Task j) const { return !(D[i][j]); }

  bool canForwardExt(Task i, Task j) const { return !((D[i][j] && !d[i][j]) || D[j][i] || i == j || sf[i][j] > c); }
  bool canBackwardExt(Task i, Task j) const { return !(D[i][j] || sb[i][j] > c); }

  double avgProcessingTime() const;
  double avgForwardSetupTime() const;
  double avgBackwardSetupTime() const;
  double orderStrength() const;

  void preprocess_dynamic(unsigned);

private:
  Time uc();
  void handleSetupSection(istream &, multi_array<Time, 2> &);
  void compute_ta_tn();
  void compute_E_T();
  int excludeSetups(unsigned);
  int increaseTimes();
  int isolatedTasks();

  void preprocess(bool=true);
  void preprocess_static();
  void compute_smallest_setups();
  void compute_ILB();
  void check_symmetry();

  bool taskTriangular(Task i) { return insertionlb[i] >= 0; }
  bool taskStrictTriangular(Task i) { return insertionlb[i] >= int(t[i]); }
  unsigned disableForward(Task i, Task j) {
    if (sf[i][j] <= c) {
      sf[i][j] = c + 1;
      return 1;
    }
    return 0;
  }
  unsigned disableBackward(Task i, Task j) {
    if (sb[i][j] <= c) {
      sb[i][j] = c + 1;
      return 1;
    }
    return 0;
  }
};
}

std::string guess_format(std::string);
