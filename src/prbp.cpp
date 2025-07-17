#include "prbp.hpp"

using namespace std;

#include "options.hpp"
#include "random.hpp"
#include "solution.hpp"
#include "util.hpp"

using namespace sualbsp;

vector<Task> randomOrder(const Instance &I) {
  vector<Task> ind(I.n, 0);
  vector<Task> pi;

  vector<Task> Q;
  for (unsigned i = 0; i < I.n; i++) {
    ind[i] = I.P[i].size();
    if (ind[i] == 0)
      Q.push_back(i);
  }
  while (Q.size() > 0) {
    short idx = getRandom(0, Q.size() - 1);
    swap(Q[idx], Q.back());
    Task i = Q.back();
    Q.pop_back();
    pi.push_back(i);
    for (auto j : I.F[i]) {
      ind[j]--;
      if (ind[j] == 0)
        Q.push_back(j);
    }
  }
  return pi;
}

struct AssignContext {
  const Instance &I;
  vector<Task> pi;

  vector<Task> A;
  unsigned notFit;

  unsigned m;

  vector<Task> ind;

  vector<int> set;
  unsigned it, ft;
  unsigned iti;
  unsigned l;

  vector<unsigned> tot_tau;
  vector<unsigned> ind_tau;

  AssignContext(const Instance &I) : I(I), notFit(0), m(1), ind(I.n, 0), set(I.n, 0), l(0), tot_tau(I.n, 0), ind_tau(I.n, 0) {
    initialize_tau();
  }

  bool isFeasible() {
    for (unsigned i = 0; i < I.n; i++)
      if (I.t[i] + I.sb[i][i] > I.c)
        return false;
    return true;
  }

#define TAU_FEASIBLE false
  void initialize_tau() {

    for (unsigned i = 0; i != I.n; ++i) {
      tot_tau[i] = 0;
      ind_tau[i] = 0;

      for (auto j : I.P[i]) {
        if (I.sf[j][i] <= I.c || TAU_FEASIBLE) {
          tot_tau[i] += I.sf[j][i];
          ind_tau[i]++;
        }
      }

      for (auto j : I.Is[i]) {
        if (I.sf[j][i] <= I.c || TAU_FEASIBLE) {
          tot_tau[i] += I.sf[j][i];
          ind_tau[i]++;
        }
      }
    }
  }

  void shadow(Task i) {
    for (auto j : I.F[i]) {
      if (I.sf[i][j] <= I.c || TAU_FEASIBLE) {
        ind_tau[j]--;
        tot_tau[j] -= I.sf[i][j];
      }
    }
    for (auto j : I.Is[i]) {
      if (I.sf[i][j] <= I.c || TAU_FEASIBLE) {
        ind_tau[j]--;
        tot_tau[j] -= I.sf[i][j];
      }
    }
  }

  void unshadow(Task i) {
    for (auto j : I.F[i]) {
      if (I.sf[i][j] <= I.c || TAU_FEASIBLE) {
        ind_tau[j]++;
        tot_tau[j] += I.sf[i][j];
      }
    }
    for (auto j : I.Is[i]) {
      if (I.sf[i][j] <= I.c || TAU_FEASIBLE) {
        ind_tau[j]++;
        tot_tau[j] += I.sf[i][j];
      }
    }
  }

  void first(Task i) {
    it = ft = i;
    l = I.t[i] + set[i];
    iti = pi.size();
    if (iti > 0)
      shadow(pi.back());
    pi.push_back(i);
  }

  void next(Task i) {
    ft = i;
    l += I.t[i] + set[i];
    shadow(pi.back());
    pi.push_back(i);
  }

  void open() {
    m++;
    l = 0;
    notFit = A.size();
    for (unsigned jb = 0, je = A.size(); jb != je; jb++) {
      Task j = A[jb];
      set[j] = I.sb[j][j];
    }
  }

  void setup() {
    A.clear();
    for (unsigned i = 0; i != I.n; ++i) {
      ind[i] = I.P[i].size();
      if (ind[i] == 0) {
        set[i] = I.sb[i][i];
        A.push_back(i);
      }
    }
    notFit = A.size();
  }

  bool queue_ok() {
    for (unsigned jb = 0, je = notFit; jb != je; ++jb) {
      Task j = A[jb];
      if (l + I.t[j] + set[j] > I.c)
        return false;
    }
    for (unsigned jb = notFit, je = A.size(); jb != je; ++jb) {
      Task j = A[jb];
      if (l + I.t[j] + set[j] <= I.c)
        return false;
    }
    return true;
  }

  Task remove(unsigned idx) {
    if (idx < notFit) {
      swap(A[idx], A[notFit - 1]);
      swap(A[notFit - 1], A.back());
      notFit--;
    } else {
      swap(A[idx], A.back());
    }
    Task i = A.back();
    A.pop_back();

    for (unsigned jb = 0, je = A.size(); jb != je; jb++) {
      Task j = A[jb];
      set[j] = int(-I.sb[ft][it]) + int(I.sf[ft][j]) + int(I.sb[j][it]);
    }
    unsigned jb = 0, je = A.size();
    while (jb != je) {
      Task j = A[jb];
      if (jb < notFit && l + I.t[j] + set[j] > I.c) {
        swap(A[jb], A[notFit - 1]);
        notFit--;
        continue;
      }
      if (jb >= notFit && l + I.t[j] + set[j] <= I.c) {
        swap(A[jb], A[notFit]);
        notFit++;
      }
      jb++;
    }

    for (auto j : I.F[i]) {
      ind[j]--;
      if (ind[j] == 0) {
        set[j] = int(-I.sb[ft][it]) + int(I.sf[ft][j]) + int(I.sb[j][it]);
        A.push_back(j);
        if (l + I.t[j] + set[j] <= I.c) {
          swap(A[notFit], A.back());
          notFit++;
        }
      }
    }

    return i;
  }
};

bool rule_precondition(AssignContext &ctx) { return ctx.notFit > 0; }

template <typename PrioT, typename Prio> vector<PrioT> getValues(AssignContext &ctx, bool first, bool fit, Prio p) {
  unsigned le = fit ? ctx.notFit : ctx.A.size();

  vector<PrioT> v(le);
  for (unsigned i = 0; i != le; i++)
    v[i] = p(ctx, first, ctx.A[i]);
  return v;
}

template <typename PrioT, typename Prio> vector<short> getOrder(AssignContext &ctx, bool first, Prio p) {
  vector<PrioT> v = getValues<PrioT>(ctx, first, false, p);
  vector<short> pi(v.size());
  iota(pi.begin(), pi.end(), 0);
  sort(pi.begin(), pi.end(), [&v](short i, short j) { return v[i] > v[j]; });
  return pi;
}

template <typename PrioT, typename Prio> vector<short> getRandomOrder(AssignContext &ctx, bool first, Prio p) {
  vector<PrioT> v = getValues<PrioT>(ctx, first, false, p);
  vector<short> pi(v.size());
  iota(pi.begin(), pi.end(), 0);

  for (unsigned i = 0, ie = v.size(); i != ie; i++) {
    auto d = discrete_distribution(v.begin() + i, v.end());
    short j = i + d(rng);
    swap(v[i], v[j]);
    swap(pi[i], pi[j]);
  }
  return pi;
}

template <typename PrioT, typename Prio> short getBest(AssignContext &ctx, bool first, Prio p) {
  vector<PrioT> v = getValues<PrioT, Prio>(ctx, first, true, p);
  short best = 0;
  for (unsigned long i = 0u, ie = v.size(); i != ie; ++i) {
    const auto j = ctx.A[i], bj = ctx.A[best];
    if ((v[i] > v[best]) || ((v[i] == v[best]) && ((int(ctx.I.t[j]) - ctx.set[j] > int(ctx.I.t[bj]) - ctx.set[bj]) || (int(ctx.I.t[j]) - ctx.set[j] == int(ctx.I.t[bj]) - ctx.set[bj] && j < bj))))
      best = i;
  }
  return best;
}

namespace MaxTMinSU {
typedef int PrioT;
const bool stochastic = false;

PrioT value(const AssignContext &ctx, bool first, Task j) { return int(ctx.I.t[j]) - ctx.set[j]; }

short best(AssignContext &ctx, bool first) { return getBest<PrioT>(ctx, first, value); }

vector<short> order(AssignContext &ctx, bool first) { return getOrder<PrioT>(ctx, first, value); }
};

struct MtSlack {
  int m, s;
  bool operator<(const MtSlack &other) const { return m * other.s < other.m * s ; }
  bool operator>(const MtSlack &other) const { return m * other.s > other.m * s ; }
  bool operator==(const MtSlack &other) const { return m * other.s == other.m * s; }
};

namespace MaxTMinSUSlack {
typedef MtSlack PrioT;
const bool stochastic = false;

PrioT value(const AssignContext &ctx, bool first, Task j) {
  const int slack = int(incumbent.m + 1) - ctx.I.T[j] - ctx.I.E[j] + 1;
  return {int(ctx.I.t[j]) - ctx.set[j], slack};
}

short best(AssignContext &ctx, bool first) { return getBest<PrioT>(ctx, first, value); }

vector<short> order(AssignContext &ctx, bool first) { return getOrder<PrioT>(ctx, first, value); }
};

namespace MaxF {
typedef unsigned PrioT;
const bool stochastic = false;

PrioT value(const AssignContext &ctx, bool first, Task j) { return ctx.I.F[j].size(); }

short best(AssignContext &ctx, bool first) { return getBest<PrioT>(ctx, first, value); }

vector<short> order(AssignContext &ctx, bool first) { return getOrder<PrioT>(ctx, first, value); }
};

namespace MinSlack {
typedef unsigned PrioT;
const bool stochastic = false;

PrioT value(AssignContext &ctx, bool first, Task j) { return ctx.I.T[j]; }

short best(AssignContext &ctx, bool first) { return getBest<PrioT>(ctx, first, value); }

vector<short> order(AssignContext &ctx, bool first) { return getOrder<PrioT>(ctx, first, value); }
};

namespace Random {
bool stochastic = true;

short best(AssignContext &ctx, bool first) { return getRandom(0, ctx.notFit - 1); }

vector<short> order(AssignContext &ctx, bool first) {
  vector<short> pi(ctx.A.size());
  iota(pi.begin(), pi.end(), 0);
  shuffle(pi.begin(), pi.end(), sualbsp::rng);
  return pi;
}
};

namespace MPRule {
typedef double PrioT;
bool stochastic = false;

PrioT value(AssignContext &ctx, bool first, Task j) {
  PrioT tv = 5.0 * ctx.I.t[j] + 0.3 * ctx.I.tn[j] - 45.3 * (first ? 0 : ctx.I.sf[ctx.ft][j]);
  for (auto s : ctx.I.Fs[j]) {
    if (ctx.ind_tau[s] > 0)
      tv += 0.3 * 3.9 * PrioT(ctx.tot_tau[s]) / PrioT(ctx.ind_tau[s]);
  }
  return tv;
}

short best(AssignContext &ctx, bool first) { return getBest<PrioT>(ctx, first, value); }

vector<short> order(AssignContext &ctx, bool first) { return getOrder<PrioT>(ctx, first, value); }
};

namespace Fs {
typedef unsigned PrioT;
bool stochastic = false;

PrioT value(AssignContext &ctx, bool first, Task j) { return ctx.I.Fs[j].size() + 1; }

short best(AssignContext &ctx, bool first) {
  vector<PrioT> Fs = getValues<PrioT>(ctx, first, true, value);
  return discrete_distribution(Fs.begin(), Fs.end())(rng);
}

vector<short> order(AssignContext &ctx, bool first) { return getRandomOrder<PrioT>(ctx, first, value); }
};

namespace availableSpacesQuotient {
typedef double PrioT;
bool stochastic = true;

PrioT value(AssignContext &ctx, bool first, Task j) {
  unsigned free = ctx.I.n - ctx.pi.size();
  return PrioT(free + ctx.I.Fs[j].size()) / PrioT(free - ctx.I.Fs[j].size());
}

short best(AssignContext &ctx, bool first) {
  vector<PrioT> asq = getValues<PrioT>(ctx, first, true, value);
  return discrete_distribution(asq.begin(), asq.end())(rng);
}

vector<short> order(AssignContext &ctx, bool first) { return getRandomOrder<PrioT>(ctx, first, value); }
};

namespace antiASQ {
typedef double PrioT;
bool stochastic = true;

PrioT value(AssignContext &ctx, bool first, Task j) {
  unsigned free = ctx.I.n - ctx.pi.size();
  return PrioT(free - ctx.I.Fs[j].size()) / PrioT(free + ctx.I.Fs[j].size());
}

short best(AssignContext &ctx, bool first) {
  vector<PrioT> asq = getValues<PrioT>(ctx, first, true, value);
  return discrete_distribution(asq.begin(), asq.end())(rng);
}

vector<short> order(AssignContext &ctx, bool first) { return getRandomOrder<PrioT>(ctx, first, value); }
};

typedef short (*rule)(AssignContext &, bool);
typedef vector<short> (*ruleOrder)(AssignContext &, bool);

using ruleDescription = tuple<string,rule,ruleOrder,bool>;


vector<ruleDescription> sbf_rules {
  { "mprule",          MPRule::best,         MPRule::order,         MPRule::stochastic,           },
  { "maxtminsu",       MaxTMinSU::best,      MaxTMinSU::order,      MaxTMinSU::stochastic,        },
  { "maxtminsuslack",  MaxTMinSUSlack::best, MaxTMinSUSlack::order, MaxTMinSUSlack::stochastic,   },
  { "maxf",            MaxF::best,           MaxF::order,           MaxF::stochastic,             },
  { "minslack",        MinSlack::best,       MinSlack::order,       MinSlack::stochastic,         },
  { "random",          Random::best,         Random::order,         Random::stochastic,           },
};

namespace sbf {
bool stochastic = true;

short best(AssignContext &ctx, bool first) {
  auto rule = sbf_rules.begin();
  std::advance(rule, getRandom(0, sbf_rules.size() - 1));
  return get<1>(*rule)(ctx, first);
}

vector<short> order(AssignContext &ctx, bool first) {
  auto rule = sbf_rules.begin();
  std::advance(rule, getRandom(0, sbf_rules.size() - 1));
  return get<2>(*rule)(ctx, first);
}

namespace rule_grasp {
bool stochastic = true;

vector<ruleDescription> rules;

short best(AssignContext &ctx, bool first) {
  auto rule = rules.begin();
  std::advance(rule, getRandom(0, rules.size() - 1));
  return get<1>(*rule)(ctx, first);
}

vector<short> order(AssignContext &ctx, bool first) {
  auto rule = rules.begin();
  std::advance(rule, getRandom(0, rules.size() - 1));
  return get<2>(*rule)(ctx, first);
}

};

};

map<string,tuple<rule,ruleOrder,bool>> rules {
  { "maxtminsu",      { MaxTMinSU::best,               MaxTMinSU::order,               MaxTMinSU::stochastic,               } },
  { "maxtminsuslack", { MaxTMinSUSlack::best,          MaxTMinSUSlack::order,          MaxTMinSUSlack::stochastic,          } },
  { "maxf",           { MaxF::best,                    MaxF::order,                    MaxF::stochastic,                    } },
  { "minslack",       { MinSlack::best,                MinSlack::order,                MinSlack::stochastic,                } },
  { "random",         { Random::best,                  Random::order,                  Random::stochastic,                  } },
  { "mprule",         { MPRule::best,                  MPRule::order,                  MPRule::stochastic,                  } },
  { "randomfs",       { Fs::best,                      Fs::order,                      Fs::stochastic                       } },
  { "sbf",            { sbf::best,                     sbf::order,                     sbf::stochastic                      } },
  { "asq",            { availableSpacesQuotient::best, availableSpacesQuotient::order, availableSpacesQuotient::stochastic  } },
  { "anti-asq",       { antiASQ::best,                 antiASQ::order,                 antiASQ::stochastic                  } },
};

short nextStation(const Instance &I, const vector<Task> &π, short i) {
  unsigned k = i, l = I.t[π[i]], j = i + 1;

  while (j < I.n) {
    l += I.sf[π[j - 1]][π[j]] + I.t[π[j]];
    if (l + I.sb[π[j]][π[i]] <= I.c)
      k = j;
    j++;
  }

  l = 0;
  return k + 1;
}

unsigned short assign(const Instance &I, const vector<Task> &π, vector<short unsigned> &a) {
  a.assign(I.n, 0);

  unsigned j = 0, m = 1;
  while (j != I.n) {
    unsigned j_ = nextStation(I, π, j);
    for (auto i = j; i != j_; ++i)
      a[i] = m;
    m++;
    j = j_;
  }

  return m - 1;
}

unsigned short assign(const Instance &I, const vector<Task> &π) {
  vector<short unsigned> a;
  return assign(I, π, a);
}

const vector<Task> hoffmann_sequence{0, 2, 3, 1, 4, 5, 6, 7, 8, 20, 12, 13, 9, 10, 11, 14, 15, 17, 18, 16, 19};

unsigned short sample(Instance &I, unsigned short lb) {
  Solution s;
  while (logging::elapsed() < opt.tlim && S.iter++ < opt.ilim && s.m > lb) {
    s.pi = randomOrder(I);
    s.m = assign(I, s.pi);
    incumbent.set(I, s);
  }
  return incumbent.m;
}

template <typename Rule> Solution assignOnce(const Instance &I, unsigned short lb, Rule rule) {

  AssignContext ctx(I);

  ctx.setup();

  unsigned ti = rule(ctx, true);
  Task t = ctx.A[ti];
  ctx.first(t);
  ctx.remove(ti);

  while (ctx.A.size() > 0) {
    bool first = false;
    if (ctx.notFit == 0) {
      ctx.open();
      first = true;
      if (ctx.m >= incumbent.m) {
        ctx.pi.clear();
        return Solution(ctx.pi, uint_inf);
      }
    }
    ti = rule(ctx, first);
    t = ctx.A[ti];
    if (first)
      ctx.first(t);
    else
      ctx.next(t);
    ctx.remove(ti);
  }

  return Solution{ctx.pi, ctx.m};
}

vector<Task> greedyOrder(const Instance &I, const vector<Task> &L) {
  vector<Task> G(L.begin(), L.end());
  const auto n0 = L.size();

  const auto n = I.n;
  vector<unsigned> ind(n, 0);
  for (auto i : L)
    for (auto j : I.F[i])
      ind[j]++;

  unsigned k = 0;
  while (ind[G[k]] != 0)
    k++;
  for (auto l = k + 1; l != n0; ++l) {
    const Task tk = G[k];
    const Task i = G[l];
    if (ind[i] == 0 && I.sb[i][i] < I.sb[tk][tk])
      k = l;
  }
  swap(G[0], G[k]);

  unsigned it = G[0];
  unsigned ft = G[0];
  for (auto j : I.F[ft])
    ind[j]--;

  for (unsigned j = 1; j != n0; j++) {
    unsigned k = j;
    while (ind[G[k]] != 0)
      k++;
    for (auto l = k + 1; l != n0; ++l) {
      Task tk = G[k];
      Task i = G[l];
      if (ind[i] == 0 && (int(I.sf[ft][i]) + int(I.sb[i][it])) < (int(I.sf[ft][tk]) + int(I.sb[tk][it])))
        k = l;
    }
    swap(G[j], G[k]);
    ft = G[j];
    for (auto j : I.F[ft])
      ind[j]--;
  }

  return G;
}

Time assignLoad(const Instance &I, const vector<Task> &L, const vector<short> &π) {
  const unsigned n = L.size();
  if (n == 0)
    return 0;

  Time load = I.t[L[π[0]]];
  for (unsigned i = 1; i != n; ++i) {
    load += I.sf[L[π[i - 1]]][L[π[i]]] + I.t[L[π[i]]];
    if (load > I.c)
      return load;
  }
  return load + I.sb[L[π[n - 1]]][L[π[0]]];
}

pair<vector<Task>, Time> reOrder(const Instance &I, const vector<Task> &L) {
  const unsigned n = L.size();

  vector<Task> G = greedyOrder(I, L);
  vector<short> π(n);
  iota(π.begin(), π.end(), 0);

  unsigned reopt_iter = 0;
  do {
    reopt_iter++;
    Time load = assignLoad(I, L, π);

    if (load <= I.c) {
      vector<Task> H(n);
      for (auto i = 0u; i != n; ++i)
        H[i] = G[π[i]];
      return {H, load};
    }
  } while (next_permutation(π.begin(), π.end()));
  return {L, I.c + 1};
}

pair<vector<Task>, Time> reOrder_prec(const Instance &I, const vector<Task> &L) {
  const unsigned n = L.size();

  vector<Task> G = greedyOrder(I, L);
  vector<short> π(n), ι(n);
  iota(π.begin(), π.end(), 0);
  iota(ι.begin(), ι.end(), 0);

  unsigned reopt_iter = 0;
  do {
    reopt_iter++;
    Time load = assignLoad(I, G, π);

    if (load <= I.c) {
      vector<Task> H(n);
      for (auto i = 0u; i != n; ++i)
        H[i] = G[π[i]];
      return {H, load};
    }

    short k = n - 1;

    while (k >= 0) {
      auto j = ι[k];
      auto l = j == 0 ? 0 : π[j - 1];
      if (j != 0 && !I.D[G[l]][G[k]]) {
        π[j - 1] = k;
        π[j] = l;
        ι[k] = j - 1;
        ι[l] = j;
        goto moved;
      }
      while (j < k) {
        l = π[j + 1];
        π[j] = l;
        ι[l] = j;
        j++;
      }
      π[k] = ι[k] = k;
      k--;
    }
    break;

  moved:;
  } while (true);
  return {L, I.c + 1};
}

template <typename Rule, typename Order> Solution assignOnceOpt(const Instance &I, unsigned short lb, Rule rule, Order order) {

  AssignContext ctx(I);

  ctx.setup();

  unsigned ti = rule(ctx, true);
  Task t = ctx.A[ti];
  ctx.first(t);
  ctx.remove(ti);

  while (ctx.A.size() > 0) {
    const bool first = false;


    vector<short> pi = order(ctx, first);
    bool extended = false;

    for (short p : pi) {
      Task j = ctx.A[p];


      if (ctx.l + ctx.I.t[j] + ctx.set[j] <= I.c) {
        ctx.next(j);
        ctx.remove(p);
        extended = true;
        break;
      }

      const unsigned max_tasks_reopt = 7;
      if (ctx.pi.size() - ctx.iti <= max_tasks_reopt) {
        unsigned light_job = ctx.pi.back();

        vector<Task> L(ctx.pi.begin() + ctx.iti, ctx.pi.end());
        L.push_back(j);
        Time load;
        tie(L, load) = reOrder_prec(ctx.I, L);

        if (load <= I.c) {
          ctx.pi.push_back(j);
          copy(L.begin(), L.end(), ctx.pi.begin() + ctx.iti);
          ctx.l = load;
          ctx.it = L[0];
          ctx.ft = L.back();
          if (L.back() == j)
            ctx.shadow(light_job);
          else if (L.back() == light_job)
            ctx.shadow(j);
          else {
            ctx.shadow(light_job);
            ctx.shadow(j);
            ctx.unshadow(L.back());
          }
          ctx.remove(p);
          extended = true;
          break;
        }
      } else {
      }
    }

    if (!extended) {
      ctx.open();
      if (ctx.m >= incumbent.m) {
        ctx.pi.clear();
        return Solution{ctx.pi, uint_inf};
      }
      pi = order(ctx, true);
      ctx.first(ctx.A[pi[0]]);
      ctx.remove(pi[0]);
    }
  }

  return Solution{ctx.pi, ctx.m};
}

tuple<rule,ruleOrder,bool> get_rule(string rule) {
  auto e = rules.find(rule);
  if (e != rules.end())
    return e->second;

  if (sbf::rule_grasp::rules.size() == 0) {
    istringstream cr(rule);
    string token;
    getline(cr, token, '(');
    if (token != "rule-grasp") throw "Unknown rule.";
    while (getline(cr, token, ',')) {
      if (token.back() == ')') token.pop_back();
      auto r = get_rule(token);
      sbf::rule_grasp::rules.push_back({ token, get<0>(r), get<1>(r), get<2>(r) });
    }
  }
  return { sbf::rule_grasp::best, sbf::rule_grasp::order, sbf::rule_grasp::stochastic };
}

unsigned short assignRule(Instance &I, unsigned short lb, string rule) {
  auto e = get_rule(rule);

  unsigned ilim = opt.ilim;
  if (get<2>(e) == false)
    ilim = S.iter + 1;

  while (logging::elapsed() < opt.tlim && S.iter < ilim && incumbent.m > lb) {
    S.iter++;
    Solution s = assignOnce(I, lb, get<0>(e));
    incumbent.set(I, s);
  }
  return incumbent.m;
}

unsigned short assignRuleOpt(Instance &I, unsigned short lb, string rule) {
  auto e = get_rule(rule);

  unsigned ilim = opt.ilim;
  if (get<2>(e) == false)
    ilim = S.iter + 1;

  while (logging::elapsed() < opt.tlim && S.iter++ < ilim && incumbent.m > lb) {
    Solution s = assignOnceOpt(I, lb, get<0>(e), get<1>(e));
    incumbent.set(I, s);
  }
  return incumbent.m;
}

unsigned short assignRuleAlt(Instance &I, Instance& Ir, unsigned short lb, string rule) {
  auto e = get_rule(rule);

  unsigned ilim = opt.ilim;
  if (get<2>(e) == false)
    ilim = S.iter + 1;

  while (logging::elapsed() < opt.tlim && S.iter < ilim && incumbent.m > lb) {
    S.iter++;
    Solution s = assignOnce(I, lb, get<0>(e));
    incumbent.set(I, s);
    s = assignOnce(Ir, lb, get<0>(e));
    if (incumbent.set(I, s)) {
      reverse(s.pi.begin(), s.pi.end());
    }
  }
  return incumbent.m;
}

unsigned short assignRuleOptAlt(Instance &I, Instance& Ir, unsigned short lb, string rule) {
  auto e = get_rule(rule);

  unsigned ilim = opt.ilim;
  if (get<2>(e) == false)
    ilim = S.iter + 1;

  while (logging::elapsed() < opt.tlim && S.iter++ < ilim && incumbent.m > lb) {
    Solution s = assignOnceOpt(I, lb, get<0>(e), get<1>(e));
    incumbent.set(I, s);
    s = assignOnceOpt(I, lb, get<0>(e), get<1>(e));
    if (incumbent.set(I, s)) {
      reverse(s.pi.begin(), s.pi.end());
    }
  }
  return incumbent.m;
}
