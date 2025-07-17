#include "bbr.hpp"
#include "logging.hpp"
#include "lowerbounds.hpp"
#include "solution.hpp"

namespace bbr {

inline unsigned round(unsigned t, unsigned c) { return (t + c - 1) / c; }

struct compare {
  bool operator()(const std::vector<short> &a, const std::vector<short> &b) const {
    for (auto i = 0u; i < a.size(); i++) {
      if ((a[i] < 0) && (b[i] >= 0))
        return false;
      if ((a[i] >= 0) && (b[i] < 0))
        return true;
    }
    return false;
  }
};

unsigned iter;

map<vector<short int>, unsigned, compare> subproblems;
typedef priority_queue<std::pair<unsigned, std::vector<short int>>, std::vector<std::pair<unsigned, std::vector<short int>>>> PriorityQueue;
vector<PriorityQueue> PriorityQueues;

vector<setupStruct> fs;
vector<setupStruct> bs;
vector<unsigned> lm2v;
vector<unsigned> lm3v;

vector<short int> deg;
vector<short int> sequence;
vector<short int> pi;
unsigned idle;
unsigned maxNodesCbfs;
unsigned nodesCbfs;
unsigned sumTaskTimes;
unsigned numCompleteTasks;
unsigned numStations;
unsigned tsum;

void initializebounds(const Instance &I) {
  setupStruct item;

  for (auto i = 0u; i < I.n; i++) {
    for (auto j = 0u; j < I.n; j++) {
      if ((i != j) && (I.sf[i][j] < I.c)) {
        item.t1 = i;
        item.t2 = j;
        item.t = I.sf[i][j];
        fs.push_back(item);
      }
    }
  }
  sort(fs.begin(), fs.end(), [](const setupStruct &a, const setupStruct &b) { return a.t < b.t; });
  for (auto i = 0u; i < I.n; i++) {
    for (auto j = 0u; j < I.n; j++) {
      if (I.sb[i][j] < I.c) {
        item.t1 = i;
        item.t2 = j;
        item.t = I.sb[i][j];
        bs.push_back(item);
      }
    }
  }
  sort(bs.begin(), bs.end(), [](const setupStruct &a, const setupStruct &b) { return a.t < b.t; });
  if (opt.usecounting) {
    lm2v.resize(I.n);
    fill(lm2v.begin(), lm2v.end(), 0);
    for (auto i = 0u; i < I.n; i++) {
      if (I.t[i] * 2 > I.c) {
        lm2v[i] = 2;
      } else {
        if (I.t[i] * 2 == I.c) {
          lm2v[i] = 1;
        }
      }
    }
    lm3v.resize(I.n);
    fill(lm3v.begin(), lm3v.end(), 0);
    for (auto i = 0u; i < I.n; i++) {
      if (I.t[i] * 3 > I.c * 2) {
        lm3v[i] = 6;
      } else {
        if (I.t[i] * 3 == I.c * 2) {
          lm3v[i] = 4;
        } else {
          if (I.t[i] * 3 > I.c) {
            lm3v[i] = 2;
          } else {
            if (I.t[i] * 3 == I.c) {
              lm3v[i] = 1;
            }
          }
        }
      }
    }
  }
}

void evaluateSetups(const vector<setupStruct> &s, vector<unsigned> &v, unsigned n) {
  unsigned numTasks = n;
  vector<bool> useful(n, true);
  vector<bool> initiallyuseful(n, true);
  v.clear();
  for (auto i = 0u; i < n; i++) {
    if (deg[i] < 0) {
      useful[i] = false;
      initiallyuseful[i] = false;
      numTasks--;
    }
  }
  for (auto i = 0u; i < s.size(); i++) {
    if ((useful[s[i].t1] == true) && (initiallyuseful[s[i].t2] == true)) {
      v.push_back(s[i].t);
      useful[s[i].t1] = false;
      numTasks--;
      if (numTasks == 0)
        break;
    }
  }
}

unsigned dobound(const Instance &I) {
  if ((numCompleteTasks + 1) == I.n)
    return 1;
  if ((!opt.uselm1p) && (!opt.usecounting) && (!opt.uselm4)) {
    return (round(tsum - sumTaskTimes, I.c));
  }
  if (round(tsum - sumTaskTimes, I.c) + numStations >= incumbent.m)
    return (round(tsum - sumTaskTimes, I.c));

  if (opt.usecounting) {
    auto m2 = 0;
    auto m3 = 0;
    for (auto i = 0u; i < I.n; i++) {
      if (deg[i] >= 0) {
        m2 += lm2v[i];
        m3 += lm3v[i];
      }
    }
    if (m2 % 2 == 1)
      m2 = (m2 + 1) / 2;
    else
      m2 = m2 / 2;
    if ((numStations + m2) >= incumbent.m) {
      return m2;
    }
    if (m3 % 6 > 0)
      m3 = (m3 + 6) / 6;
    else
      m3 = m3 / 6;
    if ((numStations + m3) >= incumbent.m) {
      return m3;
    }
  }

  unsigned m = 0;
  if (opt.uselm1p) {
    if (I.setupDirected) {
      m = round(tsum - sumTaskTimes, I.c);
      if ((numStations + m) >= incumbent.m)
        return m;
      vector<unsigned> forwards;
      evaluateSetups(fs, forwards, I.n);
      if (forwards.size() < (I.n - numCompleteTasks - m))
        return (incumbent.m);
      unsigned ta_sum = accumulate(forwards.begin(), forwards.begin() + I.n - numCompleteTasks - m, 0);

      vector<unsigned> backwards;
      evaluateSetups(bs, backwards, I.n);
      if (backwards.size() < m)
        return (incumbent.m);
      unsigned mu_sum = accumulate(backwards.begin(), backwards.begin() + m, 0);

      unsigned TC = I.c * m;
      while (tsum - sumTaskTimes + ta_sum + mu_sum > TC) {
        if (forwards.size() < (I.n - numCompleteTasks - m))
          return (incumbent.m);
        ta_sum = accumulate(forwards.begin(), forwards.begin() + I.n - numCompleteTasks - m, 0);

        TC += I.c;
        m++;
        if (backwards.size() < m)
          return (incumbent.m);
        mu_sum = accumulate(backwards.begin(), backwards.begin() + m, 0);
        if ((numStations + m) >= incumbent.m)
          return m;
      }
    } else {
      m = round(tsum - sumTaskTimes, I.c);
      if ((numStations + m) >= incumbent.m)
        return m;

      if (m == I.n)
        return m;
      vector<unsigned> forwards;
      evaluateSetups(fs, forwards, I.n);
      if ((I.n - numCompleteTasks) == 1)
        return (1);
      while (forwards.size() < (I.n - numCompleteTasks - m)) {
        forwards.push_back(0);
      }
      unsigned ta_sum = accumulate(forwards.begin(), forwards.begin() + I.n - numCompleteTasks - m, 0);

      unsigned TC = I.c * m;
      if (tsum - sumTaskTimes + ta_sum <= TC)
        return 1;

      while (tsum - sumTaskTimes + ta_sum > TC) {
        m++;
        if ((numStations + m) >= incumbent.m)
          return m;
        ta_sum = accumulate(forwards.begin(), forwards.begin() + I.n - numCompleteTasks + 1 - m, 0);
        TC += I.c;
      }
    }
  }

  if (opt.uselm4) {
    auto acc = 0u;
    auto maxacc = 0u;
    auto m4 = 0u;
    for (auto t : tails) {
      if (deg[t.id] >= 0) {
        acc += I.t[t.id];
        if ((acc + t.t) > maxacc)
          maxacc = (acc + t.t);
      }
    }
    if (maxacc % I.c > 0)
      m4 = maxacc / I.c + 1;
    else
      m4 = maxacc / I.c;
    if (m4 > m)
      m = m4;
  }
  return (m);
}

bool doEnhancedLookup(const Instance &I) {
  for (auto i = 0u; i < I.n; i++) {
    if (deg[i] == 0) {
      deg[i] = -1;
      auto it = subproblems.find(deg);
      if (it != subproblems.end()) {
        if (it->second <= numStations) {
          deg[i] = 0;
          return (true);
        }
      }
      deg[i] = 0;
    }
  }
  return (false);
}

void gen_loads(const Instance &I) {
  unsigned add;
  bool withExtension = false;

  vector<unsigned> free;
  for (auto i = 0u; i < I.n; i++)
    if (deg[i] == 0)
      free.push_back(i);

  auto mprule = [&](unsigned t) { return 5.0 * I.t[t] + 0.3 * I.tn[t] - 45.3 * (sequence.size() == 0 ? 0 : I.sf[sequence.back()][t]); };
  sort(free.begin(), free.end(), [&](unsigned a, unsigned b) { return mprule(a) > mprule(b); });

  for (auto i : free) {
    if (deg[i] == 0) {
      if (sequence.size() == 0) {
        add = I.t[i] + I.sb[i][i];
      } else {
        const auto last = sequence[sequence.size() - 1];
        add = I.t[i] - I.sb[last][sequence[0]] + I.sf[last][i] + I.sb[i][sequence[0]];
      }
      if (add <= idle) {
        withExtension = true;
        sequence.push_back(i);
        numCompleteTasks++;
        pi[i] = 0 - numCompleteTasks;
        idle -= add;
        sumTaskTimes += I.t[i];
        deg[i]--;
        for (auto j = 0u; j < I.F[i].size(); j++) {
          deg[I.F[i][j]]--;
        }
        if (nodesCbfs < maxNodesCbfs)
          gen_loads(I);
        for (auto j = 0u; j < I.F[i].size(); j++) {
          deg[I.F[i][j]]++;
        }
        deg[i]++;
        sumTaskTimes -= I.t[i];
        idle += add;
        pi[i] = 0;
        numCompleteTasks--;
        sequence.pop_back();
      }
    }
  }

  if (withExtension == false) {
    nodesCbfs++;
    if (numCompleteTasks < I.n) {
      if ((numStations + dobound(I)) < incumbent.m) {
        std::map<vector<short int>, unsigned>::iterator it = subproblems.find(pi);
        if ((it == subproblems.end()) || (it->second > numStations)) {
          if ((opt.enhancedlookup == false) || (doEnhancedLookup(I) == false)) {
            if (it == subproblems.end()) {
              PriorityQueues[numStations].push(make_pair(sumTaskTimes, deg));
              subproblems[pi] = numStations;
            } else {
              if (it->second > numStations) {
                PriorityQueues[numStations].push(make_pair(sumTaskTimes, deg));
                subproblems.erase(it);
                subproblems[pi] = numStations;
              } else {
              }
            }
          }
        }
      }
    } else {
      if (incumbent.m > numStations) {
        incumbent.iter = iter;
        incumbent.time = logging::elapsed();
        incumbent.m = numStations;
        for (auto i = 0u; i < I.n; i++) {
          incumbent.pi[0 - pi[i] - 1] = i;
        }
        auto currentm = 1;
        auto task = incumbent.pi[0];
        auto idle = I.c - I.t[task] - I.sb[task][task];
        incumbent.a[0] = currentm;
        auto first = task;
        for (auto position = 1u; position < I.n; position++) {
          task = incumbent.pi[position];
          const auto last = incumbent.pi[position - 1];
          auto remove = I.t[task] - I.sb[last][first] + I.sf[last][task] + I.sb[task][first];
          if (remove <= idle) {
            idle -= remove;
            incumbent.a[position] = currentm;
          } else {
            idle = I.c - I.t[task] - I.sb[task][task];
            currentm++;
            incumbent.a[position] = currentm;
            first = task;
          }
        }
      }
    }
  }
  return;
}

bool bbr(const Instance &I, Options &opt, unsigned short lb) {
  bool isOptimal = true;
  bool change;
  vector<short int> degrees;
  maxNodesCbfs = opt.maxnodescbfs;
  if (subproblems.empty() == false)
    subproblems.clear();
  if (PriorityQueues.empty() == false)
    PriorityQueues.clear();
  if (fs.empty() == false)
    fs.clear();
  if (bs.empty() == false)
    bs.clear();
  if (pi.empty() == false)
    pi.clear();
  if (lm2v.empty() == false)
    lm2v.clear();
  if (lm3v.empty() == false)
    lm3v.clear();
  pi.reserve(I.n);
  tsum = 0;
  for (auto i = 0u; i < I.n; i++)
    tsum += I.t[i];
  opt.tlimexact += logging::elapsed();
  for (auto i = 0u; i < incumbent.m; i++) {
    PriorityQueues.push_back(PriorityQueue());
  }
  initializebounds(I);

  for (auto i = 0u; i < I.n; i++)
    degrees.push_back(I.P[i].size());
  PriorityQueues[0].push(make_pair(0, degrees));
  for (auto i = 0u; i < I.n; i++)
    pi.push_back(0);
  subproblems[pi] = 0;

  iter = 0;
  do {
    iter++;
    change = false;
    for (auto depth = 0u; depth < incumbent.m - 1; depth++) {
      if (PriorityQueues[depth].empty())
        continue;
      numStations = depth + 1;
      change = true;
      deg = PriorityQueues[depth].top().second;
      sequence.clear();
      PriorityQueues[depth].pop();
      if ((opt.enhancedlookup == true) && (doEnhancedLookup(I) == true))
        continue;
      std::map<vector<short int>, unsigned>::iterator it = subproblems.find(deg);
      if (it->second > depth) {
        continue;
      } else {
        pi = it->first;
        nodesCbfs = 0;
        idle = I.c;
        sumTaskTimes = 0;
        numCompleteTasks = 0;
        for (auto i = 0u; i < I.n; i++) {
          if (deg[i] < 0) {
            sumTaskTimes += I.t[i];
            numCompleteTasks++;
          }
        }
        auto bound = (tsum - sumTaskTimes) / I.c;
        if ((tsum - sumTaskTimes) % I.c != 0)
          bound++;
        if ((depth + bound) < incumbent.m) {
          gen_loads(I);
          if (nodesCbfs >= maxNodesCbfs) {
            isOptimal = false;
          }
        }
      }
    }
    if (incumbent.lb == incumbent.m) {
      isOptimal = true;
      break;
    }

    if (logging::elapsed() > opt.tlimexact) {
      isOptimal = false;
      break;
    }
    if (opt.maxmemory < subproblems.size()) {
      isOptimal = false;
      break;
    }
  } while (change == true);
  return isOptimal;
}

} // namespace bbr
