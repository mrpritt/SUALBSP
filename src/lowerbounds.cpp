#include "lowerbounds.hpp"

#include <fstream>
using namespace std;

#include "lowerbounds.hpp"
#include "options.hpp"
#include "solution.hpp"

namespace sualbsp {

vector<tail> tails;
vector<tail> heads;

void new_incumbent(Instance &I) {
  I.preprocess_dynamic(incumbent.m);

  if (opt.oinstance != "") {
    ofstream out(opt.oinstance);
    I.writeSBF(out);
    out.close();
  }

  unsigned short lm1v = lm1(I.t, I.c);
  unsigned short lms1v = lms1(I);
  auto lb = max(incumbent.lb, max(lm1v, lms1v));
  if (lb > incumbent.lb)
  incumbent.lb = lb;
}

unsigned lms1_fb_subset(const Instance &I, const vector<bool> &active, bool istails) {
  vector<unsigned> fs;
  vector<unsigned> bs;
  for (auto i = 0u; i < I.n; i++) {
    if (active[i]) {
      auto smallest = I.c;
      for (auto j = 0u; j < I.n; j++) {
        if (istails) {
          if ((i != j) && (active[j]) && (smallest > I.sf[i][j])) {
            smallest = I.sf[i][j];
          }
        } else {
          if ((i != j) && (active[j]) && (smallest > I.sf[j][i])) {
            smallest = I.sf[j][i];
          }
        }
      }
      fs.push_back(smallest);
    }
  }
  sort(fs.begin(), fs.end(), [](const unsigned &a, const unsigned &b) { return a < b; });
  for (auto i = 0u; i < I.n; i++) {
    if (active[i]) {
      auto smallest = I.c;
      for (auto j = 0u; j < I.n; j++) {
        if (istails) {
          if ((active[j]) && (smallest > I.sb[i][j])) {
            smallest = I.sb[i][j];
          } else {
            if ((active[j]) && (smallest > I.sb[j][i])) {
              smallest = I.sb[j][i];
            }
          }
        }
      }
      bs.push_back(smallest);
    }
  }
  sort(bs.begin(), bs.end(), [](const unsigned &a, const unsigned &b) { return a < b; });
  unsigned t_sum = 0;
  unsigned n = 0;
  for (auto i = 0u; i < I.n; i++) {
    if (active[i]) {
      t_sum += I.t[i];
      n++;
    }
  }
  unsigned m = round(t_sum, I.c);
  unsigned ta_sum = accumulate(fs.begin(), fs.begin() + n - m, 0);
  unsigned mu_sum = accumulate(bs.begin(), bs.begin() + m, 0);
  unsigned TC = I.c * m;
  while (t_sum + ta_sum + mu_sum > TC) {
    if (fs.size() < (n - m))
      return (m);
    ta_sum = accumulate(fs.begin(), fs.begin() + n - m, 0);
    TC += I.c;
    m++;
    if (bs.size() < m)
      return (m);
    mu_sum = accumulate(bs.begin(), bs.begin() + m, 0);
  }
  return m;
}

unsigned lms1_f_subset(const Instance &I, const vector<bool> &active, bool istails) {
  vector<unsigned> fs;
  for (auto i = 0u; i < I.n; i++) {
    if (active[i]) {
      auto smallest = I.c;
      for (auto j = 0u; j < I.n; j++) {
        if (istails) {
          if ((i != j) && (active[j]) && (smallest > I.sf[i][j])) {
            smallest = I.sf[i][j];
          }
        } else {
          if ((i != j) && (active[j]) && (smallest > I.sf[j][i])) {
            smallest = I.sf[j][i];
          }
        }
      }
      fs.push_back(smallest);
    }
  }
  sort(fs.begin(), fs.end(), [](const unsigned &a, const unsigned &b) { return a < b; });
  unsigned t_sum = 0;
  unsigned n = 0;
  for (auto i = 0u; i < I.n; i++) {
    if (active[i]) {
      t_sum += I.t[i];
      n++;
    }
  }
  unsigned m = round(t_sum, I.c);
  if (m == n)
    return m;
  Time ta_sum = accumulate(fs.begin(), fs.begin() + n + 1 - m, 0);
  Time TC = I.c * m;
  while (t_sum + ta_sum > TC) {
    ta_sum -= fs[n - m];
    m++;
    TC += I.c;
  }
  return m;
}

unsigned lms1_fb(const Instance &I) {
  Time tsum = lm1Unrounded(I.t, I.c);

  unsigned m = round(tsum, I.c);
  Time ta_sum = accumulate(I.sfi.begin(), I.sfi.begin() + I.n - m, 0, plus<Time>());
  Time mu_sum = accumulate(I.sbi.begin(), I.sbi.begin() + m, 0, plus<Time>());
  Time TC = I.c * m;

  while (tsum + ta_sum + mu_sum > TC) {
    ta_sum -= I.sfi[I.n - m - 1];
    TC += I.c;
    mu_sum += I.sbi[m];
    m++;
  }
  return m;
}

unsigned lms1_f(const Instance &I) {
  Time tsum = lm1Unrounded(I.t, I.c);

  unsigned m = round(tsum, I.c);
  if (m == I.n)
    return m;

  Time ta_sum = accumulate(I.sfi.begin(), I.sfi.begin() + I.n + 1 - m, 0, plus<Time>());
  Time TC = I.c * m;

  while (tsum + ta_sum > TC) {
    ta_sum -= I.sfi[I.n - m];
    m++;
    TC += I.c;
  }
  return m;
}

unsigned lm2(const Instance &I) {
  auto count = 0;
  auto countmid = 0;
  for (auto i = 0u; i < I.n; i++) {
    if (I.t[i] * 2 > I.c)
      count++;
    else if (I.t[i] * 2 == I.c)
      countmid++;
  }
  return count + (countmid + 1) / 2;
}

unsigned lm3(const Instance &I) {
  auto cg1 = 0;
  auto cg2 = 0;
  auto cg3 = 0;
  auto cg4 = 0;
  for (auto i = 0u; i < I.n; i++) {
    if (I.t[i] * 3 > I.c * 2) {
      cg1++;
    } else {
      if (I.t[i] * 3 == I.c * 2) {
        cg2++;
      } else {
        if (I.t[i] * 3 > I.c) {
          cg3++;
        } else {
          if (I.t[i] * 3 == I.c) {
            cg4++;
          }
        }
      }
    }
  }
  if ((6 * cg1 + 4 * cg2 + 3 * cg3 + 2 * cg4) % 6 == 0) {
    return ((6 * cg1 + 4 * cg2 + 3 * cg3 + 2 * cg4) / 6);
  } else {
    return ((6 * cg1 + 4 * cg2 + 3 * cg3 + 2 * cg4) / 6 + 1);
  }
}

unsigned lm4(const Instance &I) {
  for (auto ii = 0u; ii < I.n; ii++) {
    tail i;
    i.id = ii;
    i.t = 0;
    for (auto j : I.Fs[ii]) {
      i.t += I.t[j];
    }
    tails.push_back(i);
  }
  sort(tails.begin(), tails.end(), [](const tail &a, const tail &b) { return a.t > b.t; });
  for (int task = I.n - 1; task >= 0; task--) {
    vector<bool> successors(I.n, false);
    for (auto i : I.Fs[task]) {
      successors[i] = true;
    }
    auto acc = 0u;
    auto maxacc = 0u;
    for (auto i : tails) {
      if (successors[i.id]) {
        acc += I.t[i.id];
        if ((acc + i.t) > maxacc) {
          maxacc = acc + i.t;
        }
      }
    }
    if(maxacc>0) {
      if((I.t[task]==I.c)||(maxacc%I.c==0)) {
        if(maxacc%I.c!=0) maxacc+=I.c-(maxacc%I.c);
      } else {
        auto minfs=I.c;
        auto minbs=I.c;
        if (I.isTriangular) {
          for(auto i : I.F[task])
            if(minfs>I.sf[task][i]) minfs=I.sf[task][i];
          for(auto i : I.Fs[task])
            if(minbs>I.sb[i][task]) minbs=I.sb[i][task];
        }
        if (!I.isTriangular) {
          for(auto i : I.Fs[task])
            if(minfs>I.sf[task][i]) minfs=I.sf[task][i];
          for(auto i : I.Is[task])
            if(minfs>I.sf[task][i]) minfs=I.sf[task][i];
          for(auto i : I.Fs[task])
            if(minbs>I.sb[i][task]) minbs=I.sb[i][task];
          for(auto i : I.Is[task])
            if(minbs>I.sb[i][task]) minbs=I.sb[i][task];
        }
        if((minfs+minbs>=I.c)&&(maxacc%I.c!=0)) {
          maxacc = (maxacc / I.c +1) * I.c;
        } else {
          if((maxacc + I.t[task] + minfs + minbs) > ((maxacc / I.c + 1) * I.c)) {
            maxacc = (maxacc / I.c + 1) * I.c;
          } else {
            maxacc = maxacc + minfs;
          }
        }
      }
    }
    for (auto &i : tails) {
      if (i.id == task) {
        i.t = maxacc;
        break;
      }
    }
    sort(tails.begin(), tails.end(), [](const tail &a, const tail &b) { return a.t > b.t; });
  }

  auto lm4tails = 0;
  auto acctails = 0u;
  auto maxacctails = 0u;
  for (auto i : tails) {
    acctails += I.t[i.id];
    if ((acctails + i.t) > maxacctails) {
      maxacctails = acctails + i.t;
    }
  }
  if (maxacctails % I.c == 0) {
    lm4tails = maxacctails / I.c;
  } else {
    lm4tails = maxacctails / I.c + 1;
  }
  for (auto ii = 0u; ii < I.n; ii++) {
    tail i;
    i.id = ii;
    i.t = 0;
    for (auto j : I.Ps[ii]) {
      i.t += I.t[j];
    }
    heads.push_back(i);
  }
  sort(heads.begin(), heads.end(), [](const tail &a, const tail &b) { return a.t > b.t; });
  for (auto task = 0u; task < I.n; task++) {
    vector<bool> predecessors(I.n, false);
    for (auto i : I.Ps[task]) {
      predecessors[i] = true;
    }
    auto acc = 0u;
    auto maxacc = 0u;
    for (auto i : heads) {
      if (predecessors[i.id]) {
        acc += I.t[i.id];
        if ((acc + i.t) > maxacc) {
          maxacc = acc + i.t;
        }
      }
    }
    if(maxacc>0) {
      if((I.t[task]==I.c)||(maxacc%I.c==0)) {
        if(maxacc%I.c!=0) maxacc+=I.c-(maxacc%I.c);
      } else {
        auto minfs=I.c;
        auto minbs=I.c;
        if(I.isTriangular) {
          for(auto i : I.P[task])
            if(minfs>I.sf[i][task]) minfs=I.sf[i][task];
          for(auto i : I.Ps[task])
            if(minbs>I.sb[task][i]) minbs=I.sb[task][i];
        }
        if(!I.isTriangular) {
          for(auto i : I.P[task])
            if(minfs>I.sf[i][task]) minfs=I.sf[i][task];
          for(auto i : I.Is[task])
            if(minfs>I.sf[i][task]) minfs=I.sf[i][task];
          for(auto i : I.Ps[task])
            if(minbs>I.sb[task][i]) minbs=I.sb[task][i];
          for(auto i : I.Is[task])
            if(minbs>I.sb[task][i]) minbs=I.sb[task][i];
        }
        if((minfs+minbs>=I.c)&&(maxacc%I.c!=0)) {
          maxacc = (maxacc / I.c +1) * I.c;
        } else {
          if((maxacc + I.t[task] + minfs + minbs) > ((maxacc / I.c + 1) * I.c)) {
            maxacc = (maxacc / I.c + 1) * I.c;
          } else {
            maxacc = maxacc + minfs;
          }
        }
      }
    }
    for (auto &i : heads) {
      if (i.id == task) {
        i.t = maxacc;
        break;
      }
    }
    sort(heads.begin(), heads.end(), [](const tail &a, const tail &b) { return a.t > b.t; });
  }
  auto lm4heads = 0;
  auto accheads = 0u;
  auto maxaccheads = 0u;
  for (auto i : heads) {
    accheads += I.t[i.id];
    if ((accheads + i.t) > maxaccheads) {
      maxaccheads = accheads + i.t;
    }
  }
  if (maxaccheads % I.c == 0) {
    lm4heads = maxaccheads / I.c;
  } else {
    lm4heads = maxaccheads / I.c + 1;
  }

  if (lm4tails > lm4heads)
    return lm4tails;
  return lm4heads;
}

Time lc_tasks(const Instance &I, const vector<Task> &T0) {
  Time load = 0;

  BitVector T(I.n, false);
  for (auto i : T0)
    T[i] = true;
  for (auto i : T0) {
    for (auto j : T0) {
      if (j > i)
        continue;
      for (auto k : I.O[i][j])
        T[k] = true;
    }
  }

  for (Task i = 0; i != I.n; ++i)
    if (T[i])
      load += I.t[i];

  if (!opt.heavypre)
    return load;

  BitVector valid(I.n);
  if (I.isTriangular)
    valid = T;
  else
    valid.set();

  Time Sω = 0, ωmax = 0;
  int last = -1;
  vector<Time> ωτ(I.n, 0);
  for (Task i = 0; i != I.n; ++i)
    if (T[i]) {
      ωτ[i] = I.get_ωτ(i, valid);
      if (ωτ[i] == inf_time) {
        if (last < 0)
          last = i;
        else
          return I.c + 1;
      } else {
        Sω += ωτ[i];
        ωmax = max(ωmax, ωτ[i]);
      }
    }

  Time Sα = 0, αmax = 0;
  vector<Time> ατ(I.n, 0);
  int first = -1;
  for (Task i = 0; i != I.n; ++i)
    if (T[i]) {
      ατ[i] = I.get_ατ(i, valid);
      if (ατ[i] == inf_time) {
        if (first < 0)
          first = i;
        else
          return I.c + 1;
      } else {
        Sα += ατ[i];
        αmax = max(αmax, ατ[i]);
      }
    }

  if (I.isTriangular) {
    Time SUα = inf_time, SUω = inf_time;
    int fbegin = 0, fend = I.n;
    int lbegin = 0, lend = I.n;
    if (first >= 0) {
      fbegin = first;
      fend = first + 1;
    }
    if (last >= 0) {
      lbegin = last;
      lend = last + 1;
    }
    for (Task i = fbegin; i != fend; ++i) {
      if (!T[i])
        continue;
      for (Task j = lbegin; j != lend; ++j) {
        if (i == j || !T[j] || I.D[j][i])
          continue;
        SUα = min(SUα, Sα - (int(i) == first ? 0 : ατ[i]) + I.sb[j][i]);
        SUω = min(SUω, Sω - (int(j) == last ? 0 : ωτ[j]) + I.sb[j][i]);
      }
    }
    load += max(SUα, SUω);
  } else {
    if (last < 0)
      Sω -= ωmax;
    if (first < 0)
      Sα -= αmax;
    load += max(Sω, Sα);

    load += I.sbi[0];
  }
  return load;
}
}
