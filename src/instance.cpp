#include "instance.hpp"

#include <cassert>
#include <filesystem>
#include <iostream>
#include <regex>
#include <utility>
namespace fs = std::filesystem;
using namespace std;

#include "logging.hpp"
#include "lowerbounds.hpp"
#include "options.hpp"
#include "util.hpp"

namespace sualbsp {

regex blank_line("^\\s*\\r?$");
regex section_start("<(.*)>\\r?");

constexpr bool disable_preprocessing = false;

Instance::Instance(Instance &&other) { swap(other); }

Instance::Instance(const Instance &other) {
  n = other.n;
  c = other.c;
  t = other.t;
  t0 = other.t0;
  P = other.P;
  F = other.F;
  d = other.d;
  Ps = other.Ps;
  Fs = other.Fs;
  Is = other.Is;
  O = other.O;
  D = other.D;
  sf.resize(boost::extents[other.sf.shape()[0]][other.sf.shape()[1]]);
  sf = other.sf;
  sf0.resize(boost::extents[other.sf0.shape()[0]][other.sf0.shape()[1]]);
  sf0 = other.sf0;
  sb.resize(boost::extents[other.sb.shape()[0]][other.sb.shape()[1]]);
  sb = other.sb;
  sb0.resize(boost::extents[other.sb0.shape()[0]][other.sb0.shape()[1]]);
  sb0 = other.sb0;
  sfi = other.sfi;
  sbi = other.sbi;
  setupDirected = other.setupDirected;
  fwdSymmetric = other.fwdSymmetric;
  bwdSymmetric = other.bwdSymmetric;
  isTriangular = other.isTriangular;
  isStrictTriangular = other.isStrictTriangular;
  insertionlb = other.insertionlb;
  ta = other.ta;
  tn = other.tn;
  E = other.E;
  T = other.T;
}

Instance &Instance::operator=(Instance other) {
  swap(other);
  return *this;
}

void Instance::swap(Instance &other) {
  using std::swap;
  swap(n, other.n);
  swap(c, other.c);
  swap(t, other.t);
  swap(t0, other.t0);
  P.swap(other.P);
  F.swap(other.F);
  d.swap(other.d);
  Ps.swap(other.Ps);
  Fs.swap(other.Fs);
  Is.swap(other.Is);
  O.swap(other.O);
  D.swap(other.D);
  sf.resize(boost::extents[other.sf.shape()[0]][other.sf.shape()[1]]);
  swap(sf, other.sf);
  sb.resize(boost::extents[other.sb.shape()[0]][other.sb.shape()[1]]);
  swap(sb, other.sb);
  sf0.resize(boost::extents[other.sf0.shape()[0]][other.sf0.shape()[1]]);
  swap(sf0, other.sf0);
  sb.resize(boost::extents[other.sb.shape()[0]][other.sb.shape()[1]]);
  swap(sb, other.sb);
  sb0.resize(boost::extents[other.sb0.shape()[0]][other.sb0.shape()[1]]);
  swap(sb0, other.sb0);
  sfi.swap(other.sfi);
  sbi.swap(other.sbi);
  swap(setupDirected, other.setupDirected);
  swap(fwdSymmetric, other.fwdSymmetric);
  swap(bwdSymmetric, other.bwdSymmetric);
  swap(isTriangular, other.isTriangular);
  swap(isStrictTriangular, other.isStrictTriangular);
  insertionlb.swap(other.insertionlb);
  ta.swap(other.ta);
  tn.swap(other.tn);
  E.swap(other.E);
  T.swap(other.T);
}

void Instance::allocate() {
  t.resize(n);
  t0.resize(n);
  d.resize(n, BitVector(n));
  sf0.resize(extents[n][n]);
  sb0.resize(extents[n][n]);
  sf.resize(extents[n][n]);
  sb.resize(extents[n][n]);
}

void Instance::readMP(istream &in, bool do_preprocess) {
  char comma;

  in >> n;
  allocate();

  unsigned p;
  in >> p;
  in >> c;
  for (Task i = 0; i != n; ++i) {
    unsigned id;
    in >> id >> comma >> t0[i];
  }
  while (p > 0) {
    unsigned pred, succ;
    in >> pred >> comma >> succ;
    d[pred][succ] = true;
    p--;
  }
  for (Task i = 0; i != n; ++i)
    for (Task j = 0; j != n; ++j) {
      in >> sf0[i][j];
      sb0[i][j] = sf0[i][j];
      if (j < n - 1) {
        in >> comma;
      }
    }
  if (opt.type2)
    c = uc();
  setupDirected = false;

  t = t0;
  sf = sf0;
  sb = sb0;

  preprocess(do_preprocess);
}

void skipTo(istream &in, string tag) {
  regex thisSection("<" + tag + ">\\r?");
  while (!in.eof()) {
    string line;
    getline(in, line);
    if (regex_match(line, thisSection))
      return;
  }
  throw tag;
}

string skipToNext(istream &in) {
  while (!in.eof()) {
    string line;
    getline(in, line);
    smatch tag;
    if (regex_match(line, tag, section_start)) {
      return tag[1].str();
    }
  }
  return "end";
}

bool isEmpty(string line) { return regex_match(line, blank_line); }

Time Instance::uc() {
  Time r = accumulate(t.begin(), t.end(), 0);
  Time maxsu = 0;
  for (Task i = 0; i != n; ++i)
    for (Task j = 0; j != n; ++j)
      maxsu = max(maxsu, max(sf[i][j], sb[i][j]));
  return r + n * maxsu;
}

void Instance::handleSetupSection(istream &in, multi_array<Time, 2> &s) {
  char comma, colon;
  string line;
  unsigned pred, succ;

  std::fill_n(s.data(), s.num_elements(), 0);

  while (true) {
    getline(in, line);
    if (isEmpty(line))
      return;
    stringstream sin(line);
    pred = succ = 1000;
    sin >> pred >> comma >> succ >> colon;
    sin >> s[pred - 1][succ - 1];
  }
}

void Instance::readSBF(istream &in, bool do_preprocess) {
  try {

    char comma;

    skipTo(in, "number of tasks");
    in >> n;
    allocate();

    skipTo(in, "cycle time");
    in >> c;

    skipTo(in, "task times");
    for (Task i = 0; i != n; ++i) {
      unsigned id;
      in >> id >> t0[i];
    }

    string line;
    unsigned pred, succ;

    skipTo(in, "precedence relations");
    while (true) {
      getline(in, line);
      if (isEmpty(line))
        break;
      stringstream sin(line);
      sin >> pred >> comma >> succ;
      d[pred - 1][succ - 1] = true;
    }

    bool seenForward = false, seenBackward = false;
    string section = "";

    while (section != "end") {
      section = skipToNext(in);
      if (section == "setup times forward") {
        if (seenForward)
          fmt::print(cerr, "Warning: input contains multiple forward sections.\n");
        else
          seenForward = true;
        handleSetupSection(in, sf0);
      } else if (section == "setup times backward") {
        if (seenBackward)
          fmt::print(cerr, "Warning: input contains multiple backward sections.\n");
        else
          seenBackward = true;
        handleSetupSection(in, sb0);
      }
    }

    section = skipToNext(in);
    if (section == "optimal SALBP-1 value")
      in >> optm;

    if (!seenForward)
      if (!seenBackward)

        if (opt.type2)
          c = uc();

    setupDirected = true;
    t = t0;
    sf = sf0;
    sb = sb0;

    preprocess(do_preprocess);
  } catch (string &s) {
    fmt::print(cerr, "Section {}\n", s);
  }
}

void Instance::read(istream &in, string format, bool do_preprocess) {
  if (format == "mp")
    readMP(in, do_preprocess);
  else
    readSBF(in, do_preprocess);
}

int Instance::increaseTimes() {
  unsigned inco = 0, inal = 0;

  vector<Time> ω(n, inf_time);
  for (Task i = 0; i != n; ++i)
    ω[i] = get_ω(i);

  for (Task i = 0; i != n; ++i) {
    if (0 < ω[i] && ω[i] <= c) {
      inco++;
      t[i] += ω[i];
      for (Task j = 0; j != n; ++j) {
        if (sf[i][j] <= c)
          sf[i][j] -= ω[i];
        if (sb[i][j] <= c)
          sb[i][j] -= ω[i];
      }
    }
  }
  vector<Time> α(n, inf_time);
  for (Task i = 0; i != n; ++i)
    α[i] = get_α(i);

  for (Task i = 0; i != n; ++i) {
    if (0 < α[i] && α[i] <= c) {
      inal++;
      t[i] += α[i];
      for (Task h = 0; h != n; ++h) {
        if (sf[h][i] <= c)
          sf[h][i] -= α[i];
        if (sb[h][i] <= c)
          sb[h][i] -= α[i];
      }
    }
  }

  return inco + inal;
}

int Instance::isolatedTasks() {
  unsigned nit = 0;
  for (Task i = 0; i != n; ++i) {
    unsigned j = 0;
    while (j != n) {
      if (isTriangular) {
        if (t[i] + sf[i][j] + t[j] + sb[j][i] <= c)
          break;
        if (t[i] + sf[j][i] + t[j] + sb[i][j] <= c)
          break;
      } else {
        auto idle = c - t[i] - t[j];
        if (t[i] <= idle || sf[j][i] <= idle)
          break;
      }
      ++j;
    }
    if (j == n && t[i] + sb[i][i] <= c) {
      if (t[i] + sb[i][i] < c) {
        nit++;
        t[i] += c - t[i] - sb[i][i];
      }
    }
  }
  return nit;
}

void Instance::preprocess(bool do_preprocess) {
  preprocess_static();
  if (do_preprocess)
    preprocess_dynamic(n);
}

void Instance::preprocess_static() {

  P.resize(n);
  F.resize(n);
  for (Task i = 0; i != n; ++i) {
    P[i].clear();
    F[i].clear();
  }

  for (Task i = 0; i != n; ++i)
    for (Task j = 0; j != n; ++j)
      if (d[i][j]) {
        P[j].push_back(i);
        F[i].push_back(j);
      }

  D = d;
  for (Task g1 = 0; g1 != n; g1++)
    for (Task g2 = 0; g2 != n; g2++)
      for (Task g3 = 0; g3 != n; g3++)
        D[g2][g3] = D[g2][g3] || (D[g2][g1] && D[g1][g3]);

  Ps.resize(n);
  Fs.resize(n);
  Is.resize(n);
  for (Task i = 0; i != n; ++i) {
    Ps[i].clear();
    Fs[i].clear();
    Is[i].clear();
  }
  for (Task i = 0; i != n; ++i)
    for (Task j = 0; j != n; ++j)
      if (D[i][j]) {
        Ps[j].push_back(i);
        Fs[i].push_back(j);
      } else if (!D[j][i] && i < j) {
        Is[i].push_back(j);
        Is[j].push_back(i);
      }
  maxPs = Ps[0].size();
  for (Task i = 1; i != n; ++i)
    maxPs = max(maxPs, unsigned(Ps[i].size()));
  maxFs = Fs[0].size();
  for (Task i = 1; i != n; ++i)
    maxFs = max(maxFs, unsigned(Fs[i].size()));

  O.resize(n);
  for (Task i = 0; i != n; ++i) {
    O[i].resize(n);
    for (Task j = 0; j != n; ++j) {
      O[i][j].clear();
      for (Task k = 0; k != n; k++)
        if ((D[i][k] && D[k][j]) || (D[j][k] && D[k][i]))
          O[i][j].push_back(k);
    }
  }

  if (!disable_preprocessing) {
    for (Task i = 0; i != n; ++i)
      for (Task j = 0; j != n; ++j) {
        if (!canForward(i, j))
          disableForward(i, j);
        if (!canBackward(i, j))
          disableBackward(i, j);
      }
  }
}

void Instance::compute_ta_tn() {
  ta.assign(n, 0);
  tn.assign(n, 0);

  for (Task i = 0; i != n; ++i) {
    for (Task j = 0; j != n; ++j) {
      if (D[j][i])
        ta[i] += t[j];
      if (D[i][j])
        tn[i] += t[j];
    }
  }
}

void Instance::compute_E_T() {
  E.resize(n);
  T.resize(n);
  for (Task i = 0; i != n; ++i)
    E[i] = (ta[i] + t[i] + c - 1) / c;
  for (Task i = 0; i != n; ++i) {
    T[i] = (t[i] + tn[i] + c - 1) / c;
  }
}

int Instance::excludeSetups(unsigned um) {
  unsigned nsb = 0;
  for (Task i = 0; i != n; ++i)
    for (Task j = i + 1; j != n; ++j)
      if (i != j && !canShare(i, j, um)) {
        nsb += disableForward(i, j) + disableForward(j, i);
        nsb += disableBackward(i, j) + disableBackward(j, i);
      } else {
        if (!isTriangular && t[i] + sf[i][j] + t[j] + sb[j][i] > c && t[j] + sf[j][i] + t[i] + sb[i][j] > c) {
          Time lc = inf_time;
          for (Task k = 0; k != n; ++k)
            if (k != i && k != j)
              lc = min(lc_tasks(*this, {i, j, k}), lc);
          if (lc > c) {
            nsb += disableForward(i, j) + disableForward(j, i);
            nsb += disableBackward(i, j) + disableBackward(j, i);
          }
        }
      }
  return nsb;
}

void Instance::preprocess_dynamic(unsigned um) {

  unsigned it = 1, is = 0, ns = 0;
  while (it + is + ns > 0) {

    compute_ta_tn();
    compute_E_T();
    compute_smallest_setups();
    compute_ILB();

    if (disable_preprocessing) {
      for (Task i = 0; i != n; ++i)
        sf[i][i] = c + 1;
      break;
    }
    ns = excludeSetups(um);
    it = increaseTimes();
    is = isolatedTasks();
  }

  check_symmetry();
}

void Instance::check_symmetry() {
  fwdSymmetric = true;
  for (Task i = 0; i != n; ++i) {
    if (sf[i][i] != c + 1)
      throw "Forward self setup time.";

    if (fwdSymmetric) {
      for (Task j = 0; j != n; ++j)
        if (!D[i][j] && !D[j][i] && sf[i][j] != sf[j][i]) {
          fwdSymmetric = false;
        }
    }
  }

  if (!setupDirected)
    bwdSymmetric = fwdSymmetric;
  else {
    bwdSymmetric = true;
    for (Task i = 0; i != n; ++i)
      for (Task j = 0; j != n; ++j)
        if (!D[i][j] && !D[j][i] && sb[i][j] != sb[j][i]) {
          bwdSymmetric = false;
          return;
        }
  }
}

void Instance::reverse() {
  t = t0;
  sf = sf0;
  sb = sb0;

  for (Task i = 0; i != n; ++i)
    for (Task j = i + 1; j != n; ++j) {
      std::swap(d[i][j], d[j][i]);
      std::swap(sf[i][j], sf[j][i]);
      std::swap(sb[i][j], sb[j][i]);
    }
  for (Task i = 0, ie = (n + 1) / 2; i != ie; ++i) {
    std::swap(t[i], t[n - 1 - i]);
    for (Task j = 0, je = (i == n - 1 - i ? n / 2 : n); j != je; ++j) {
      std::swap(d[i][j], d[n - 1 - i][n - 1 - j]);
      std::swap(sf[i][j], sf[n - 1 - i][n - 1 - j]);
      std::swap(sb[i][j], sb[n - 1 - i][n - 1 - j]);
    }
  }
  preprocess();
}

void Instance::compute_ILB() {
  isTriangular = true;
  isStrictTriangular = true;
  insertionlb.assign(n, numeric_limits<int>::max());
  for (Task i = 0; i != n; ++i) {
    insertionlb[i] = t[i] + sb[i][i];
    for (Task j = 0; j != n; ++j) {
      if (i == j)
        continue;

      if (canForwardExt(i, j) && canBackwardExt(j, i)) {
        int insertion = int(sf[i][j] + t[j] + sb[j][i]) - int(sb[i][i]);
        insertionlb[j] = min(insertionlb[j], insertion);
        if (insertion < 0 && isTriangular) {
          isTriangular = false;
        }
        if (sf[i][j] + sb[j][i] < sb[i][i] && isStrictTriangular) {
          isStrictTriangular = false;
        }
      }
      if (canForwardExt(j, i) && canBackwardExt(i, j)) {
        int insertion = int(sf[j][i] + t[j] + sb[i][j]) - int(sb[i][i]);
        insertionlb[j] = min(insertionlb[j], insertion);
        if (insertion < 0 && isTriangular) {
          isTriangular = false;
        }
        if (sf[j][i] + sb[i][j] < sb[i][i] && isStrictTriangular) {
          isStrictTriangular = false;
        }
      }

      for (Task h = 0; h < n; h++) {
        if (h == i || h == j)
          continue;
        if (canForwardExt(i, h) && canForwardExt(h, j) && canForwardExt(i, j)) {
          int insertion = int(sf[i][h] + t[h] + sf[h][j]) - int(sf[i][j]);
          insertionlb[h] = min(insertionlb[h], insertion);
          if (insertion < 0 && isTriangular) {
            isTriangular = false;
          }
          if (sf[i][h] + sf[h][j] < sf[i][j] && isStrictTriangular) {
            isStrictTriangular = false;
          }
        }
        if (canForwardExt(i, h) && canBackwardExt(h, j) && canBackwardExt(i, j)) {
          int insertion = int(sf[i][h] + t[h] + sb[h][j]) - int(sb[i][j]);
          insertionlb[h] = min(insertionlb[h], insertion);
          if (insertion < 0 && isTriangular) {
            isTriangular = false;
          }
          if (sf[i][h] + sb[h][j] < sb[i][j] && isStrictTriangular) {
            isStrictTriangular = false;
          }
        }
      }
    }
  }
}

void Instance::compute_smallest_setups() {
  sfi.clear();
  unsigned nfwd = 0;
  for (Task i = 0; i != n; ++i) {
    Time sfmin = inf_time;
    for (Task j = 0; j != n; ++j)
      if (i != j) {
        sfmin = min(sfmin, sf[i][j]);
        if (sf[i][j] <= c)
          nfwd++;
      }
    if (sfmin <= c)
      sfi.push_back(sfmin);
  }

  sort(sfi.begin(), sfi.end());

  sbi.clear();
  unsigned nbwd = 0;
  for (Task i = 0; i != n; ++i) {
    Time sbmin = inf_time;
    for (Task j = 0; j != n; ++j) {
      sbmin = min(sbmin, sb[i][j]);
      if (sb[i][j] <= c)
        nbwd++;
    }
    if (sbmin <= c)
      sbi.push_back(sbmin);
  }

  sort(sbi.begin(), sbi.end());
}

void Instance::writeSBF(ostream &out) {
  fmt::print(out, "<number of tasks>\n");
  fmt::print(out, "{}\n\n", n);

  fmt::print(out, "<cycle time>\n");
  fmt::print(out, "{}\n\n", c);

  fmt::print(out, "<task times>\n");
  for (Task i = 0; i != n; ++i)
    fmt::print(out, "{} {}\n", i + 1, t[i]);
  fmt::print(out, "\n");

  fmt::print(out, "<precedence relations>\n");
  for (Task i = 0; i != n; ++i)
    for (Task j = 0; j != n; ++j)
      if (d[i][j])
        fmt::print(out, "{},{}\n", i + 1, j + 1);
  fmt::print(out, "\n");

  fmt::print(out, "<setup times forward>\n");
  for (Task i = 0; i != n; ++i)
    for (Task j = 0; j != n; ++j)
      if (sf[i][j] > 0)
        fmt::print(out, "{},{}:{}\n", i + 1, j + 1, sf[i][j]);
  fmt::print(out, "\n");

  fmt::print(out, "<setup times backward>\n");
  for (Task i = 0; i != n; ++i)
    for (Task j = 0; j != n; ++j)
      if (sb[i][j] > 0)
        fmt::print(out, "{},{}:{}\n", i + 1, j + 1, sb[i][j]);
  fmt::print(out, "\n");

  fmt::print(out, "<end>\n");

  if (optm > 0) {
    fmt::print(out, "<optimal SALBP-1 value>\n");
    fmt::print(out, "{}\n\n", optm);
  }
}
double Instance::avgProcessingTime() const {
  Time T = 0;

  for (Task i = 0; i != n; ++i)
    T += t[i];

  return double(T) / n;
}

double Instance::avgForwardSetupTime() const {
  Time T = 0;
  unsigned nfwd = 0;
  for (Task i = 0; i != n; ++i) {
    for (Task j = 0; j != n; ++j)
      if (i != j && sf[i][j] <= c) {
        T += sf[i][j];
        nfwd++;
      }
  }
  return double(T) / nfwd;
}

double Instance::avgBackwardSetupTime() const {
  Time T = 0;
  unsigned nback = 0;
  for (Task i = 0; i != n; ++i) {
    for (Task j = 0; j != n; ++j)
      if (i != j && sb[i][j] <= c) {
        T += sb[i][j];
        nback++;
      }
  }
  return double(T) / nback;
}

double Instance::orderStrength() const {
  unsigned prec = 0;
  for (Task i = 0; i != n; ++i)
    for (Task j = 0; j != n; ++j)
      if (D[i][j])
        prec++;
  return double(prec) / (n * (n - 1) / 2);
}

Time Instance::get_ωτ(Task i, const BitVector &valid) const {
  Time ω = inf_time;
  for (Task j = 0; j != n; ++j)
    if (valid[j] && canForwardExt(i, j))
      ω = min(ω, sf[i][j]);
  return ω;
}

Time Instance::get_ωμ(Task i, const BitVector &valid) const {
  Time ω = inf_time;
  for (Task j = 0; j != n; ++j)
    if (valid[j] && canBackwardExt(i, j))
      ω = min(ω, sb[i][j]);
  return ω;
}

Time Instance::get_ατ(Task i, const BitVector &valid) const {
  Time α = inf_time;
  for (Task h = 0; h != n; ++h)
    if (valid[h] && canForwardExt(h, i))
      α = min(α, sf[h][i]);
  return α;
}

Time Instance::get_αμ(Task i, const BitVector &valid) const {
  Time α = inf_time;
  for (Task h = 0; h != n; ++h)
    if (valid[h] && canBackwardExt(h, i))
      α = min(α, sb[h][i]);
  return α;
}

bool Instance::canShare(Task i, Task j, unsigned um) const {
  if (um + 1 - T[i] < E[j] || um + 1 - T[j] < E[i])
    return false;
  if (lc_tasks(*this, {i, j}) > c)
    return false;
  return true;
}

bool Instance::canAssignFwd(Task i, Task j, unsigned k, unsigned um) const {
  Time ωhead = t[i] + sf[i][j] + t[j], ωtail = ωhead;
  for (Task k = 0; k != n; ++k) {
    if ((D[k][i] || D[k][j]) && k != i && k != j)
      ωhead += t[k];
    if ((D[i][k] || D[j][k]) && k != i && k != j)
      ωtail += t[k];
  }
  return (ωhead + c - 1) / c <= k && k <= um + 1 - (ωtail + c - 1) / c;
}

bool Instance::canAssignBwd(Task i, Task j, unsigned k, unsigned um) const {
  if (i == j)
    return true;
  Time αhead = t[i] + sb[i][j] + t[j], αtail = αhead;
  for (Task k = 0; k != n; ++k) {
    if ((D[k][i] || D[k][j]) && k != i && k != j)
      αhead += t[k];
    if ((D[i][k] || D[j][k]) && k != i && k != j)
      αtail += t[k];
  }
  return (αhead + c - 1) / c <= k && k <= um + 1 - (αtail + c - 1) / c;
}
} // namespace sualbsp

string guess_format(string instance) {
  if (fs::path(instance).extension() == ".alb")
    return "sbf";
  else
    return "mp";
}
