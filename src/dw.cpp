#include <bit>
#include <cassert>
#include <iostream>
#include <string>
using namespace std;

#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

#include "dw.hpp"
#include "hoffmann.hpp"
#include "options.hpp"
#include "solution.hpp"
#include "sop.hpp"
#include "util.hpp"

#define NUM_ITERS_TO_TEST 100
#define RC_EPS 0.000001

void dw_t::clear(unsigned n) {
  tasks.clear();

  order.assign(n, 0);
  iota(order.begin(), order.end(), 0);

  bestSols.clear();
  state.assign(n, UNKNOWN);

  pi.clear();

  obj = 0.0;
  bestObj = 1.0;

  itCount = tAcc = 0;
}

vector<unsigned> support(const vector<bool> &v) {
  vector<unsigned> s;
  for (auto i = 0ul, ie = v.size(); i != ie; ++i)
    if (v[i])
      s.push_back(i);
  return s;
}

void support(const vector<short int> &v, vector<Task> &s) {
  for (auto i = 0ul, ie = v.size(); i != ie; ++i)
    if (v[i] == IN)
      s.push_back(i);
}

double solveEnhancedKnapsack(const Instance &I, const vector<double> &price, const ModelOptions &mopt) {
  IloEnv env;
  IloModel model(env);
  IloNumVarArray v(env, I.n, 0, 1, ILOINT);
  IloNumVarArray b(env, I.n, 0, 1, ILOINT);
  IloNumVarArray f(env, I.n, 0, 1, ILOINT);
  IloExpr objective(env);

  for (auto i = 0u; i != I.n; ++i)
    objective += price[i] * v[i];
  model.add(IloMaximize(env, objective));
  objective.end();

  IloExpr capacity(env);
  for (auto i = 0u; i != I.n; ++i) {
    capacity += int(I.t[i]) * v[i];
    if (I.isTriangular) {
      Time minτ = I.c;
      for (auto j = 0u; j != I.n; ++j)
        minτ = min(minτ, I.sf[i][j]);
      capacity += int(minτ) * f[i];

      Time minμ = I.c;
      for (auto j = 0u; j != I.n; ++j)
        minμ = min(minμ, I.sb[i][j]);
      capacity += int(minμ) * b[i];
    }
  }
  model.add(capacity <= int(I.c));
  capacity.end();

  if (I.isTriangular) {
    for (auto i = 0u; i != I.n; ++i)
      model.add(b[i] + f[i] >= v[i]);

    IloExpr oneback(env);
    for (auto i = 0u; i != I.n; ++i)
      oneback += b[i];

    model.add(oneback == 1);
    oneback.end();
  }
  IloCplex solver(model);
  setupSolver(solver, mopt);

  solver.solve();
  auto objValue = solver.getObjValue();
  env.end();
  return objValue - 1.0;
}

void dw_t::allocate_memory(unsigned n) {
  assigncost = new cost *[n];
  for (auto i = 0u; i < n; i++)
    assigncost[i] = new cost[n];
  rowsol = new col[n];
  colsol = new row[n];
  u = new cost[n];
  v = new cost[n];
  map = new Task[n];
}

void dw_t::free_memory(unsigned n) {
  for (auto i = 0u; i < n; i++)
    delete[] assigncost[i];
  delete[] assigncost;
  delete[] rowsol;
  delete[] colsol;
  delete[] u;
  delete[] v;
  delete[] map;
}

bool checkFeasibility(const Instance *I, dw_t *d) {
  vector<Task> whichTasks;
  vector<Task> initialTasks;

  whichTasks.reserve(bit_ceil(d->pi.size()));
  support(d->state, whichTasks);

  for (auto ii = 0ul, ie = whichTasks.size(); ii != ie; ++ii) {
    Task i = whichTasks[ii];
    for (Task j : whichTasks)
      if (I->D[j][i])
        goto nextTask;
    initialTasks.push_back(i);
  nextTask:;
  }

  if (whichTasks.size() == 1)
    if (I->sb[whichTasks[0]][whichTasks[0]] + d->tAcc <= I->c)
      return true;

  for (Task i = 0; i < whichTasks.size(); i++)
    d->map[i] = whichTasks[i];
  for (Task i = 0; i < whichTasks.size(); i++)
    for (Task j = 0; j < whichTasks.size(); j++)
      d->assigncost[i][j] = i == j ? BIG : I->sf[whichTasks[i]][whichTasks[j]];

  for (Task initial = 0; initial < initialTasks.size(); initial++) {
    for (Task i = 0; i < whichTasks.size(); i++)
      d->assigncost[i][initial] = i == initial ? BIG : I->sb[whichTasks[i]][initialTasks[initial]];

    cost extraTime;
    extraTime = lap(whichTasks.size(), d->assigncost, d->rowsol, d->colsol, d->u, d->v);
    if (extraTime + d->tAcc <= I->c)
      return true;
    for (Task i = 0; i < whichTasks.size(); i++)
      d->assigncost[i][initial] = i == initial ? BIG : I->sf[whichTasks[i]][initialTasks[initial]];
  }

  return false;
}

double boundPricing(const Instance *I, dw_t *d, int pos) {
  Time idle = I->c - d->tAcc;
  double bound = d->obj;
  Task id;
  if ((opt.dwidlebound == true) && (I->isTriangular)) {
    vector<Task> whichTasks;
    whichTasks.reserve(bit_ceil(d->pi.size()));
    support(d->state, whichTasks);
    vector<unsigned> minForward(whichTasks.size(), I->c);
    vector<unsigned> minBackward(whichTasks.size(), I->c);
    if (whichTasks.size() == 1) {
      minBackward[0] = I->sb[whichTasks[0]][whichTasks[0]];
    }
    for (auto pos = 0u; pos < whichTasks.size(); pos++) {
      auto task = whichTasks[pos];
      for (auto i : I->Ps[task]) {
        short int valid = d->state[i];
        if (valid == UNKNOWN) {
          valid = IN;
          if (I->t[i] > idle)
            valid = OUT;
          if (I->E[task] > d->lSt)
            valid = OUT;
          if ((incumbent.m - I->T[task]) < d->eSt)
            valid = OUT;
        }
        if (valid == IN) {
          if (minBackward[pos] > I->sb[task][i])
            minBackward[pos] = I->sb[task][i];
        }
      }
      for (auto i : I->Fs[task]) {
        short int valid = d->state[i];
        if (valid == UNKNOWN) {
          valid = IN;
          if (I->t[i] > idle)
            valid = OUT;
          if (I->E[task] > d->lSt)
            valid = OUT;
          if ((incumbent.m - I->T[task]) < d->eSt)
            valid = OUT;
        }
        if (valid == IN) {
          if (minForward[pos] > I->sf[task][i])
            minForward[pos] = I->sf[task][i];
        }
      }
      for (auto i : I->Is[task]) {
        short int valid = d->state[i];
        if (valid == UNKNOWN) {
          valid = IN;
          if (I->t[i] > idle)
            valid = OUT;
          if (I->E[task] > d->lSt)
            valid = OUT;
          if ((incumbent.m - I->T[task]) < d->eSt)
            valid = OUT;
        }
        if (valid == IN) {
          if (minBackward[pos] > I->sb[task][i])
            minBackward[pos] = I->sb[task][i];
          if (minForward[pos] > I->sf[task][i])
            minForward[pos] = I->sf[task][i];
        }
      }
    }
    auto sumSetup = accumulate(minForward.begin(), minForward.end(), 0u);
    auto minSetup = I->c;
    for (auto i = 0ul; i < minForward.size(); i++) {
      auto currentSetup = sumSetup - minForward[i] + minBackward[i];
      if (currentSetup < minSetup)
        minSetup = currentSetup;
    }
    if (minSetup >= idle) {
      return bound;
    }
    idle = idle - minSetup;
  }

  Time maxTaskTime = idle;
  for (Task i = pos + 1; i < d->tasks.size(); i++) {
    id = d->tasks[i].id;
    bool stCompatible = true;
    if (I->E[id] > d->lSt)
      stCompatible = false;
    if ((incumbent.m - I->T[id]) < d->eSt)
      stCompatible = false;
    if (d->state[id] == UNKNOWN && I->t[id] <= maxTaskTime && stCompatible) {
      if (I->t[id] > idle) {
        bound += d->tasks[i].priceRatio * idle;
        idle = 0;
        break;
      } else {
        bound += d->tasks[i].price;
        idle -= I->t[id];
      }
    }
  }
  return bound;
}

bool canAddCompact(const Instance *I, dw_t *d, int pos, vector<Task> &load) {

  d->state[d->tasks[pos].id] = IN;
  if (d->eSt < I->E[d->tasks[pos].id])
    d->eSt = I->E[d->tasks[pos].id];
  if (d->lSt > incumbent.m - I->T[d->tasks[pos].id])
    d->lSt = incumbent.m - I->T[d->tasks[pos].id];
  d->obj += d->tasks[pos].price;
  d->tAcc += I->t[d->tasks[pos].id];

  if (d->tAcc > I->c)
    return false;
  if (d->eSt > d->lSt)
    return false;

  vector<short int> state(d->state);
  load.push_back(d->tasks[pos].id);
  unsigned k = load.size();

  for (Task i = 0; i < k; ++i) {
    Task ii = load[i];
    for (Task j = i + 1; j < k; ++j) {
      Task jj = load[j];

      Task t1 = jj, t2 = ii;
      if (ii < jj)
        swap(t1, t2);

      for (Task k : I->O[t1][t2]) {
        if (d->order[k] < pos) {
          if (state[k] == OUT) {
            load.pop_back();
            return false;
          }
        } else if (d->order[k] > pos) {
          if (d->state[k] == UNKNOWN) {
            d->state[k] = IN;
            d->obj += d->tasks[d->order[k]].price;
            d->tAcc += I->t[k];
            if (d->tAcc > I->c) {
              load.pop_back();
              return false;
            }
          }
        }
      }
    }
  }
  load.pop_back();
  return true;
}

bool branchPricing(const Instance *I, dw_t *d, int pos) {
  if (d->itCount % NUM_ITERS_TO_TEST == 0 && d->bestObj > 1.1) {
    return false;
  }
  if (logging::elapsed() > opt.maxtimemodel)
    throw "timeout";

  int tAcc = d->tAcc;
  double obj = d->obj;
  vector<short int> state(d->state);
  auto eSt = d->eSt;
  auto lSt = d->lSt;

  d->itCount++;

  vector<Task> load;
  load.reserve(bit_ceil(d->pi.size()));
  support(d->state, load);

  for (Task i = pos + 1; i < I->n; i++) {
    const auto t = d->tasks[i];
    d->tAcc = tAcc;
    d->obj = obj;
    d->state = state;
    d->eSt = eSt;
    d->lSt = lSt;

    if (d->state[t.id] == UNKNOWN) {
      if (canAddCompact(I, d, i, load)) {
        if (boundPricing(I, d, i) > d->bestObj + RC_EPS) {
          if (d->obj > d->bestObj) {
            if (opt.dwtype == "salbp") {
              d->bestObj = d->obj;
              d->bestSols.push_back(d->state);
            } else {
              bool ls = checkFeasibility(I, d);
              if (ls) {
                d->bestObj = d->obj;
                d->bestSols.push_back(d->state);
              }
            }
          }
          d->pi.push_back(t.id);
          state[t.id] = IN;
          d->state[t.id] = IN;
          bool opt = branchPricing(I, d, i);
          d->state[t.id] = UNKNOWN;
          state[t.id] = OUT;
          d->pi.pop_back();
          if (!opt)
            return false;
        } else {
          state[t.id] = OUT;
        }
      } else {
        state[t.id] = OUT;
      }
    } else {
    }
  }
  return true;
}

double pricingProblem(const Instance *I, vector<double> price, bool &isOpt, double &bound, dw_t *d) {
  isOpt = true;

  d->clear(I->n);

  for (Task i = 0; i < I->n; i++)
    d->tasks.push_back({i, price[i], price[i] / double(I->t[i])});
  sort(d->tasks.begin(), d->tasks.end(), [](auto const &a, auto const &b) { return a.priceRatio > b.priceRatio; });
  for (auto i = 0u; i < I->n; i++)
    d->order[d->tasks[i].id] = i;

  for (Task i = 0; i < I->n; i++) {
    const auto t = d->tasks[i];
    d->eSt = I->E[t.id];
    d->lSt = incumbent.m - I->T[t.id];
    d->pi.push_back(t.id);
    d->obj = t.price;
    d->tAcc = I->t[t.id];
    d->state[t.id] = IN;
    if (d->obj > d->bestObj) {
      d->bestObj = d->obj;
      d->bestSols.push_back(d->state);
    }
    if (boundPricing(I, d, i) > d->bestObj + RC_EPS) {
      isOpt = branchPricing(I, d, i);
    } else

      d->tAcc = 0;
    d->pi.pop_back();
    d->state[t.id] = OUT;
    for (Task pos = i + 1; pos < I->n; pos++)
      d->state[d->tasks[pos].id] = UNKNOWN;

    if (!isOpt)
      break;
  }
  if (d->bestObj > 1.0 + RC_EPS) {
    if (isOpt) {
      bound = d->bestObj - 1.0;
    } else {
      Time idle = I->c;
      bound = -1.0;
      for (Task i = 0; i < I->n; i++) {
        const auto t = d->tasks[i];
        if (I->t[t.id] >= idle) {
          idle -= I->t[t.id];
          bound += t.price;
        } else {
          bound += t.priceRatio * idle;
          break;
        }
      }
    }
  } else
    bound = 0.0;

  return d->bestObj;
}

unsigned short int dw(Instance &I, const ModelOptions &mopt) {
  opt.maxtimemodel += logging::elapsed();
  int num_pricings = 0;
  bool optimalPricing;
  double bestColumn, boundValue, objValue = 0.0;
  vector<double> price(I.n);
  dw_t d;
  unsigned currentlb = 0;

  d.allocate_memory(I.n);

  IloEnv env;
  IloModel m(env);
  IloObjective obj = IloAdd(m, IloMinimize(env));

  IloNumArray RHS = IloNumArray(env, I.n);
  for (auto i = 0u; i < I.n; i++)
    RHS[i] = 1.0;
  IloRangeArray constraints(env, RHS, IloInfinity);
  m.add(constraints);
  IloNumVarArray v(env);
  IloNumArray newColumn(env, I.n);
  for (auto i = 0u; i < I.n; i++) {
    for (auto j = 0u; j < I.n; j++)
      newColumn[j] = i == j;
    v.add(IloNumVar(obj(1) + constraints(newColumn), 0, IloInfinity, ILOFLOAT));
  }
  if (incumbent.m < I.n) {
    for (auto i = 0u; i < incumbent.m; i++) {
      for (auto j = 0u; j < I.n; j++)
        newColumn[incumbent.pi[j]] = i + 1 == incumbent.a[j];
      if (count(incumbent.a.begin(), incumbent.a.end(), i + 1) > 1)
        v.add(IloNumVar(obj(1) + constraints(newColumn), 0, IloInfinity, ILOFLOAT));
    }
  }

  IloCplex solver(m);
  if (mopt.wrelaxation != "/dev/null")
    solver.exportModel(mopt.wrelaxation.c_str());
  setupSolver(solver, mopt);

  bool early_stop = true;

  while (logging::elapsed() < opt.maxtimemodel) {
    num_pricings++;
    solver.solve();
    for (auto i = 0u; i < I.n; i++) {
      price[i] = solver.getDual(constraints[i]);
    }

    objValue = solver.getObjValue();

    optimalPricing = false;
    try {
      bestColumn = pricingProblem(&I, price, optimalPricing, boundValue, &d);
    } catch (...) {
      break;
    }
    if (bestColumn > 1.0 + RC_EPS) {
      for (auto bestSol : d.bestSols) {
        for (auto j = 0u; j < I.n; j++)
          newColumn[j] = bestSol[j] == IN;
        v.add(IloNumVar(obj(1) + constraints(newColumn), 0, 1, ILOFLOAT));
      }
    } else {
      early_stop = false;
      break;
    }
    if (abs(ceil(objValue - RC_EPS) - ceil(objValue / (1 + boundValue) - RC_EPS)) < RC_EPS) {
      currentlb = ceil(objValue - RC_EPS);
      early_stop = false;
      break;
    }
    if (currentlb < ceil(objValue / (1 + boundValue) - RC_EPS)) {
      currentlb = ceil(objValue / (1 + boundValue) - RC_EPS);
    }
  }

  d.free_memory(I.n);

  unsigned lb = ceil(objValue - RC_EPS);

  if (early_stop) {
    solver.solve();
    for (auto i = 0u; i < I.n; i++)
      price[i] = solver.getDual(constraints[i]);
    auto reducedCost = solveEnhancedKnapsack(I, price, mopt);
    lb = ceil(objValue / (1.0 + reducedCost) - RC_EPS);
  }
  if (lb < currentlb)
    lb = currentlb;

  env.end();
  if (incumbent.m < lb) {
    lb = incumbent.m;
  }
  return lb;
}
