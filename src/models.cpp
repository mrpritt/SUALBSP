#include "models.hpp"

using namespace std;

#include "logging.hpp"
#include "prbp.hpp"

void setupSolver(IloCplex solver, const ModelOptions &m) {
  IloEnv env = solver.getEnv();
  solver.setOut(env.getNullStream());
  solver.setWarning(env.getNullStream());
  solver.setParam(IloCplex::Param::Threads, m.threads);
  solver.setParam(IloCplex::Param::TimeLimit, m.timelimit);
}

IloCplex solveRelaxed(IloModel model, const vector<IloNumVarArray> &relax, const ModelOptions &m) {
  IloEnv env = model.getEnv();
  IloModel relaxed(env);
  relaxed.add(model);
  for (auto r : relax)
    relaxed.add(IloConversion(env, r, ILOFLOAT));

  IloCplex solver(relaxed);
  setupSolver(solver, m);
  if (m.wrelaxation != "/dev/null")
    solver.exportModel(m.wrelaxation.c_str());
  solver.solve();

  return solver;

}

void SSBF::addVars() {
  for (unsigned i = 0; i != I.n; i++) {
    for (auto [k, ke] = I.FS(i, um); k != ke; ++k) {
      x.add(IloNumVar(env, 0.0, 1.0, ILOINT));
      x[x.getSize() - 1].setName(fmt::format("x[{},{}]", i + 1, k).c_str());
      ix[{i, k}] = x.getSize() - 1;
    }
    z.add(IloNumVar(env, 1.0, IloInfinity, ILOFLOAT)); 
    z[z.getSize() - 1].setName(fmt::format("z[{}]", i + 1).c_str());
  }
  nx = ix.size();

  for (auto k = lm + 1; k <= um; ++k) {
    u.add(IloNumVar(env, 0.0, 1.0, ILOINT));
    u[u.getSize() - 1].setName(fmt::format("u[{}]", k).c_str());
    obj.setLinearCoef(u[u.getSize() - 1], 1.0);
  }

  for (unsigned i = 0; i != I.n; ++i)
    for (auto [k, ke] = I.FS(i, um); k != ke; ++k)
      for (unsigned j = 0; j != I.n; ++j)
        if (I.canForwardExt(i, j) && I.canAssign(j, k, um) && I.canAssignFwd(i, j, k, um)) {
          g.add(IloNumVar(env, 0.0, 1.0, ILOINT));
          g[g.getSize() - 1].setName(fmt::format("g[{},{},{}]", i + 1, j + 1, k).c_str());
          ig[{i, j, k}] = g.getSize() - 1;
        }
  ng = ig.size();

  for (unsigned i = 0; i != I.n; ++i)
    for (auto [k, ke] = I.FS(i, um); k != ke; ++k)
      for (unsigned j = 0; j != I.n; ++j)
        if (I.canBackwardExt(i, j) && I.canAssign(j, k, um) && I.canAssignBwd(i, j, k, um)) {
          h.add(IloNumVar(env, 0.0, 1.0, ILOINT));
          h[h.getSize() - 1].setName(fmt::format("h[{},{},{}]", i + 1, j + 1, k).c_str());
          ih[{i, j, k}] = h.getSize() - 1;
        }
  nh = ih.size();

  for (unsigned i = 0; i != I.n; ++i) {
    r.add(IloNumVar(env, I.Ps[i].size() + 1, I.n - I.Fs[i].size(), ILOFLOAT));
    r[r.getSize() - 1].setName(fmt::format("r[{}]", i + 1).c_str());
  }
}

void SSBF::addAssignmentStations() {
  IloRangeArray assign(env);
  IloRangeArray stations(env);

  for (unsigned i = 0; i != I.n; ++i) {
    IloExpr assign_i(env);
    IloExpr station = -z[i];

    for (auto [k, ke] = I.FS(i, um); k != ke; ++k) {
      assign_i += x[ix[{i, k}]];
      station += double(k) * x[ix[{i, k}]];
    }
    assign.add(IloRange(env, 1.0, assign_i, 1.0, fmt::format("a{}", i + 1).c_str()));
    stations.add(IloRange(env, 0.0, station, 0.0, fmt::format("z{}", i + 1).c_str()));
  }
  model.add(assign);
  model.add(stations);
  assign.end();
  stations.end();
}

void SSBF::addPrefix() {
  for (auto k = lm + 1; k < um; ++k)
    model.add(IloRange(env, -IloInfinity, u[k - lm] - u[k - lm - 1], 0.0, fmt::format("u{}", k).c_str()));
}

void SSBF::addSequence() {
  IloRangeArray forward(env);

  for (unsigned i = 0; i != I.n; ++i)
    for (auto [k, ke] = I.FS(i, um); k != ke; ++k) {
      IloExpr succ = -x[ix[{i, k}]];

      for (unsigned j = 0; j != I.n; ++j) {
        if (I.canForwardExt(i, j) && I.canAssign(j, k, um) && I.canAssignFwd(i, j, k, um))
          succ += g[ig[{i, j, k}]];
        if (I.canBackwardExt(i, j) && I.canAssign(j, k, um) && I.canAssignBwd(i, j, k, um))
          succ += h[ih[{i, j, k}]];
      }
      forward.add(IloRange(env, 0.0, succ, 0.0, fmt::format("fwd{},{}", i + 1, k).c_str()));
      succ.end();
    }
  model.add(forward);
  forward.end();

  IloRangeArray backward(env);

  for (unsigned i = 0; i != I.n; ++i)
    for (auto [k, ke] = I.FS(i, um); k != ke; ++k) {
      IloExpr pred = -x[ix[{i, k}]];

      for (unsigned j = 0; j != I.n; ++j) {
        if (I.canForwardExt(j, i) && I.canAssign(j, k, um) && I.canAssignFwd(j, i, k, um))
          pred += g[ig[{j, i, k}]];
        if (I.canBackwardExt(j, i) && I.canAssign(j, k, um) && I.canAssignBwd(j, i, k, um))
          pred += h[ih[{j, i, k}]];
      }

      backward.add(IloRange(env, 0.0, pred, 0.0, fmt::format("bwd{},{}", i + 1, k).c_str()));
      pred.end();
    }
  model.add(backward);
  backward.end();
}

void SSBF::addStationsOcurrence() {
  IloRangeArray ocurrence(env);

  IloExprArray last_first(env, um + 1);

  for (auto k = 1u; k <= um; ++k) {
    last_first[k] = IloExpr(env);

    for (unsigned i = 0; i != I.n; ++i) {
      if (!I.canAssign(i, k, um))
        continue;
      for (unsigned j = 0; j != I.n; ++j) {
        if (!I.canAssign(j, k, um))
          continue;
        if (I.canBackwardExt(i, j) && I.canAssignBwd(i, j, k, um))
          last_first[k] += h[ih[{i, j, k}]];
      }
    }

    if (k <= lm) {
      ocurrence.add(IloRange(env, 1.0, last_first[k], 1.0, fmt::format("ol{}", k).c_str()));
    } else {
      last_first[k] -= u[k - lm - 1];
      ocurrence.add(IloRange(env, 0.0, last_first[k], 0.0, fmt::format("ou{}", k).c_str()));
    }
  }
  model.add(ocurrence);
  ocurrence.end();
}

void SSBF::addPrecedences() {
  IloRangeArray precedence(env);
  int back = -1;

  for (unsigned i = 0; i != I.n; ++i)
    for (unsigned j = 0; j != I.n; ++j) {
      if (I.canForwardExt(i, j) && !I.d[i][j]) {
        const auto max_between = I.n - I.Fs[i].size() - I.Ps[j].size();

        precedence.add(IloRange(env, -IloInfinity, max_between - 1.0, fmt::format("prec1-{},{}", i + 1, j + 1).c_str()));
        back++;
        precedence[back].setLinearCoef(r[i], 1.0);
        precedence[back].setLinearCoef(r[j], -1.0);

        auto [fi, li] = I.FS(i, um);
        auto [fj, lj] = I.FS(j, um);

        for (auto k = max(fi, fj), ke = min(li, lj); k != ke; ++k)
          if (I.canAssignFwd(i, j, k, um))
            precedence[back].setLinearCoef(g[ig[{i, j, k}]], max_between);
      }
      if (I.d[i][j]) {
        model.add(IloRange(env, -IloInfinity, r[i] - r[j], -1.0, fmt::format("prec2-{},{}", i + 1, j + 1).c_str()));
        model.add(IloRange(env, -IloInfinity, z[i] - z[j], (I.canShare(i, j, um) ? 0.0 : -1.0), fmt::format("prec3-{},{}", i, j).c_str())); 
      }
    }
  model.add(precedence);
  precedence.end();
}

void SSBF::addCycle() {
  IloRangeArray cycle(env);

  IloExprArray ctime(env, um + 1);

  for (auto k = 1u; k <= um; ++k) {
    ctime[k] = IloExpr(env);

    for (unsigned i = 0; i != I.n; ++i) {
      if (!I.canAssign(i, k, um))
        continue;
      ctime[k] += double(I.t[i]) * x[ix[{i, k}]];
      for (unsigned j = 0; j != I.n; ++j) {
        if (!I.canAssign(j, k, um))
          continue;
        if (I.canForwardExt(i, j) && I.canAssignFwd(i, j, k, um))
          ctime[k] += double(I.sf[i][j]) * g[ig[{i, j, k}]];
        if (I.canBackwardExt(i, j) && I.canAssignBwd(i, j, k, um))
          ctime[k] += double(I.sb[i][j]) * h[ih[{i, j, k}]];
      }
    }

    if (k <= lm)
      cycle.add(IloRange(env, -IloInfinity, ctime[k], I.c, fmt::format("cycle1-{}", k).c_str()));
    else {
      ctime[k] += -double(I.c) * u[k - lm - 1];
      cycle.add(IloRange(env, -IloInfinity, ctime[k], 0.0, fmt::format("cycle2-{}", k).c_str()));
    }
    ctime[k].end();
  }
  model.add(cycle);
  cycle.end();
}

void SSBF::addSingleTask() {
  IloRangeArray single(env);
  int back = -1;

  const auto max_tasks = I.n - lm + 1;

  for (auto k = 1u; k <= um; ++k)
    for (unsigned j = 0; j != I.n; ++j)
      if (I.canAssign(j, k, um)) {
        single.add(IloRange(env, -IloInfinity, max_tasks, fmt::format("single{},{}", k, j + 1).c_str()));
        back++;
        single[back].setLinearCoef(h[ih[{j, j, k}]], max_tasks);
        for (unsigned i = 0; i != I.n; ++i)
          if (i != j && I.canAssign(i, k, um))
            single[back].setLinearCoef(x[ix[{i, k}]], 1.0);
      }
  model.add(single);
  single.end();
}

void SSBF::build() {

  addVars();
  addAssignmentStations();
  addPrefix();
  addSequence();
  addStationsOcurrence();
  addPrecedences();
  addCycle();
  addSingleTask();

}

double SSBF::solveRelaxed() {
  IloCplex solver = ::solveRelaxed(model, {x, u, g, h}, m);

  IloAlgorithm::Status status = solver.getStatus();
  if (status == IloAlgorithm::Feasible || status == IloAlgorithm::Optimal) {
    fmt::print("LP-SSBF: {}\n", lm + solver.getObjValue());

    return lm + solver.getObjValue();
  } else
    throw "Relaxation of SSBF model is infeasible";
}

std::pair<Solution, IloAlgorithm::Status> SSBF::solve(const Solution &S) {
  IloCplex solver(model);
  setupSolver(solver, m);

  if (S.isValid())
    setSolution(solver, S);

  solver.solve();

  Solution T;
  IloAlgorithm::Status status = solver.getStatus();
  if (status == IloAlgorithm::Feasible || status == IloAlgorithm::Optimal)
    T = getSolution(solver);
  else
  solver.end();
  return {T, status};
}

void SSBF::setSolution(IloCplex cplex, const Solution &S) {

  vector<short unsigned> a;
  assign(I, S.pi, a);

  IloNumVarArray startVar(env);
  IloNumArray startVal(env);

  unsigned cm = 0, ft = I.n, pt = I.n;
  for (unsigned i = 0; i != I.n; ++i) {
    unsigned m = a[i];
    Task ct = S.pi[i];
    if (m != cm) {
      startVar.add(h[ih[{ft, pt, cm}]]);
      startVal.add(1.0);
      cm = m;
      ft = ct;
    } else {
      startVar.add(g[ig[{pt, ct, cm}]]);
      startVal.add(1.0);
    }
    startVar.add(x[ix[{ct, cm}]]);
    startVal.add(1.0);
    pt = ct;
  }
  cplex.addMIPStart(startVar, startVal);
  startVar.end();
  startVal.end();
}

Solution SSBF::getSolution(IloCplex solver) {

  IloNumArray xvalue(env);
  solver.getValues(x, xvalue);
  vector<unsigned> a(I.n, 0);
  for (auto [el, idx] : ix) {
    auto i = get<0>(el);
    auto k = get<1>(el);
    if (xvalue[idx] > 0.5) {
      a[i] = k;
    }
  }

  IloNumArray gvalue(env);
  solver.getValues(g, gvalue);
  vector<int> p(I.n, -1), s(I.n, -1);
  for (auto [el, idx] : ig) {
    auto i = get<0>(el);
    auto j = get<1>(el);
    [[maybe_unused]] auto k = get<2>(el);
    if (gvalue[idx] > 0.5) {
      p[j] = i;
      s[i] = j;
    }
  }

  Solution S;

  unsigned k = 1;
  while (true) {
    auto ktask = find(a.begin(), a.end(), k);
    if (ktask == a.end())
      break;
    Task i = *ktask;
    while (p[i] != -1)
      i = p[i];
    while (i != -1) {
      S.pi.push_back(i);
      i = s[i];
    }
    k++;
  }
  S.m = k - 1;
  return S;
}
