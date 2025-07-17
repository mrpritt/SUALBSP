#pragma once

#include <map>
#include <string>
#include <tuple>

#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

#include "instance.hpp"
#include "solution.hpp"

IloCplex solveRelaxed(IloModel model, IloNumVarArray x);

struct ModelOptions {
  unsigned threads;
  unsigned timelimit;
  bool verbose;
  std::string wrelaxation;

  ModelOptions() : threads(1), timelimit(30), verbose(false), wrelaxation("/dev/null") {}
  ModelOptions(unsigned threads, unsigned timelimit, bool verbose, std::string wrelaxation) : threads(threads), timelimit(timelimit), verbose(verbose), wrelaxation(wrelaxation) {}
};

void setupSolver(IloCplex solver, const ModelOptions&);

struct ModelStat {
  unsigned rows, cols, nnz;
};

struct SSBF {
  const sualbsp::Instance &I;
  IloEnv env;
  IloModel model;
  IloObjective obj;

  ModelOptions m;
  ModelStat mstat;

  std::map<std::tuple<unsigned, unsigned>, unsigned> ix;
  std::map<std::tuple<unsigned, unsigned, unsigned>, unsigned> ig, ih;

  IloNumVarArray x, u, g, h; 
  IloNumVarArray z, r;       
  unsigned lm, um;

  unsigned nx, ng, nh;

  SSBF(const sualbsp::Instance &I, const ModelOptions &m, unsigned lm, unsigned um) : I(I), model(env), obj(IloAdd(model, IloMinimize(env))), m(m), x(env), u(env), g(env), h(env), z(env), r(env), lm(lm), um(um), nx(0), ng(0), nh(0) {}

  ~SSBF() { env.end(); }

  void addVars();
  void addAssignmentStations();
  void addPrefix();
  void addSequence();
  void addStationsOcurrence();
  void addPrecedences();
  void addCycle();
  void addSingleTask();

  void build();

  void setSolution(IloCplex, const Solution &);
  Solution getSolution(IloCplex);
  ModelStat getStatistics() const;

  std::pair<Solution, IloAlgorithm::Status> solve(const Solution &);
  double solveRelaxed();
};
