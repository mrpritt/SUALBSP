#include <algorithm>
#include <cassert>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <queue>
#include <vector>
namespace fs = std::filesystem;
using namespace std;

#include "bbr.hpp"
#include "dw.hpp"
#include "hoffmann.hpp"
#include "instance.hpp"
#include "logging.hpp"
#include "lowerbounds.hpp"
#include "models.hpp"
#include "options.hpp"
#include "prbp.hpp"
#include "random.hpp"
#include "solution.hpp"
#include "util.hpp"
using namespace sualbsp;

string iname;

bool run(Instance &I, const Options &opt, unsigned short lb) {
  unsigned short int ms;
  if (opt.algorithm == "sampler")
    ms = sample(I, lb);
  else if (opt.algorithm == "hoffmann")
    ms = Hoffmann(I);
  else if (opt.algorithm == "SW_hoffmann")
    ms = SW_Hoffmann(I, lb);
  else if (opt.algorithm == "S_hoffmann")
    ms = S_Hoffmann(I, lb);
  else if (opt.algorithm == "FH_hoffmann")
    ms = FH_Hoffmann(I, lb);
  else if (opt.algorithm == "rule")
    ms = assignRule(I, lb, opt.rule);
  else if (opt.algorithm == "ruleopt")
    ms = assignRuleOpt(I, lb, opt.rule);
  else if (opt.algorithm == "none")
    ms = incumbent.m;
  else {
    return false;
  }

  if (ms < incumbent.m)
    incumbent.m = ms;

  if (incumbent.m < lb) {
    return false;
  }

  return true;
}

void create_reverse_instance(const Instance &I, Instance &Ir) {
  Ir = I;
  Ir.reverse();
}

void rule_fwd_bwd(Instance &I, Instance &Ir) {
  if (opt.rilim > 0 && (opt.algorithm == "rule" || opt.algorithm == "ruleopt")) {

    if (!opt.mprule_once)
      create_reverse_instance(I, Ir);

    unsigned short int ms;
    opt.ilim += opt.rilim;

    if (opt.algorithm == "rule")
      ms = assignRuleAlt(I, Ir, incumbent.lb, opt.rule);
    else {
      ms = assignRuleOptAlt(I, Ir, incumbent.lb, opt.rule);
    }
    if (ms < incumbent.m)
      incumbent.m = ms;

    if (incumbent.m < incumbent.lb) {
      return throw "LB/UB mismatch.";
    }
  } else {
    bool success;
    unsigned bwdtime = opt.tlim / 2;

    if (opt.rilim > 0)
      opt.tlim -= bwdtime;
    success = run(I, opt, incumbent.lb);

    if (success && opt.rilim > 0 && incumbent.m > incumbent.lb) {
      if (!opt.mprule_once)
        create_reverse_instance(I, Ir);
      opt.ilim += opt.rilim; 
      opt.tlim += bwdtime;
      success = run(Ir, opt, incumbent.lb);
    }
  }
}

int main(int argc, char *argv[]) {
  S.iter = 0;

  po::variables_map vm;
  if (!process_options(argc, argv, opt, vm))
    return 1;

  iname = string(fs::path(opt.instance).parent_path().filename()) + "/" + string(fs::path(opt.instance).stem());
  Instance I, Ir;
  ifstream in(opt.instance);
  if (!in.is_open()) {
    fmt::print(cerr, "Cannot open {}.\n", iname);
    exit(1);
  }
  if (opt.format == "guess")
    opt.format = guess_format(opt.instance);
  I.read(in, opt.format);
  in.close();

  string alg = opt.algorithm;
  if (alg == "rule" || alg == "ruleopt")
    alg += "/" + opt.rule;
  if (opt.tag != "")
    alg += "/" + opt.tag;
  if (opt.cbfs)
    alg += "/cbfs";

  try {
    unsigned short lm1v = lm1(I.t, I.c);
    unsigned short lms1v = lms1(I);
    unsigned short lm2v = lm2(I);
    unsigned short lm3v = lm3(I);
    unsigned short lm4v = lm4(I);

    incumbent.m = I.n;
    incumbent.pi.resize(I.n);
    incumbent.a.resize(I.n);
    iota(incumbent.pi.begin(), incumbent.pi.end(), 0);
    incumbent.lb = max({lm1v, lms1v, lm2v, lm3v, lm4v});

    if (opt.mprule_once) {
      assignRuleOpt(I, incumbent.lb, "mprule");

      if (opt.rilim > 0 && incumbent.m > incumbent.lb) {
        create_reverse_instance(I, Ir);
        auto ci = incumbent.m;
        assignRuleOpt(Ir, incumbent.lb, "mprule");
        if (ci > incumbent.m)
          incumbent.reverse(I);
      }
    }

    ModelOptions m{1, Time(ceil(logging::elapsed())) + opt.maxtimemodel, opt.trace_models, opt.wrelaxation}; 
    unsigned short model_lb = opt.test_lb ? max({lm1v, lms1v}) : incumbent.lb;

    unsigned short lms1dw = 0;
    double timedw = 0.0;

    if (opt.dw) {
      timer t;
      lms1dw = dw(I, m);
      incumbent.lb = max(model_lb, lms1dw);
      timedw = t.elapsed();
    }

    unsigned short lms1ssbf = 0;
    double timessbf = 0.0;
    if (opt.ssbf && I.n <= opt.maxtasks) {
      timer t;
      SSBF ssbf(I, m, model_lb, opt.ssbf_ub_opt ? I.optm : incumbent.m);
      ssbf.build();
      try {
        lms1ssbf = ceil(ssbf.solveRelaxed());
      } catch (char const *msg) {
        lms1ssbf = 0;
        fmt::print(fg(fmt::color::red), "Solving SSBF relaxation of {} failed.\n", iname);
      }

      incumbent.lb = max(incumbent.lb, lms1ssbf);
      timessbf = t.elapsed();
    }

    fmt::print("LOWERBOUNDS {} {} {} {} {} {} {} {} {} {} {} {} {}\n", get_date_string(Clock::now()), logging::elapsed(), iname, I.n, lm1v, lms1v, lm2v, lm3v, lm4v, lms1dw, timedw, lms1ssbf, timessbf);
    if (opt.onlylb)
      return 0;

    opt.tlim += logging::elapsed();
    if (incumbent.m > incumbent.lb && !opt.onlylb) {
      rule_fwd_bwd(I, Ir);
    } else if (!opt.onlylb) {
    }

    if (incumbent.m > incumbent.lb && opt.cbfs) {
      bool isOptimal;
      isOptimal = bbr::bbr(I, opt, incumbent.lb);
      if (isOptimal) {
        incumbent.lb = incumbent.m;
      } else {
        if (opt.secondpass && !isOptimal && incumbent.m > incumbent.lb) {
          if (logging::elapsed() < opt.tlimexact) {
            opt.tlimexact -= unsigned(logging::elapsed());
            opt.maxnodescbfs = opt.maxmemory;
            isOptimal = bbr::bbr(I, opt, incumbent.lb);
            if (isOptimal)
              incumbent.lb = incumbent.m;
          }
        }
      }
    }

    fmt::print("SUMMARY {} {} {} {} {} {} {} {} {} {} {} {}\n", get_date_string(Clock::now()), logging::elapsed(), opt.tlim, opt.ilim, alg, iname, opt.seed, I.n, incumbent.lb, incumbent.m, incumbent.time, incumbent.iter);
  } catch (char const *_msg) {
    auto msg = fmt::format("ERROR {} {} {} {} {}", _msg, alg, iname, opt.seed, opt.rule);
    return 1;
  }
}
