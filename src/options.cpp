#include "options.hpp"

#include <cassert>
#include <cstdlib>
#include <iostream>
using namespace std;

#include <boost/program_options/value_semantic.hpp>

#include <sys/ioctl.h>

#include "dimensions.hpp"
#include "logging.hpp"
#include "random.hpp"
#include "util.hpp"
#include "version.hpp"
using namespace sualbsp;

Options opt;

void add_std_options(po::options_description &desc, standardOptions &opt) { desc.add_options()("help", "Show help.")("version", "Show version.")("instance", po::value<string>(&opt.instance), "Instance.")("format", po::value<string>(&opt.format)->default_value("guess"), "Input format, either Martino&Pastor (mp), Scholl,Boysen,Fliedner (sbf), or guess based on the extension.")("seed", po::value<unsigned>(&opt.seed)->default_value(1), "Random seed (0 for random value)."); }

unsigned get_terminal_width() {
  struct winsize wsize;
  ioctl(STDOUT_FILENO, TIOCGWINSZ, &wsize);
  return wsize.ws_col;
}

std_description::std_description(const string &caption, standardOptions &opt) : Base(caption, get_terminal_width()) { add_std_options(*this, opt); }

bool process_options(int argc, char *argv[], Options &opt, po::variables_map &vm) {

  std_description all("General options", opt);
  all.add_options()("type2", po::bool_switch(&opt.type2)->default_value(false), "Type 2 instance (used only for pre-processing).");

  po::options_description algdesc("Algorithmic options", get_terminal_width());
  algdesc.add_options()("algorithm", po::value<string>(&opt.algorithm)->default_value("sampler"), "Algorithm (sampler,rule,ruleopt,hoffmann,SW_hoffmann,FH_hoffmann,S_hoffmann,none).")("rule,R", po::value<string>(&opt.rule)->default_value("norule"), "Rule for rule-based assignment (maxtminsu,maxtminsuslack,maxf,minslack,random,mprule,sbf,randomfs,asq); implies algorithm \"rule\".")("mprule", po::bool_switch(&opt.mprule_once)->default_value(false), "Run M&P's rule once as initial heuristic (in both directions).")("maxloads", po::value<unsigned>(&opt.maxloads)->default_value(1000), "Maximum number of loads per station in Hoffmann.")("iterlimit,i", po::value<unsigned>(&opt.ilim)->default_value(uint_inf), "Maximum number of iterations.")("riterlimit,r", po::value<unsigned>(&opt.rilim)->default_value(0)->implicit_value(uint_inf), "Maximum number of iterations on reverse instance (no option: none, only option: same number as forward).")("timelimit,t", po::value<unsigned>(&opt.tlim)->default_value(600), "Time limit in seconds.")("onlylb", po::bool_switch(&opt.onlylb)->default_value(false), "Show only lower bounds.")("heavy-preprocess", po::bool_switch(&opt.heavypre)->default_value(false), "Enable heavy pre-processing.");

  po::options_description mod("Solver options", get_terminal_width());
  mod.add_options()("dwlb", po::bool_switch(&opt.dw)->default_value(false), "Evaluate Dantzig-Wolfe reformulation lower bound.")("dwidlebound", po::bool_switch(&opt.dwidlebound)->default_value(false), "Bound idle time during pricing branch and bound.")("pricing,p", po::value<string>(&opt.dwtype)->default_value("ap"), "Type of pricing for the Dantzig-Wolfe reformulation (salbp,ap).")("maxtimemodel", po::value<unsigned>(&opt.maxtimemodel)->default_value(86400), "Maximum time for DW and SSBF (default limit 24 hours).")("ssbf", po::bool_switch(&opt.ssbf)->default_value(false), "Evaluate SSBF lower bound.")("ssbf-ub-opt", po::bool_switch(&opt.ssbf_ub_opt)->default_value(false), "Use the optimal value when evaluating the SSBF lower bound.")("test-lb", po::bool_switch(&opt.test_lb)->default_value(false), "For testing lower bounds: models are not built using best simple lower bounds.")("maxtasks", po::value<unsigned>(&opt.maxtasks)->default_value(100), "Upper bound for number of tasks to apply SSBF.");

  po::options_description outdesc("Output options", get_terminal_width());
  outdesc.add_options()("tag", po::value<string>(&opt.tag)->default_value(""), "Optional tag which will be added to algorithm description.")("oinstance", po::value<string>(&opt.oinstance)->default_value(""), "Write pre-processed instance to file.")("output", po::value<string>(&opt.output)->default_value("/dev/null"), "Output file name for solution.")("wrelaxation", po::value<string>(&opt.wrelaxation)->default_value("/dev/null"), "Output file name for linear relaxation of SSBF model.")("trace-models", po::bool_switch(&opt.trace_models)->default_value(false), "Trace models output.");

  po::options_description bbrdesc("BBR options", get_terminal_width());
  bbrdesc.add_options()("cbfs", po::bool_switch(&opt.cbfs)->default_value(false), "Use cbfs (false by default).")("secondpass", po::bool_switch(&opt.secondpass)->default_value(false), "Second pass (no maxnodescbfs in second step).")("maxnodescbfs", po::value(&opt.maxnodescbfs)->default_value(1500), "Maximum extensions.")("exacttimelimit", po::value<unsigned>(&opt.tlimexact)->default_value(600), "Time limit in seconds for the exact method.")("maxmemory", po::value<unsigned>(&opt.maxmemory)->default_value(20000000), "Max states in memory.")("uselm1p", po::bool_switch(&opt.uselm1p)->default_value(false), "Use lm1 with setup times (false by default).")("enhancedlookup", po::bool_switch(&opt.enhancedlookup)->default_value(false), "Enhanced lookup (time intensive).")("usecounting", po::bool_switch(&opt.usecounting)->default_value(false), "Use counting bounds in cbfs (false by default).")("uselm4", po::bool_switch(&opt.uselm4)->default_value(false), "Use lm4 bound in cbfs (false by default).");

  all.add(algdesc).add(mod).add(outdesc).add(bbrdesc);

  po::positional_options_description pod;
  pod.add("instance", 1);

  po::store(po::command_line_parser(argc, argv).options(all).positional(pod).run(), vm);
  po::notify(vm);

  if (vm.count("version")) {
    fmt::print("{}\n", version);
    return 0;
  }

  if (vm.count("help") || !vm.count("instance")) {
    cout << all << endl;
    return false;
  }

  opt.seed = setupRandom(opt.seed);

  if (opt.rule != "norule" && opt.algorithm == "sampler")
    opt.algorithm = "rule";

  if (opt.rilim > 0 && opt.rilim == uint_inf)
    opt.rilim = opt.ilim;

  return true;
}
