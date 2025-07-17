#pragma once

#include <boost/program_options.hpp>
namespace po = boost::program_options;

struct standardOptions {
  std::string instance;
  std::string format;
  unsigned seed;
};
struct Options : public standardOptions {
  std::string algorithm;
  std::string tag;
  std::string oinstance;

  std::string output;

  unsigned tlim;
  unsigned ilim;
  unsigned rilim;
  std::string rule;
  unsigned maxloads;
  unsigned maxtasks;

  bool onlylb;
  bool dwidlebound;
  bool dw, ssbf, ssbf_ub_opt;
  std::string dwtype;
  unsigned maxtimemodel;
  bool mprule_once;
  bool type2;

  bool cbfs;
  bool secondpass;
  unsigned maxnodescbfs;
  unsigned tlimexact;
  unsigned maxmemory;
  bool uselm1p;
  bool enhancedlookup;
  bool usecounting;
  bool uselm4;
  bool heavypre;
  bool trace_models;
  bool test_lb;

  std::string wrelaxation;
};

extern Options opt;

void add_std_options(boost::program_options::options_description &, standardOptions &);
unsigned get_terminal_width();

struct std_description : public boost::program_options::options_description {
  typedef boost::program_options::options_description Base;
  std_description(const std::string &caption, standardOptions &);
};

bool process_options(int, char *[], Options &, po::variables_map &);
