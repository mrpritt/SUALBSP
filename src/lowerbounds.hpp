#pragma once

#include "instance.hpp"

struct tail {
  Task id;
  unsigned t;
};

namespace sualbsp {

extern vector<tail> tails; 
extern vector<tail> heads; 

inline unsigned round(Time t, Time c) { return (t + c - 1) / c; }

inline unsigned lm1Unrounded(const vector<Time> &t, Time c) { return accumulate(t.begin(), t.end(), 0, plus<Time>()); }

void new_incumbent(Instance &I);

inline unsigned lm1(const vector<Time> &tm, Time c) { return round(lm1Unrounded(tm, c), c); }

unsigned lms1_fb(const Instance &I);
unsigned lms1_fb_subset(const Instance &I, const vector<bool> &active, bool istails);

unsigned lm2(const Instance &I);
unsigned lm3(const Instance &I);

unsigned lm4(const Instance &I);

unsigned lms1_f(const Instance &I);
unsigned lms1_f_subset(const Instance &I, const vector<bool> &active, bool istails);

inline unsigned lms1(const Instance &I) {
  if (I.setupDirected)
    return lms1_fb(I);
  else
    return lms1_f(I);
}

Time lc_tasks(const Instance &, const vector<Task> &);
} 
