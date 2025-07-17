#pragma once

#include "instance.hpp"
#include "options.hpp"
#include "solution.hpp"
#include <map>
#include <queue>

using namespace std;
using namespace sualbsp;

namespace bbr {

typedef struct setupStruct {
  Task t1;
  Task t2;
  unsigned t;
} setupStruct;

bool bbr(const Instance &I, Options &opt, unsigned short lb);

} // namespace bbr
