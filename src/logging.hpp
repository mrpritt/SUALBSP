#pragma once

#include <iostream>
#include <vector>

#include <boost/any.hpp>

#include "timer.hpp"

#define FMT_HEADER_ONLY
#include "fmt/color.h"
#include "fmt/format.h"
#include "fmt/ostream.h"
#include "fmt/ranges.h"

namespace logging {
extern timer start;
inline double elapsed() { return start.elapsed(); }
} // namespace logging

std::string get_date_string(Timepoint t);
