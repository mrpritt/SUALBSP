#include "timer.hpp"

#include <stdexcept>

using namespace std;

string get_date_string(Timepoint t) {
  auto as_time_t = chrono::system_clock::to_time_t(t);
  struct tm tm;
  char some_buffer[128];
  if (::gmtime_r(&as_time_t, &tm))
    if (strftime(some_buffer, sizeof(some_buffer), "%F %T", &tm))
      return string{some_buffer};
  throw runtime_error("Failed to get current date as string");
}
