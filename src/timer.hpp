#pragma once

#include <chrono>
#include <string>

typedef std::chrono::system_clock Clock;
typedef std::chrono::duration<double> Duration;
typedef Clock::time_point Timepoint;

struct timer {
  using clock = std::chrono::steady_clock;
  using time_point = clock::time_point;

  time_point tp_start;
  double tm_lim;

  timer(double timeLimitSecs, const timer &parent) { reset(std::min(timeLimitSecs, parent.left())); }

  timer(double time_lim_secs = std::numeric_limits<double>::max()) { reset(time_lim_secs); }

  void reset(double time_lim_secs = std::numeric_limits<double>::max()) { tm_lim = time_lim_secs, tp_start = clock::now(); }

  double elapsed() const { return std::chrono::duration_cast<std::chrono::duration<double>>(clock::now() - tp_start).count(); }

  double left() const { return tm_lim - elapsed(); }
  bool timed_out() const { return elapsed() >= tm_lim; }
};

std::string get_date_string(Timepoint t);
