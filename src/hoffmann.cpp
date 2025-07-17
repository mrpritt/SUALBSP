#include <cassert>
#include <iostream>
using namespace std;

#include "hoffmann.hpp"
#include "options.hpp"
#include "random.hpp"
#include "solution.hpp"

using namespace sualbsp;

void detail_solution(const Instance *I, vector<Task> hoffmann_sequence, vector<Task> best_solution) {
  unsigned short int currentStation = 0;
  unsigned short int first = 0;
  unsigned short int last;
  unsigned short int task;
  unsigned int acc;
  vector<Task> counter(I->n, 0);

  for (unsigned i = 1; i < I->n; i++) {
  }
  for (unsigned i = 0; i < I->n; i++) {
    task = hoffmann_sequence[i];
    counter[task] = 1;
    if (currentStation != best_solution[task]) {
      last = i - 1;
      acc = I->sb[hoffmann_sequence[last]][hoffmann_sequence[first]];
      for (int j = first; j <= last; j++)
        acc += I->t[hoffmann_sequence[j]];
      for (int j = first; j < last; j++)
        acc += I->sf[hoffmann_sequence[j]][hoffmann_sequence[j + 1]];
      first = i;
      currentStation = best_solution[task];
    }
  }
  last = I->n - 1;
  acc = I->sb[hoffmann_sequence[last]][hoffmann_sequence[first]];
  for (int j = first; j <= last; j++)
    acc += I->t[hoffmann_sequence[j]];
  for (int j = first; j < last; j++)
    acc += I->sf[hoffmann_sequence[j]][hoffmann_sequence[j + 1]];
}

void hoffmann_t::compute_standard_weights(const Instance &I) { set_weights(I.t); }

void hoffmann_t::compute_sw_weights(const Instance &I, double alpha, double beta, double gamma) {
  for (auto i = 0u; i < I.n; i++)
    weights[i] = double(I.t[i]) + alpha * I.get_pw(i) + beta * double(I.F[i].size()) - gamma;
  best = total_weight();
}

void hoffmann_t::compute_s_weights(const Instance &I, double alpha, double beta, double gamma, bool fwd) {
  double df = fwd ? double(I.Fs.size()) / double(I.maxFs) : double(I.Ps.size()) / double(I.maxPs);
  for (auto i = 0u; i < I.n; i++) {
    weights[i] = alpha * double(I.t[i]) / double(I.c) + beta * df - gamma;
    if (weights[i] < 0.0)
      weights[i] = s_epsilon;
  }
  best = total_weight();
}

void hoffmann_t::set_degrees(const std::vector<std::vector<Task>> &adj) {
  const auto n = adj.size();
  for (auto i = 0u; i < n; i++)
    deg[i] = adj[i].size();
}

void hoffmann_t::set_weights(const std::vector<Time> &t) {
  const auto n = t.size();
  for (auto i = 0u; i < n; i++)
    weights[i] = t[i];
  best = total_weight();
}

hoffmann_t::hoffmann_t(const Instance &I) : eligible(I.n, -1), tasks(I.n, -1), best_tasks(I.n, -1), weights(I.n, 0.0), pending(I.n, true) {
  for (unsigned i = 0; i < I.n; i++)
    deg.push_back(I.P[i].size());
}

void hoffmann_t::compute_eligible(unsigned &n_eligible) {
  const auto n = deg.size();
  n_eligible = 0;
  for (auto i = 0u; i != n; ++i)
    if (deg[i] == 0)
      eligible[n_eligible++] = i;
}

void gen_load_f(int depth, int idle, const Instance *I, unsigned n_eligible, hoffmann_t *current, double eval) {
  int i, j, add, n_sub_eligible;

  if (current->n_loads >= opt.maxloads)
    return;

  for (auto ii = 0u; ii < n_eligible; ii++) {
    i = current->eligible[ii];
    if (current->deg[i] == 0 && current->pending[i] == true) {
      if (depth == 0)
        add = I->t[i] + I->sb[i][i];
      else {
        const auto last = current->tasks[depth - 1];
        add = I->t[i] - I->sb[last][current->tasks[0]] + I->sf[last][i] + I->sb[i][current->tasks[0]];
      }
      if (add <= idle) {
        current->tasks[depth] = i;
        eval -= current->weights[i];
        current->deg[i] = -1;
        current->pending[i] = false;
        n_sub_eligible = n_eligible;
        for (auto jj = 0u; jj < I->F[i].size(); jj++) {
          j = I->F[i][jj];
          current->deg[j]--;
          if (current->deg[j] == 0 && current->pending[j] == true)
            current->eligible[n_sub_eligible++] = j;
        }
        if (eval < current->best) {
          for (j = 0; j <= depth; j++)
            current->best_tasks[j] = current->tasks[j];
          current->best = eval;
          current->n_best_tasks = depth + 1;
        }
        current->n_loads++;
        gen_load_f(depth + 1, idle - add, I, n_sub_eligible, current, eval);
        for (auto jj = 0u; jj < I->F[i].size(); jj++) {
          j = I->F[i][jj];
          current->deg[j]++;
        }
        current->pending[i] = true;
        current->deg[i] = 0;
        current->tasks[depth] = 0;
        eval += current->weights[i];
      }
    }
  }
}

void gen_load_b(int depth, int idle, const Instance *I, unsigned n_eligible, hoffmann_t *current, double eval) {
  int i, j, add, n_sub_eligible;
  if (current->n_loads >= opt.maxloads)
    return;
  for (auto ii = 0u; ii < n_eligible; ii++) {
    i = current->eligible[ii];
    if (current->deg[i] == 0 && current->pending[i] == true) {
      if (depth == 0)
        add = (I->t[i] + I->sb[i][i]);
      else {
        const auto last = current->tasks[depth - 1];
        add = I->t[i] - I->sb[current->tasks[0]][last] + I->sf[i][last] + I->sb[current->tasks[0]][i];
      }
      if (add <= idle) {
        current->tasks[depth] = i;
        eval -= current->weights[i];
        current->deg[i] = -1;
        current->pending[i] = false;
        n_sub_eligible = n_eligible;
        for (auto jj = 0u; jj < I->P[i].size(); jj++) {
          j = I->P[i][jj];
          current->deg[j]--;
          if (current->deg[j] == 0 && current->pending[j] == true)
            current->eligible[n_sub_eligible++] = j;
        }
        if (eval < current->best) {
          for (j = 0; j <= depth; j++)
            current->best_tasks[j] = current->tasks[j];
          current->best = eval;
          current->n_best_tasks = depth + 1;
        }
        current->n_loads++;
        gen_load_b(depth + 1, idle - add, I, n_sub_eligible, current, eval);
        for (auto jj = 0u; jj < I->P[i].size(); jj++) {
          j = I->P[i][jj];
          current->deg[j]++;
        }
        current->pending[i] = true;
        current->deg[i] = 0;
        current->tasks[depth] = 0;
        eval += current->weights[i];
      }
    }
  }
}

unsigned short int partial_hoffmann_f(const Instance *I, hoffmann_t *c, vector<Task> &partial_assignment, vector<Task> &partial_solution, unsigned n_assigned) {
  unsigned short int n_stations = 0;
  unsigned short int i, ii, j;
  vector<bool> pending(I->n, 0);
  vector<short int> degrees(I->n, 0);
  double eval = 0;
  int assigned = 0;

  for (i = 0; i < I->n; i++)
    if (c->pending[i] == true)
      eval += c->weights[i];

  c->best = eval;
  unsigned n_eligible = 0;

  pending = c->pending;
  degrees = c->deg;
  for (i = 0; i < I->n; i++)
    if (c->deg[i] == 0 && c->pending[i] == true) {
      c->eligible[n_eligible] = i;
      n_eligible++;
    }

  while (n_assigned < I->n) {
    c->n_best_tasks = 0;
    c->n_loads = 0;
    gen_load_f(0, I->c, I, n_eligible, c, eval);
    for (ii = 0; ii < c->n_best_tasks; ii++) {
      i = c->best_tasks[ii];
      partial_solution[i] = n_stations;
      partial_assignment[assigned] = i;
      n_assigned++;
      assigned++;
      degrees[i] = -1;
      pending[i] = false;
      for (auto jj = 0u; jj < I->F[i].size(); jj++) {
        j = I->F[i][jj];
        degrees[j]--;
      }
    }
    if (n_assigned < I->n) {
      c->deg = degrees;
      c->pending = pending;
      c->compute_eligible(n_eligible);
      eval = c->best;
    }
    n_stations++;
  }
  return n_stations;
}

unsigned short int partial_hoffmann_b(const Instance *I, hoffmann_t *c, vector<Task> &partial_assignment, vector<Task> &partial_solution, unsigned n_assigned) {
  unsigned short int n_stations = 0;
  unsigned short int i, j;
  vector<bool> pending(I->n, 0);
  vector<short int> degrees(I->n, 0);
  double eval = 0;

  for (auto i = 0u; i < I->n; i++)
    if (c->pending[i] == true)
      eval += c->weights[i];

  c->best = eval;
  unsigned n_eligible = 0;

  pending = c->pending;
  degrees = c->deg;
  for (auto i = 0u; i < I->n; i++)
    if (c->deg[i] == 0 && c->pending[i] == true) {
      c->eligible[n_eligible] = i;
      n_eligible++;
    }

  while (n_assigned < I->n) {
    c->n_best_tasks = 0;
    c->n_loads = 0;
    gen_load_b(0, I->c, I, n_eligible, c, eval);
    for (short ii = 0; ii < c->n_best_tasks; ii++) {
      i = c->best_tasks[ii];
      partial_solution[i] = n_stations;
      partial_assignment[n_assigned] = i;
      n_assigned++;
      degrees[i] = -1;
      pending[i] = false;
      for (auto jj = 0u; jj < I->P[i].size(); jj++) {
        j = I->P[i][jj];
        degrees[j]--;
      }
    }
    if (n_assigned < I->n) {
      c->deg = degrees;
      c->pending = pending;
      c->compute_eligible(n_eligible);
      eval = c->best;
    }
    n_stations++;
  }
  return n_stations;
}

unsigned short int hoffmann_f(const Instance *I, hoffmann_t *f, vector<unsigned short int> &best_assignment, vector<unsigned short int> &best_solution) {
  unsigned short int n_stations = 0;
  unsigned short int i, j;
  unsigned short int n_assigned = 0;
  vector<short int> degrees;
  vector<bool> pending(I->n, true);
  double eval = f->best;
  unsigned n_eligible = 0;

  f->set_degrees(I->P);
  f->pending = pending;
  f->compute_eligible(n_eligible);
  degrees = f->deg;

  while (n_assigned < I->n) {
    f->n_best_tasks = 0;
    f->n_loads = 0;
    gen_load_f(0, I->c, I, n_eligible, f, eval);
    for (short ii = 0; ii < f->n_best_tasks; ii++) {
      i = f->best_tasks[ii];
      best_solution[i] = n_stations;
      best_assignment[n_assigned] = i;
      n_assigned++;
      degrees[i] = -1;
      pending[i] = false;
      for (auto jj = 0u; jj < I->F[i].size(); jj++) {
        j = I->F[i][jj];
        degrees[j]--;
      }
    }
    if (n_assigned < I->n) {
      f->deg = degrees;
      f->pending = pending;
      f->compute_eligible(n_eligible);
      eval = f->best;
    }
    n_stations++;
  }
  return n_stations;
}

unsigned short int hoffmann_b(const Instance *I, hoffmann_t *f, vector<unsigned short int> &best_assignment, vector<unsigned short int> &best_solution) {
  unsigned short int n_stations = 0;
  unsigned short int i, j;
  unsigned short int n_assigned = 0;
  vector<short int> degrees;
  vector<bool> pending(I->n, true);
  double eval = f->best;
  unsigned n_eligible = 0;

  f->set_degrees(I->F);
  f->pending = pending;
  f->compute_eligible(n_eligible);
  degrees = f->deg;

  while (n_assigned < I->n) {
    f->n_best_tasks = 0;
    f->n_loads = 0;
    gen_load_b(0, I->c, I, n_eligible, f, eval);
    for (short ii = 0u; ii < f->n_best_tasks; ii++) {
      i = f->best_tasks[ii];
      best_solution[i] = n_stations;
      best_assignment[n_assigned] = i;
      n_assigned++;
      degrees[i] = -1;
      pending[i] = false;
      for (auto jj = 0u; jj < I->P[i].size(); jj++) {
        j = I->P[i][jj];
        degrees[j]--;
      }
    }
    if (n_assigned < I->n) {
      f->deg = degrees;
      f->pending = pending;
      f->compute_eligible(n_eligible);
      eval = f->best;
    }
    n_stations++;
  }
  vector<Task> hoffmann_sequence(best_assignment);
  reverse(hoffmann_sequence.begin(), hoffmann_sequence.end());
  best_assignment = hoffmann_sequence;
  for (i = 0; i < I->n; i++)
    best_solution[i] = n_stations - 1 - best_solution[i];
  detail_solution(I, hoffmann_sequence, best_solution);
  return n_stations;
}

unsigned short int prepare_partial_forward(const Instance *I, hoffmann_t *c, vector<Task> order, vector<Task> station, unsigned short int *n_stations) {
  unsigned short int repeatUpToStation = 0;
  int j;

  c->set_degrees(I->P);
  for (auto i = 0u; i < I->n; i++)
    if (c->pending[order[i]] == false) {
      repeatUpToStation = station[order[i]];
      break;
    }

  *n_stations = repeatUpToStation;
  for (int t = 0; station[order[t]] < repeatUpToStation; t++) {
    auto i = order[t];
    c->deg[i] = -1;
    c->pending[i] = false;
    for (auto jj = 0u; jj < I->F[i].size(); jj++) {
      j = I->F[i][jj];
      c->deg[j]--;
    }
  }

  return c->num_assigned();
}

unsigned short int prepare_partial_backward(const Instance *I, hoffmann_t *c, vector<Task> order, vector<Task> station, unsigned short int *n_stations) {
  unsigned short int repeatUpToStation = 0;
  int i, j, small = 0, big = 0;

  c->set_degrees(I->F);

  for (i = I->n - 1; i > 0; i--) {
    if (station[order[i]] > big)
      big = station[order[i]];
    if (c->pending[order[i]] == false) {
      repeatUpToStation = station[order[i]];
      small = station[order[i]];
      break;
    }
  }
  *n_stations += (big - small);
  for (int t = I->n - 1; station[order[t]] != repeatUpToStation; t--) {
    i = order[t];
    c->deg[i] = -1;
    c->pending[i] = false;
    for (auto jj = 0u; jj < I->P[i].size(); jj++) {
      j = I->P[i][jj];
      c->deg[j]--;
    }
  }

  return c->num_assigned();
}

unsigned short int hoffmann_bidirectional(const Instance *I, hoffmann_t *f, hoffmann_t *b, vector<Task> &best_assignment, vector<Task> &best_solution) {
  int workload;
  unsigned short int n_best;
  unsigned short int n_stations, n_stations2;
  unsigned short int n_assigned;
  hoffmann_t c(*I);

  vector<unsigned short int> forward_solution(I->n, -1);
  vector<unsigned short int> forward_assignment(I->n, -1);
  unsigned short int n_forward = hoffmann_f(I, f, forward_assignment, forward_solution);

  vector<unsigned short int> backward_solution(I->n, -1);
  vector<unsigned short int> backward_assignment(I->n, -1);
  unsigned short int n_backward = hoffmann_b(I, b, backward_assignment, backward_solution);
  n_best = min(n_forward, n_backward);
  if (n_best == n_forward) {
    best_assignment = forward_assignment;
    best_solution = forward_solution;
  } else {
    best_assignment = backward_assignment;
    best_solution = backward_solution;
  }
  for (auto i = 0u; i < I->n; i++)
    for (auto i = 0u; i < I->n; i++)

      for (int station = n_forward - 2; station > 0; station--) {
        fill(c.pending.begin(), c.pending.end(), false);
        workload = 0;
        for (auto i = 0u; i < I->n; i++)
          if (forward_solution[i] >= station) {
            c.pending[i] = true;
            workload += I->t[i];
          }

        {
          n_stations = station;
          n_stations2 = 0;
          n_assigned = prepare_partial_backward(I, &c, backward_assignment, backward_solution, &n_stations2);
          workload = 0;
          for (auto i = 0u; i < I->n; i++)
            if (c.pending[i] == true)
              workload += I->t[i];

          if ((n_stations + n_stations2 + workload / I->c + (workload % I->c != 0)) < n_best) {
            vector<Task> partial_solution(I->n, numeric_limits<Task>::max());
            vector<Task> partial_assignment(I->n, numeric_limits<Task>::max());
            c.weights = b->weights;

            unsigned short int extra_stations = partial_hoffmann_b(I, &c, partial_assignment, partial_solution, n_assigned);
            if (extra_stations + n_stations + n_stations2 < n_best) {
              n_best = extra_stations + n_stations + n_stations2;
              for (auto i = 0u; i < I->n; i++)
                best_solution[i] = numeric_limits<Task>::max();
              for (auto i = 0u; i < I->n; i++)
                best_assignment[i] = numeric_limits<Task>::max();
              int ass = 0;
              for (auto i = 0u; i < I->n; i++) {
                if (forward_solution[forward_assignment[i]] < station) {
                  best_assignment[ass] = forward_assignment[ass];
                  best_solution[forward_assignment[ass]] = forward_solution[forward_assignment[ass]];
                  ass++;
                } else
                  break;
              }
              for (int i = I->n - 1; i >= 0; i--) {
                if (partial_assignment[i] == numeric_limits<Task>::max())
                  break;
                best_assignment[ass] = partial_assignment[i];
                best_solution[partial_assignment[i]] = n_stations + (extra_stations - partial_solution[partial_assignment[i]] - 1);
                ass++;
              }
              for (auto i = 0u; i < I->n; i++) {
                if (best_solution[backward_assignment[i]] == numeric_limits<Task>::max()) {
                  best_assignment[ass] = backward_assignment[i];
                  best_solution[best_assignment[ass]] = extra_stations + n_stations + n_stations2 - 1 - (n_backward - 1 - backward_solution[backward_assignment[i]]);
                  ass++;
                }
              }
              detail_solution(I, best_assignment, best_solution);
            }
          }
        }
      }

  for (int station = n_backward - 2; station > 1; station--) {
    workload = 0;
    for (auto i = 0u; i < I->n; i++) {
      workload += I->t[i];
      c.pending[i] = true;
    }
    for (auto i = 0u; i < I->n; i++) {
      if (backward_solution[i] >= station) {
        c.pending[i] = false;
        workload -= I->t[i];
      }
    }
    {
      n_stations = n_backward - station;
      n_stations2 = 0;
      n_assigned = prepare_partial_forward(I, &c, forward_assignment, forward_solution, &n_stations2);
      workload = 0;
      for (auto i = 0u; i < I->n; i++)
        if (c.pending[i] == true)
          workload += I->t[i];

      if ((n_stations + n_stations2 + workload / I->c + (workload % I->c != 0)) < n_best) {
        vector<Task> partial_solution(I->n, numeric_limits<Task>::max());
        vector<Task> partial_assignment(I->n, numeric_limits<Task>::max());
        c.weights = b->weights;

        unsigned short int extra_stations = partial_hoffmann_f(I, &c, partial_assignment, partial_solution, n_assigned);
        if (extra_stations + n_stations + n_stations2 < n_best) {
          n_best = extra_stations + n_stations + n_stations2;
          for (auto i = 0u; i < I->n; i++)
            best_solution[i] = numeric_limits<Task>::max();
          for (auto i = 0u; i < I->n; i++)
            best_assignment[i] = numeric_limits<Task>::max();
          unsigned ass = 0;
          for (auto i = 0u; i < I->n; i++) {
            if (forward_solution[forward_assignment[i]] < n_stations2) {
              best_assignment[ass] = forward_assignment[ass];
              best_solution[forward_assignment[ass]] = forward_solution[forward_assignment[ass]];
              ass++;
            } else
              break;
          }
          for (auto i = 0u; i < I->n; i++) {
            if (partial_assignment[i] == numeric_limits<Task>::max())
              break;
            best_assignment[ass] = partial_assignment[i];
            best_solution[partial_assignment[i]] = n_stations2 + partial_solution[partial_assignment[i]];
            ass++;
          }
          for (auto i = ass; i < I->n; i++) {
            best_assignment[ass] = backward_assignment[i];
            best_solution[backward_assignment[i]] = backward_solution[backward_assignment[i]] - n_backward + n_best;
            ass++;
          }
          detail_solution(I, best_assignment, best_solution);
        }
      }
    }
  }

  return n_best;
}

unsigned short int Hoffmann(Instance &I) {
  vector<unsigned short int> best_solution(I.n, -1);
  vector<unsigned short int> best_assignment(I.n, -1);

  hoffmann_t f(I);

  f.compute_standard_weights(I);

  unsigned short int n_stations = hoffmann_f(&I, &f, best_assignment, best_solution);

  if (incumbent.m > n_stations)
    incumbent.set(I, best_assignment, best_solution, n_stations);
  return n_stations;
}

unsigned short int SW_Hoffmann(Instance &I, unsigned short lb) {
  unsigned short ms = numeric_limits<unsigned short>::max();
  vector<unsigned short int> best_solution(I.n, -1);
  vector<unsigned short int> best_assignment(I.n, -1);

  hoffmann_t f(I);

  for (double alpha = 0.000; alpha <= 0.02; alpha += 0.005) {
    for (double beta = 0.000; beta <= 0.02; beta += 0.005) {
      for (double gamma = 0; gamma <= 0.03; gamma += 0.01) {
        S.iter++;
        f.compute_sw_weights(I, alpha, beta, gamma);

        unsigned short int n_stations = hoffmann_f(&I, &f, best_assignment, best_solution);

        if (n_stations < ms)
          ms = n_stations;
        incumbent.set(I, best_assignment, best_solution, ms);
      }
    }
  }
  return ms;
}

unsigned short int FH_Hoffmann(Instance &I, unsigned short lb) {
  vector<Task> best_solution(I.n, 0);
  vector<Task> best_assignment(I.n, 0);

  hoffmann_t f(I), b(I);

  f.compute_standard_weights(I);
  b.compute_standard_weights(I);

  unsigned short int n_stations = hoffmann_bidirectional(&I, &f, &b, best_assignment, best_solution);

  detail_solution(&I, best_assignment, best_solution);
  incumbent.set(I, best_assignment, best_solution, n_stations);
  return n_stations;
}

unsigned short int S_Hoffmann(Instance &I, unsigned short lb) {
  double alpha, beta, gamma;
  unsigned short ms = numeric_limits<unsigned short>::max();
  unsigned short m;
  vector<unsigned short int> best_solution(I.n, -1);
  vector<unsigned short int> best_assignment(I.n, -1);
  vector<unsigned short int> current_best_solution(I.n, -1);
  vector<unsigned short int> current_best_assignment(I.n, -1);

  hoffmann_t f(I), b(I);

  while (logging::elapsed() < opt.tlim && S.iter++ < opt.ilim && ms > lb) {
    alpha = getRandom(0.0, 1.0);
    beta = getRandom(0.0, 1.0);
    gamma = getRandom(0.0, 1.0 / (double)(I.n));

    f.compute_s_weights(I, alpha, beta, gamma, true);
    b.compute_s_weights(I, alpha, beta, gamma, false);

    m = hoffmann_bidirectional(&I, &f, &b, best_assignment, best_solution);

    if (m < ms)
      ms = m;
    incumbent.set(I, best_assignment, best_solution, ms);
  }
  return ms;
}
