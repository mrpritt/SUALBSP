#include <cassert>

#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

#include "dw.hpp"
#include "logging.hpp"
#include "options.hpp"
#include "sop.hpp"

#define RC_EPS 0.00001

ILOLAZYCONSTRAINTCALLBACK6(myLazyCallback, IloBoolVarArray, x, IloNumArray, values, const Instance *, I, dw_t *, d, int, nVertices, int, initial) {
  bool found;
  int subtourLength, numVisited = 0;
  int currentCity = initial;
  int initialCity = initial;
  vector<bool> visited(nVertices, false);
  vector<bool> subtour(nVertices, false);
  vector<int> tour;

  getValues(values, x);
  subtour[currentCity] = true;
  visited[currentCity] = true;
  subtourLength = 1;
  numVisited++;
  tour.push_back(currentCity);

  for (;;) {
    for (;;) {
      found = false;
      for (int j = 0; j < nVertices; j++) {
        if (values[currentCity * nVertices + j] > (1 - RC_EPS)) {
          if (j == initialCity)
            break;
          currentCity = j;
          tour.push_back(currentCity);
          subtour[j] = true;
          visited[j] = true;
          found = true;
          subtourLength++;
          numVisited++;
          break;
        }
      }
      if (found == false) {
        if (subtourLength != nVertices) {
          IloExpr expr1(getEnv());
          IloExpr expr2(getEnv());
          for (int i = 0; i < nVertices; i++) {
            for (int j = 0; j < nVertices; j++) {
              if ((i != j) && (subtour[i] == true) && (subtour[j] == false)) {
                expr1 += x[i * nVertices + j];
                expr2 += x[j * nVertices + i];
              }
            }
          }
          add(expr1 >= 1.0).end();
          expr1.end();
          add(expr2 >= 1.0).end();
          expr2.end();
        }
        break;
      }
    }

    if (numVisited < nVertices) {
      tour.clear();
      fill(subtour.begin(), subtour.end(), false);
      currentCity = -1;
      for (int i = 0; i < nVertices; i++) {
        if (visited[i] == false) {
          currentCity = i;
          tour.push_back(currentCity);
          visited[i] = true;
          subtour[i] = true;
          numVisited++;
          subtourLength = 1;
          initialCity = currentCity;
          break;
        }
      }
    } else
      break;
  }
  if (opt.dwtype == "atsp")
    return;
  if (subtourLength != nVertices)
    return;

  for (int origin = 0; origin < nVertices; origin++) {
    for (int dest = (origin + 1); dest < nVertices; dest++) {
      if (I->D[d->map[tour[dest]]][d->map[tour[origin]]]) {
        IloExpr expr(getEnv());
        for (int i = origin; i < dest; i++)
          expr += x[tour[i] * nVertices + tour[i + 1]];
        add(expr <= dest - origin - 1).end();
        expr.end();
      }
    }
  }
}

int sop(const Instance *I, dw_t *d, int nVertices, Task initial) {

  int retValue;
  IloEnv env;
  IloModel m(env);
  IloNumArray values(env, nVertices * nVertices);
  IloBoolVarArray x(env, nVertices * nVertices);

  IloExpr expr(env);
  for (int i = 0; i < nVertices; i++)
    for (int j = 0; j < nVertices; j++)
      expr += d->assigncost[i][j] * x[i * nVertices + j];
  m.add(IloMinimize(env, expr));
  expr.end();
  for (int i = 0; i < nVertices; i++) {
    IloExpr expr(env);
    for (int j = 0; j < nVertices; j++)
      expr += x[i * nVertices + j];
    m.add(expr == 1);
    expr.end();
  }
  for (int j = 0; j < nVertices; j++) {
    IloExpr expr(env);
    for (int i = 0; i < nVertices; i++)
      expr += x[i * nVertices + j];
    m.add(expr == 1);
    expr.end();
  }
  IloCplex solver(m);
  solver.setParam(IloCplex::Param::Threads, 1);
  solver.use(myLazyCallback(env, x, values, I, d, nVertices, initial));

  solver.solve();
  retValue = floor(solver.getObjValue() + RC_EPS);
  env.end();
  return retValue;
}
