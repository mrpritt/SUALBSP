#include "lap.hpp"
#include <stdio.h>
#include <stdlib.h>

cost lap(int dim, cost **assigncost, col *rowsol, row *colsol, cost *u, cost *v) {
  boolean unassignedfound;
  row i, imin, numfree = 0, prvnumfree, f, i0, k, freerow, *pred, *free;
  col j, j1, j2 = 0, endofpath, last = 0, low, up, *collist, *matches;
  cost min = 0, h, umin, usubmin, v2, *d;

  free = new row[dim];    
  collist = new col[dim]; 
  matches = new col[dim]; 
  d = new cost[dim];      
  pred = new row[dim];    

  for (i = 0; i < dim; i++)
    matches[i] = 0;

  for (j = dim - 1; j >= 0; j--) 
  {
    min = assigncost[0][j];
    imin = 0;
    for (i = 1; i < dim; i++)
      if (assigncost[i][j] < min) {
        min = assigncost[i][j];
        imin = i;
      }
    v[j] = min;

    if (++matches[imin] == 1) {
      rowsol[imin] = j;
      colsol[j] = imin;
    } else
      colsol[j] = -1; 
  }

  for (i = 0; i < dim; i++)
    if (matches[i] == 0) 
      free[numfree++] = i;
    else if (matches[i] == 1) 
    {
      j1 = rowsol[i];
      min = BIG;
      for (j = 0; j < dim; j++)
        if (j != j1)
          if (assigncost[i][j] - v[j] < min)
            min = assigncost[i][j] - v[j];
      v[j1] = v[j1] - min;
    }

  int loopcnt = 0; 
  do {
    loopcnt++;

    k = 0;
    prvnumfree = numfree;
    numfree = 0; 
    while (k < prvnumfree) {
      i = free[k];
      k++;

      umin = assigncost[i][0] - v[0];
      j1 = 0;
      usubmin = BIG;
      for (j = 1; j < dim; j++) {
        h = assigncost[i][j] - v[j];
        if (h < usubmin) {
          if (h >= umin) {
            usubmin = h;
            j2 = j;
          } else {
            usubmin = umin;
            umin = h;
            j2 = j1;
            j1 = j;
          }
        }
      }

      i0 = colsol[j1];
      if (umin < usubmin)
        v[j1] = v[j1] - (usubmin - umin);
      else             
          if (i0 >= 0) 
      {
        j1 = j2;
        i0 = colsol[j2];
      }

      rowsol[i] = j1;
      colsol[j1] = i;

      if (i0 >= 0) { 
        if (umin < usubmin)
          free[--k] = i0;
        else
          free[numfree++] = i0;
      }
    }
  } while (loopcnt < 2); 

  for (f = 0; f < numfree; f++) {
    freerow = free[f]; 

    for (j = 0; j < dim; j++) {
      d[j] = assigncost[freerow][j] - v[j];
      pred[j] = freerow;
      collist[j] = j; 
    }

    low = 0; 
    up = 0;  
    unassignedfound = FALSE;
    do {
      if (up == low) 
      {
        last = low - 1;

        min = d[collist[up++]];
        for (k = up; k < dim; k++) {
          j = collist[k];
          h = d[j];
          if (h <= min) {
            if (h < min) 
            {
              up = low; 
              min = h;
            }
            collist[k] = collist[up];
            collist[up++] = j;
          }
        }

        for (k = low; k < up; k++)
          if (colsol[collist[k]] < 0) {
            endofpath = collist[k];
            unassignedfound = TRUE;
            break;
          }
      }

      if (!unassignedfound) {
        j1 = collist[low];
        low++;
        i = colsol[j1];
        h = assigncost[i][j1] - v[j1] - min;

        for (k = up; k < dim; k++) {
          j = collist[k];
          v2 = assigncost[i][j] - v[j] - h;
          if (v2 < d[j]) {
            pred[j] = i;
            if (v2 == min) { 
              if (colsol[j] < 0) {
                endofpath = j;
                unassignedfound = TRUE;
                break;
              }
              else {
                collist[k] = collist[up];
                collist[up++] = j;
              }
            }
            d[j] = v2;
          }
        }
      }
    } while (!unassignedfound);

    for (k = 0; k <= last; k++) {
      j1 = collist[k];
      v[j1] = v[j1] + d[j1] - min;
    }

    do {
      i = pred[endofpath];
      colsol[endofpath] = i;
      j1 = endofpath;
      endofpath = rowsol[i];
      rowsol[i] = j1;
    } while (i != freerow);
  }

  cost lapcost = 0;
  for (i = 0; i < dim; i++) {
    j = rowsol[i];
    u[i] = assigncost[i][j] - v[j];
    lapcost = lapcost + assigncost[i][j];
  }

  delete[] pred;
  delete[] free;
  delete[] collist;
  delete[] matches;
  delete[] d;

  return lapcost;
}

int checklap(int dim, cost **assigncost, col *rowsol, row *colsol, cost *u, cost *v) {
  row i;
  col j;
  cost redcost = 0;
  boolean *matched;

  matched = new boolean[dim];

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      if ((redcost = assigncost[i][j] - u[i] - v[j]) < 0) {
        printf("\n");
        printf("negative reduced cost i %d j %d redcost %d\n", i, j, redcost);
        printf("\n\ndim %5d - press key\n", dim);
        return 0;
      }
    }
  }
  for (i = 0; i < dim; i++)
    if ((redcost = assigncost[i][rowsol[i]] - u[i] - v[rowsol[i]]) != 0) {
      printf("\n");
      printf("non-null reduced cost i %d soli %d redcost %d\n", i, rowsol[i], redcost);
      printf("\n\ndim %5d - press key\n", dim);
      return 0;
    }

  for (j = 0; j < dim; j++)
    matched[j] = FALSE;

  for (i = 0; i < dim; i++)
    if (matched[rowsol[i]]) {
      printf("\n");
      printf("column matched more than once - i %d soli %d\n", i, rowsol[i]);
      printf("\n\ndim %5d - press key\n", dim);
      return 0;
    } else
      matched[rowsol[i]] = TRUE;

  for (i = 0; i < dim; i++)
    if (colsol[rowsol[i]] != i) {
      printf("\n");
      printf("error in row solution i %d soli %d solsoli %d\n", i, rowsol[i], colsol[rowsol[i]]);
      printf("\n\ndim %5d - press key\n", dim);
      return 0;
    }
  for (j = 0; j < dim; j++)
    if (rowsol[colsol[j]] != j) {
      printf("\n");
      printf("error in col solution j %d solj %d solsolj %d\n", j, colsol[j], rowsol[colsol[j]]);
      printf("\n\ndim %5d - press key\n", dim);
      return 0;
    }
  delete[] matched;
  return 1;
}

void examplelap(void) {
  cost **assigncost, *u, *v;
  row i, *colsol;
  col j, *rowsol;
  int dim = 10;

  assigncost = new cost *[dim]; 
  for (i = 0; i < dim; i++)
    assigncost[i] = new cost[dim];

  rowsol = new col[dim];
  colsol = new row[dim];
  u = new cost[dim];
  v = new cost[dim];

  for (i = 0; i < dim; i++)
    for (j = 0; j < dim; j++)
      assigncost[i][j] = (i * 975 + j * 331 + 641) % 10000;

  lap(dim, assigncost, rowsol, colsol, u, v);
  checklap(dim, assigncost, rowsol, colsol, u, v);

  delete[] assigncost;
  delete[] rowsol;
  delete[] colsol;
  delete[] u;
  delete[] v;
}
