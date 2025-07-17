
#define BIG 1000000
#if !defined TRUE
#define TRUE 1
#endif
#if !defined FALSE
#define FALSE 0
#endif


typedef int row;
#define ROW_TYPE NPY_INT
typedef int col;
#define COL_TYPE NPY_INT
typedef int cost;
#define COST_TYPE NPY_DOUBLE
typedef int boolean;


cost lap(int dim, cost **assigncost, col *rowsol, row *colsol, cost *u, cost *v);
void examplelap(void);
int checklap(int dim, cost **assigncost, col *rowsol, row *colsol, cost *u, cost *v);
