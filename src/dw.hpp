#pragma once

#include <vector>

#include "instance.hpp"
#include "lap.hpp"
#include "models.hpp"
#include "util.hpp"

#define UNKNOWN -1
#define OUT 0
#define IN 1

struct dwTask_t {
  Task id;
  double price;
  double priceRatio;
  void clear();
};

struct dw_t {
  std::vector<dwTask_t> tasks;    
  std::vector<Task> order;        
  std::vector<std::vector<short int>> bestSols; 
  std::vector<short int> state;   
  std::vector<Task> pi;           
  double obj;                     
  double bestObj;                 
  int itCount;                    
  Time tAcc;                      
  Task *map; 
  cost **assigncost, *u, *v;
  row *colsol;
  col *rowsol;

  void clear(unsigned);
  void allocate_memory(unsigned);
  void free_memory(unsigned);
  unsigned eSt; 
  unsigned lSt;
};

struct sched_t {   
  Task first;      
  Task last;       
  int acc;         
  int accS;        
  int depth;       
  bool isFeasible; 
};

using namespace sualbsp;

unsigned short int dw(Instance &, const ModelOptions&);
