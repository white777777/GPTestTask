#ifndef SOLVER_H
#define SOLVER_H

#include <cmath>
#include <array>
#include <limits>

#include "regressionmodels.h"

/// Regression task solver
//template<class TRegr>
class Solver
{
private:
  typedef RegressionModelLn TRegr;
  //Solver params
  double eps = 1e-4;
  size_t nMaxIter = 10000;
  
  //Solver state;
  bool _isInited = false;
  
  //input data
  TRegr _regressionModel;
  
public:
  //working set
  WorkingSet _ws;
  Eigen::VectorXd _modelParams;
  
public:
  Solver(RegressionModelLn&& rm);
  
  void SolverInit();
  /// One solve step. Genereates and solves SLE from regression model
  Eigen::VectorXd  SolveStep();
  /// Returns result model params
  Eigen::VectorXd GetResult() const;
  /// Solve problem
  /// returns true if solution found
  bool Solve();
};


#endif // SOLVER_H
