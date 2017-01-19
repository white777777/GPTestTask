#ifndef SOLVER_H
#define SOLVER_H

#include "regressionmodels.h"
#include <memory>

/// Regression task solver
class Solver
{
private:
  //Solver params
  double eps = 1e-4;
  size_t nMaxIter = 10000;
  
  //Solver state;
  bool _isInited = false;

  //input data
  std::unique_ptr<IRegressionModel> _regressionModel;
  
  //working set
  WorkingSet _ws;
  Eigen::VectorXd _modelParams;
public:
  Solver(std::unique_ptr<IRegressionModel> rm);
  
  void SolverInit();
  /// One solve step. Genereates and solves SLE from regression model
  Eigen::VectorXd  SolveStep();
  /// Returns result model params
  Eigen::VectorXd GetResult() const;
  
  WorkingSet GetWorkingSet() const;

  /// Solve problem
  /// returns true if solution found
  bool Solve();
};


#endif // SOLVER_H
