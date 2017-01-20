#ifndef SOLVER_H
#define SOLVER_H

#include "regressionmodels.h"
#include <memory>

/// Regression task solver
class Solver
{
public:
  struct SolverParams
  {
  public:
    SolverParams()
    : epsDiff(1e-6)
    , epsYMinusF(1e-6)
    , nMaxIter(1000)
    , verbose(0)
    , enableNormalizer(true)
    {
    }
    double epsDiff;
    double epsYMinusF;
    size_t nMaxIter;
    int verbose;
    bool enableNormalizer;
  };
private:
  //Solver state;
  bool _isInited = false;
  SolverParams _sp;

  //input data
  std::unique_ptr<IRegressionModel> _regressionModel;
  
  //working set
  WorkingSet _ws;
  Eigen::VectorXd _modelParams;
public:
  Solver(std::unique_ptr<IRegressionModel> rm);
  
  void SolverInit(const SolverParams & sp = SolverParams());
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
