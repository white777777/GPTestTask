#ifndef REGRESSIONMODELS_H
#define REGRESSIONMODELS_H

#include <vector>
#include <utility>

#include "functions.h"
#include "taskdata.h"
#include <eigen3/Eigen/Dense>

/// struct for interacting with solver
struct WorkingSet
{
  WorkingSet(){};
  WorkingSet(size_t taskSize, size_t nParams)
  : J(taskSize, nParams)
  , yMinusF(taskSize)
  {
  }
  /// Task jacobi matrix \f$ J \f$. See (1)
  Eigen::MatrixXd J;
  /// \f$ y-f \f$ in terms of (1). Where \f$ Y=dQ/dTt = Q_{ij}/t_{ij} \f$
  Eigen::VectorXd yMinusF;
};

/// Regression model
/// \f$ ln q_i (t_j) = ln q_{0i} + ln(t_j)+\ksi_{ij} \f$
/// See also
/// (1) http://www.machinelearning.ru/wiki/index.php?title=Нелинейная_регрессия
//TODO: template<class TFunc>
class RegressionModelLn
{
public:
  typedef FunctionRef TFunc;  
private:
  struct OptimizedHoleData
  {
    std::vector<double> sumT;
    std::vector<double> qDivT;
  };
  
  struct OptimizedTaskData
  {
    std::vector<OptimizedHoleData> holes;
  };

  TFunc _func;
  /// Input statistical data
  TaskData _taskData;
  /// Task data optimized for our purposes
  OptimizedTaskData _oTD;
  /// Task size: size of statistical data
  const size_t _taskSize = 0;
  const size_t _nQParams = 0;
  const size_t _nFuncParams = 0;
  const size_t _nParams = 0;
public:  
  RegressionModelLn() = delete;
  RegressionModelLn(const TaskData& taskData);
  bool IsReady() const;
  
  Eigen::VectorXd GenParams0Vec();
  WorkingSet InitWorkingSet()
  {
    return WorkingSet(_taskSize, _nParams);
  }
  void CalcValue(const Eigen::VectorXd & params, WorkingSet & ws);
private:
  OptimizedTaskData fillOptimizedHoleData(const TaskData& taskData);
  friend class Tester;
};

#endif // REGRESSIONMODELS_H
