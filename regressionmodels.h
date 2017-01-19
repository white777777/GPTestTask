#ifndef REGRESSIONMODELS_H
#define REGRESSIONMODELS_H

#include <vector>

#include "functions.h"
#include "taskdata.h"
#include <eigen3/Eigen/Dense>

/// struct for interacting with solver
struct WorkingSet
{
  WorkingSet(){};
  WorkingSet(size_t taskSize, size_t nParams)
  : J(Eigen::MatrixXd::Zero(taskSize, nParams))
  , yMinusF(taskSize)
  {
  }
  /// Task jacobi matrix \f$ J \f$. See (1)
  Eigen::MatrixXd J;
  /// \f$ y-f \f$ in terms of (1). Where \f$ Y=dQ/dTt = Q_{ij}/t_{ij} \f$
  Eigen::VectorXd yMinusF;
};

class IRegressionModel
{
public:
  virtual bool IsReady() const = 0;
  virtual Eigen::VectorXd GenParams0Vec() = 0;
  virtual WorkingSet InitWorkingSet() = 0;
  virtual void CalcValue(const Eigen::VectorXd & params, WorkingSet & ws) = 0;
  virtual void NormalizeParams(Eigen::VectorXd & params) = 0;
};

/// Regression model
/// \f$ ln q_i (t_j) = ln q_{0i} + ln(t_j)+\ksi_{ij} \f$
/// See also
/// (1) http://www.machinelearning.ru/wiki/index.php?title=Нелинейная_регрессия
template<class TFunc>
class RegressionModelLn: public IRegressionModel
{
public:
  //typedef TemplateFunc TFunc;
  struct OptimizedHoleData
  {
    std::vector<double> sumT;
    std::vector<double> qDivT;
  };
  
  struct OptimizedTaskData
  {
    std::vector<OptimizedHoleData> holes;
  };
  
private:

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
  
  RegressionModelLn(const TaskData& taskData)
  //please be carefull with initialization order
  : _taskData(taskData)
  , _oTD(fillOptimizedHoleData(_taskData))
  , _taskSize(TaskDataHelper::GetTaskSize(_taskData))
  , _nQParams(_taskData.holes.size())
  , _nFuncParams(_func.nParams)
  , _nParams(_nQParams + _nFuncParams)
  {
  }
  
  bool IsReady() const
  {
    if(_taskData.holes.size()==0
      && _taskData.holes.size() == _oTD.holes.size()
    )
      return false;
      return true;
  }
  
  Eigen::VectorXd GenParams0Vec()
  {
    Eigen::VectorXd params(_nParams);
    for(size_t i = 0; i<_nQParams; ++i)
      params[i] = _oTD.holes[i].qDivT[0];
    for(size_t i = 0; i<_nFuncParams; ++i)
      params[_taskData.holes.size() + i] = _func.GetDefaultParam(i);
    return params;
  }
  
  WorkingSet InitWorkingSet()
  {
    return WorkingSet(_taskSize, _nParams);
  }

  void CalcValue(const Eigen::VectorXd& params, WorkingSet& ws)
  {
    //const Eigen::Map<const TFunc::VParams> funcParams(&params[_nQParams], _nFuncParams);
    const Eigen::Map<const Eigen::VectorXd> funcParams(&params[_nQParams], _nFuncParams);
    size_t it = 0;
    for(size_t i = 0; i<_taskData.holes.size(); ++i)
    {
      const HoleData& holeData = _taskData.holes[i];
      const OptimizedHoleData& optHoleData = _oTD.holes[i];
      for(size_t j = 0; j<holeData.ts.size();++j)
      {
        // Find the way to minimize task size
        ws.J(it, i) = 1.0/params[i];
        
        const double tFromStart = optHoleData.sumT[j];
        for(size_t iParam = 0; iParam<_nFuncParams; ++iParam)
          ws.J(it, _nQParams+iParam) = _func.calcDFDIParamDivFT(iParam, funcParams, tFromStart);
        
        const double qDivTVal = optHoleData.qDivT[j];
        if(abs(qDivTVal) > 0)
          ws.yMinusF[it] = log(qDivTVal) - log(params[i]) - _func.calcLnFT(funcParams, tFromStart);
        else
        {
          // TODO: filter this line from input task data
          ws.yMinusF[it] = 0.0;
          for(size_t l = 0; l<_nParams; ++l)
            ws.J(it, l) = 0.0;
        }
        ++it;
      }
    }
  }
  
  void NormalizeParams(Eigen::VectorXd& params)
  {
    size_t nClip = 0;
    for(size_t i = 0; i< _nQParams; ++i)
    {
      const double eps = 1e-16;
      if(params[i] < eps)
      {
        params[i] = eps;
        ++nClip;
      }
    }
    for(size_t i = 0; i< _nFuncParams; ++i)
    {
      if(params[i+_nQParams] < _func.GetParamLowerLimits(i)) 
      {
        params[i+_nQParams] = _func.GetParamLowerLimits(i);
        ++nClip;
      }
      if(params[i+_nQParams] > _func.GetParamUpperLimits(i)) 
      {
        params[i+_nQParams] = _func.GetParamUpperLimits(i);
        ++nClip;
      }
    }
    //if(nClip>0)
    //  std::cout<<"Warning: "<<nClip<<" params out of range"<< std::endl;
  }
private:
  OptimizedTaskData fillOptimizedHoleData(const TaskData& taskData)
  {
    OptimizedTaskData oTD;
    oTD.holes.resize(taskData.holes.size());
    for(size_t i = 0; i< taskData.holes.size(); ++i)
    {
      OptimizedHoleData& hole = oTD.holes[i];
      const HoleData& tHole = taskData.holes[i];
      hole.sumT.resize(tHole.ts.size());
      hole.qDivT.resize(tHole.ts.size());
      for(size_t j = 0; j<tHole.ts.size(); ++j)
      {
        hole.qDivT[j] = tHole.qOils[j]/tHole.ts[j];
        if(j==0)
          hole.sumT[j] = tHole.ts[j]/2;
        else
          //We use half of time (ts) to get more precice Q derivative
          hole.sumT[j] = hole.sumT[j-1] +tHole.ts[j-1]/2 + tHole.ts[j]/2;
      }
    }
    return oTD;
  };
  friend class Tester;
};

typedef RegressionModelLn<Function1> RegressionModelLn1;
typedef RegressionModelLn<Function2> RegressionModelLn2;
typedef RegressionModelLn<Function3> RegressionModelLn3;
typedef RegressionModelLn<Function4> RegressionModelLn4;

#endif // REGRESSIONMODELS_H
