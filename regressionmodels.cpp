#include "regressionmodels.h"

RegressionModelLn::OptimizedTaskData RegressionModelLn::fillOptimizedHoleData(const TaskData& taskData)
{
  RegressionModelLn::OptimizedTaskData oTD;
  oTD.holes.resize(taskData.holes.size());
  for(size_t i = 0; i< taskData.holes.size(); ++i)
  {
    RegressionModelLn::OptimizedHoleData& hole = oTD.holes[i];
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


RegressionModelLn::RegressionModelLn(const TaskData& taskData)
//please be carefull with initialization order
: _taskData(taskData)
, _oTD(fillOptimizedHoleData(_taskData))
, _taskSize(TaskDataHelper::GetTaskSize(_taskData))
, _nQParams(_taskData.holes.size())
, _nFuncParams(_func.nParams)
, _nParams(_nQParams + _nFuncParams)
{
}
  
Eigen::VectorXd RegressionModelLn::GenParams0Vec()
{
  Eigen::VectorXd params(_nParams);
  for(size_t i = 0; i<_nQParams; ++i)
    params[i] = _oTD.holes[i].qDivT[0];
  for(size_t i = 0; i<_nFuncParams; ++i)
    params[_taskData.holes.size() + i] = _func.GetDefaultParam(i);
  return params;
}

void RegressionModelLn::CalcValue(const Eigen::VectorXd& params, WorkingSet& ws)
{
  const Eigen::Map<const TFunc::VParams> funcParams(&params[_nQParams], _nFuncParams);
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

bool RegressionModelLn::IsReady() const
{
  if(_taskData.holes.size()==0
    && _taskData.holes.size() == _oTD.holes.size()
  )
    return false;
  return true;
}

void RegressionModelLn::NormalizeParams(Eigen::VectorXd& params)
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
  if(nClip>0)
    std::cout<<"Warning: "<<nClip<<" params out of range"<< std::endl;
}

WorkingSet RegressionModelLn::InitWorkingSet()
{
  return WorkingSet(_taskSize, _nParams);
}
