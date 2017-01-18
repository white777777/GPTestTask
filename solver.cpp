#include "solver.h"
void Solver::SolverInit()
{
  if(!_regressionModel.IsReady())
    throw std::invalid_argument("Empty task data");
  _modelParams = _regressionModel.GenParams0Vec();
  _ws = _regressionModel.InitWorkingSet();
  _isInited = true;
}

Eigen::VectorXd Solver::SolveStep()
{
  _regressionModel.CalcValue(_modelParams, _ws);
  return _ws.J.colPivHouseholderQr().solve(_ws.yMinusF);
}

Eigen::VectorXd Solver::GetResult() const
{
  return _modelParams;
}

bool Solver::Solve()
{
  if(!_isInited)
    throw std::invalid_argument("Solver not initialized");
  for(size_t nIter = 0; nIter<nMaxIter; ++nIter)
  {
    Eigen::VectorXd deltaParams = SolveStep();
    _modelParams+= deltaParams;
    double diff = sqrt(deltaParams.transpose() *deltaParams);
    if(diff<eps)
    {
      _isInited = false;
      return true;
    }
  }
  return false;
}

Solver::Solver(RegressionModelLn&& rm)
: _regressionModel(std::move(rm))
{
}
