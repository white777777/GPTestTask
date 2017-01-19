#include "solver.h"
#include <eigen3/Eigen/IterativeLinearSolvers>
#include <cmath>
#include <limits>
#include <iostream>

void Solver::SolverInit()
{
  if(!_regressionModel->IsReady())
    throw std::invalid_argument("Empty task data");
  _modelParams = _regressionModel->GenParams0Vec();
  _ws = _regressionModel->InitWorkingSet();
  _isInited = true;
}

Eigen::VectorXd Solver::SolveStep()
{
  using namespace Eigen;
  _regressionModel->NormalizeParams(_modelParams);
  _regressionModel->CalcValue(_modelParams, _ws);
  //Eigen::MatrixXd A = _ws.J.transpose()*_ws.J;
  //Eigen::MatrixXd b = _ws.J.transpose()*_ws.yMinusF;
  //return A.householderQr().solve(b);
  
  MatrixXd A = _ws.J.transpose()*_ws.J;
  MatrixXd b = _ws.J.transpose()*_ws.yMinusF;
  ConjugateGradient<MatrixXd, Lower|Upper> cg;
  cg.compute(A);
  return cg.solve(b);
  
  //return _ws.J.householderQr().solve(_ws.yMinusF);
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
    _modelParams += deltaParams;

    double diff1 = deltaParams.lpNorm<Eigen::Infinity>();
    double diff2 = _ws.yMinusF.lpNorm<Eigen::Infinity>();
    std::cout<<"Step: "<<nIter<<" diff1: "<<diff1<<" Y-F ln: "<<diff2<<std::endl;
    if(diff1<eps || diff2<eps)
    {
      _isInited = false;
      return true;
    }
  }
  return false;
}

Solver::Solver(std::unique_ptr<IRegressionModel> rm)
: _regressionModel(std::move(rm))
{
}

WorkingSet Solver::GetWorkingSet() const
{
  return _ws;
}
