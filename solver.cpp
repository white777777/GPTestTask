#include "solver.h"
#include <eigen3/Eigen/IterativeLinearSolvers>
#include <cmath>
#include <limits>
#include <iostream>

void Solver::SolverInit(const Solver::SolverParams & sp)
{
  if(!_regressionModel->IsReady())
    throw std::invalid_argument("Empty task data");
  _sp = sp;
  _modelParams = _regressionModel->GenParams0Vec();
  _ws = _regressionModel->InitWorkingSet();
  _isInited = true;
}

Eigen::VectorXd Solver::SolveStep()
{
  using namespace Eigen;

  _regressionModel->CalcValue(_modelParams, _ws);
  MatrixXd A = _ws.J.transpose()*_ws.J;
  MatrixXd b = _ws.J.transpose()*_ws.yMinusF;
  ConjugateGradient<MatrixXd, Lower|Upper> cg;
  cg.compute(A);
  return cg.solve(b);
}

Eigen::VectorXd Solver::GetResult() const
{
  return _modelParams;
}

bool Solver::Solve()
{
  if(!_isInited)
    throw std::invalid_argument("Solver not initialized");
  
  for(size_t nIter = 0; nIter<_sp.nMaxIter; ++nIter)
  {
    Eigen::VectorXd deltaParams = SolveStep();
    _modelParams += deltaParams;

    if(_sp.verbose > 2)
      std::cout<<"params: "<<_modelParams<<std::endl;
    if(_sp.verbose > 3)
      std::cout<<"y-f: "<< _ws.yMinusF<<std::endl;
    
    if(_sp.enableNormalizer)
    {
      size_t nClip =_regressionModel->NormalizeParams(_modelParams);
      if( nClip>0  && _sp.verbose>1)
        std::cout<<"Warning: "<<nClip<<" params out of range"<< std::endl;
    }
    
    double diff1 = deltaParams.lpNorm<Eigen::Infinity>();
    double diff2 = _ws.yMinusF.lpNorm<Eigen::Infinity>();
    if(_sp.verbose > 0)
      std::cout<<"Step: "<<nIter<<" diff1: "<<diff1<<" Y-F: "<<diff2<<std::endl;
    if(diff1<_sp.epsDiff || diff2<_sp.epsYMinusF)
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
