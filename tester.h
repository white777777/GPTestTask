#ifndef TESTER_H
#define TESTER_H

#include <iostream>
#include <cstdlib>
#include "solver.h"
#include "dataimporter.h"
#include <memory>
#include <ostream>

/**
 * Overloaded output operator for vectors.
 * taken from http://stackoverflow.com/a/31394136,
 * last access 11/26/2016
 */
template<typename T>
std::ostream &operator<<(std::ostream &ostr, const std::vector<T> &v) {
  ostr << "[";
  size_t last = v.size() - 1;
  for (size_t i = 0; i < v.size(); ++i) {
    ostr << v[i];
    if (i != last)
      ostr << ", ";
  }
  ostr << "]";
  return ostr;
}

class Tester
{
  void tassert(bool val);

  template<class TFunc>
  TaskData generateTaskData(const std::vector<size_t> sizes, const Eigen::VectorXd& params, const Eigen::VectorXd& funcParams)
  {
    tassert(sizes.size() == size_t(params.size()));
    srand (123455);
    TaskData taskData;
    taskData.holes.resize(sizes.size());
    for(size_t iHole = 0; iHole<sizes.size(); ++iHole)
    {
      const size_t n = sizes[iHole];
      HoleData& hole = taskData.holes[iHole];
      hole.qOils.resize(n);
      hole.ts.resize(n);
      double t = 0.0;
      for(size_t i = 0; i<n; ++i)
      {
        double deltaT = 500;// rand()%100+500;
        //double f1 = TFunc::calcIntF(funcParams, t);
        //double f2 = TFunc::calcIntF(funcParams, t+deltaT);
        //hole.qOils[i] = params[iHole]*(f2-f1);
        double f1 = TFunc::CalcFT(funcParams, t);
        double f2 = TFunc::CalcFT(funcParams, t+deltaT);
        hole.qOils[i] = params[iHole]*(f2+f1)/2.0*deltaT;
        hole.ts[i] = deltaT;
        t+=deltaT;
      }
    }
    return taskData;
  }
  
  template<class TFunc>
  void testExactSolution(const Eigen::VectorXd & funcParams, bool p = false)
  {
    using namespace Eigen;
    std::cout<<"testExactSolution"<<std::endl;

    const std::vector<size_t> sizes{10, 20, 10};
    VectorXd q0iParams(sizes.size());
    q0iParams<<1, 8, 3;
    VectorXd params(q0iParams.size() + funcParams.size());
    params<<q0iParams, funcParams;

    if(p)
      std::cout<<"params:"<<params.transpose()<<std::endl<<std::endl;
    
    RegressionModelLn<TFunc> rm(generateTaskData<TFunc>(sizes, q0iParams, funcParams));
    for(size_t i = 1; i<rm._oTD.holes[0].qDivT.size(); ++i)
    {
      const double eps = 1e-8;
      tassert(rm._oTD.holes[0].sumT[i] > rm._oTD.holes[0].sumT[i-1]);
      tassert(rm._oTD.holes[0].qDivT[i] - eps<rm._oTD.holes[0].qDivT[i-1]);
      tassert(rm._oTD.holes[0].qDivT[i] + eps>0);
    }
    if(p)
    {
      std::cout<<"sumT: "<<rm._oTD.holes[0].sumT<<std::endl<<std::endl;
      std::cout<<"qDivT: "<<rm._oTD.holes[0].qDivT<<std::endl<<std::endl;
    }
    
    WorkingSet ws = rm.InitWorkingSet();
    rm.CalcValue(params, ws);
    if(p)
    {
      std::cout<<"y-f:"<<ws.yMinusF<<std::endl<<std::endl;
      std::cout<<"J: "<<ws.J<<std::endl;
    }
    tassert(ws.J(0,0) == 1/params[0]);
    tassert(ws.J(0,1) == 0.0);
    for(int i = 0; i< ws.yMinusF.size(); ++i)
      tassert(ws.yMinusF[i]<0.1);
    
    VectorXd paramsDelta = ws.J.colPivHouseholderQr().solve(ws.yMinusF);
    if(p)
      std::cout<<"paramsDelta: "<<paramsDelta.transpose()<<std::endl<<std::endl;
    //for(int i = 0; i< paramsDelta.size(); ++i)
      //tassert(paramsDelta[i]<0.1);
    std::cout<<"test passed"<<std::endl;
  }
  
  void testExactSolutionPrint();
  void testSolver();
  void testRealWorldIterative();
  void testRealWorld();
public:
  void Test();
};


#endif // TESTER_H
