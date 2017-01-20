#include "tester.h"

void Tester::testSolver()
{
  const std::vector<size_t> sizes{300, 200};
  Eigen::VectorXd funcParams(1);
  funcParams<<0.001;
  Eigen::VectorXd q0iParams(sizes.size());
  q0iParams<<2, 4;
  Eigen::VectorXd params(q0iParams.size() + funcParams.size());
  params<<q0iParams, funcParams;
  Solver solver(std::make_unique<RegressionModelLn1>(generateTaskData<Function1>(sizes, q0iParams, funcParams)));
  solver.SolverInit();
  if(!solver.Solve())
  {
    std::cout<<"Solution not found"<<std::endl;
    tassert(false);
  }
  
  Eigen::VectorXd res = solver.GetResult();
  std::cout<<"Aaand solution iiiis..."<<res.transpose()<<std::endl;
  Eigen::VectorXd delta = res - Eigen::VectorXd(params);
  //works only with integrate defined //tassert(delta.lpNorm<Eigen::Infinity>()<1e-1);
  std::cout<<"test passed"<<std::endl;
}

void Tester::testRealWorld()
{
  CSVDataImporter dataImporter;
  Solver solver(std::make_unique<RegressionModelLn1>(dataImporter.read("/mnt/windows/Users/user/Documents/projects/GPTestTask/taskData.csv")));
  solver.SolverInit();
  if(!solver.Solve())
  {
    std::cout<<"Solution not found"<<std::endl;
    tassert(false);
  }
  std::cout<<"Result model params"<<solver.GetResult().transpose()<<std::endl;
  std::cout<<"test passed"<<std::endl;
}



template<class TFunc>
class fhlp
{
public:
  static Eigen::VectorXd GetDefaultParams()
  {
    size_t n = TFunc::nParams;
    Eigen::VectorXd params(n);
    for(size_t i=0; i< n; ++i)
      params[i] = TFunc::GetDefaultParam(i);
    return params;
  }
};

void Tester::Test()
{
  try
  {
    testExactSolution<Function1>(fhlp<Function1>::GetDefaultParams(), true);
    testExactSolution<Function2>(fhlp<Function2>::GetDefaultParams());
    testExactSolution<Function3>(fhlp<Function3>::GetDefaultParams());
    testExactSolution<Function4>(fhlp<Function4>::GetDefaultParams());
    testSolver();
    testRealWorld();
  }
  catch(...)
  {
    std::cout<<"Test error"<<std::endl;
  }
}
void Tester::tassert(bool val)
  {
    if(!val)
    {
      std::cout<<"error"<<std::endl;
      throw std::logic_error("Assertion failed");
    }
  }
