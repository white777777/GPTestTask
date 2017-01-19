#include "tester.h"
void tassert(bool val)
{
  if(!val)
  {
    std::cout<<"error"<<std::endl;
    throw;
  }
}

TaskData generateTaskData(const FunctionRef& f, const std::vector<size_t> sizes, const Eigen::VectorXd& params)
{
  tassert(sizes.size() == size_t(params.size()-1));
  srand (123455);
  const FunctionRef::VParams funcParams{params[sizes.size()]};
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
      double deltaT = rand()%100+500;
      double f1 = f.calcIntF(funcParams, t);
      double f2 = f.calcIntF(funcParams, t+deltaT);
      hole.qOils[i] = params[iHole]*(f2-f1);
      hole.ts[i] = deltaT;
      t+=deltaT;
    }
  }
  return taskData;
}

void Tester::testExactSolution()
{
  std::cout<<"testExactSolution"<<std::endl;
  FunctionRef f;
  Eigen::Matrix<double,1, 4> params{1, 8, 3, 0.0001};
  const std::vector<size_t> sizes{100, 200, 100};
  RegressionModelLn rm(generateTaskData(f, sizes, params));
  const double eps = 1e-8;
  for(size_t i = 1; i<rm._oTD.holes[0].qDivT.size(); ++i)
  {
    tassert(rm._oTD.holes[0].qDivT[i] - eps<rm._oTD.holes[0].qDivT[i-1]);
    tassert(rm._oTD.holes[0].qDivT[i] + eps>0);
  }
  WorkingSet ws = rm.InitWorkingSet();
  rm.CalcValue(params, ws);
  tassert(ws.J(0,0) == 1/params[0]);
  tassert(ws.J(0,1) == 0.0);
  for(int i = 0; i< ws.yMinusF.size(); ++i)
    tassert(ws.yMinusF[i]<0.1);
  if(false)
  {
    std::cout<<ws.yMinusF<<std::endl<<std::endl;
    std::cout<<ws.J<<std::endl<<std::endl;
  }
  
  Eigen::VectorXd paramsDelta = ws.J.colPivHouseholderQr().solve(ws.yMinusF);
  for(int i = 0; i< paramsDelta.size(); ++i)
    tassert(paramsDelta[i]<0.1);
  std::cout<<"test passed"<<std::endl;
}

void Tester::testSolver()
{
  FunctionRef f;
  Eigen::Matrix<double,1, 3> params{2, 4, 0.001};
  const std::vector<size_t> sizes{300, 200};
  Solver solver(RegressionModelLn(generateTaskData(f, sizes, params)));
  solver.SolverInit();
  if(!solver.Solve())
  {
    std::cout<<"Solution not found"<<std::endl;
    tassert(false);
  }
  
  Eigen::VectorXd res = solver.GetResult();
  Eigen::VectorXd delta = res - Eigen::VectorXd(params);
  tassert(sqrt(delta.transpose()*delta)<1e-1);
  std::cout<<"Aaand solution iiiis..."<<res.transpose()<<std::endl;
  std::cout<<"test passed"<<std::endl;
}

void Tester::testRealWorldIterative()
{
  CSVDataImporter dataImporter;
  TaskData taskDataOrig = dataImporter.read("/mnt/windows/Users/user/Documents/projects/GPTestTask/taskData.csv");
  
  Eigen::VectorXd prevResult;
  {
    TaskData taskData = taskDataOrig;
    TaskDataHelper::StripTaskData(taskData, 0, 20000000, 3);
    RegressionModelLn rm(taskData);
    Solver solver(std::move(rm));
    solver.SolverInit();
    if(!solver.Solve())
    {
      std::cout<<"Solution not found"<<std::endl;
    }
    prevResult = solver.GetResult();
    std::cout<<"Result model params"<<prevResult.transpose()<<std::endl;
  }
  {
    TaskData taskData = taskDataOrig;
    TaskDataHelper::StripTaskData(taskData, 0, 20000000, 20);
    RegressionModelLn rm(taskData);
    Solver solver(std::move(rm));
    solver.SolverInit();
    solver._modelParams = prevResult;
    if(!solver.Solve())
    {
      std::cout<<"Solution not found"<<std::endl;
    }
    prevResult = solver.GetResult();
    std::cout<<"Result model params"<<prevResult.transpose()<<std::endl;
  }
  {
    TaskData taskData = taskDataOrig;
    TaskDataHelper::StripTaskData(taskData, 0, 20000000, 2000000000);
    RegressionModelLn rm(taskData);
    Solver solver(std::move(rm));
    solver.SolverInit();
    solver._modelParams = prevResult;
    if(!solver.Solve())
    {
      std::cout<<"Solution not found"<<std::endl;
    }
    prevResult = solver.GetResult();
    std::cout<<"Result model params"<<prevResult.transpose()<<std::endl;
  }
}

void Tester::testRealWorld()
{
  CSVDataImporter dataImporter;
  TaskData taskDataOrig = dataImporter.read("/mnt/windows/Users/user/Documents/projects/GPTestTask/taskData.csv");
  RegressionModelLn rm(taskDataOrig);
  Solver solver(std::move(rm));
  solver.SolverInit();
  if(!solver.Solve())
    std::cout<<"Solution not found"<<std::endl;
    std::cout<<"Result model params"<<solver.GetResult().transpose()<<std::endl;
}

void Tester::Test()
{
  try
  {
    //testExactSolution();
    //testSolver();
    testRealWorld();
  }
  catch(...)
  {
    std::cout<<"Test error"<<std::endl;
  }
}
