#include <iostream>
#include "solver.h"
#include "dataimporter.h"
#include "tester.h"

int main(int argc, char **argv) 
{
  //TODO: move test to anotother executable
  Tester tester;
  tester.Test();
  return 0;
  try
  {
    CSVDataImporter dataImporter;
    RegressionModelLn rm(dataImporter.read("/mnt/windows/Users/user/Documents/projects/GPTestTask/taskData.csv"));
    Solver solver(std::move(rm));
    solver.SolverInit();
    if(!solver.Solve())
    {
      std::cout<<"Solution not found"<<std::endl;
    }
    std::cout<<"Result model params"<<solver.GetResult().transpose()<<std::endl;
  }
  catch(std::exception &e)
  {
    std::cout<<"Exception:"<<e.what()<<std::endl;
    return 1;
  }
  catch(...)
  {
    return 1;
  }
  return 0;
}
