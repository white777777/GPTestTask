#include <iostream>
#include "solver.h"
#include "dataimporter.h"
#include "tester.h"
#include "analyze.h"
#include <boost/program_options.hpp>

int main(int argc, char **argv) 
{
  std::string filename;
  bool isTest;
  bool isAnalyze;
  bool isSolve;
  
  boost::program_options::options_description desc("General options");
  desc.add_options()
  ("help,h", "Show help")
  ("filepath,f", boost::program_options::value<std::string>(&filename)->default_value("../taskData.csv"), "filename")
  ("test,t"    , boost::program_options::bool_switch(&isTest)->default_value(false), "run test")
  ("analyze,a" , boost::program_options::bool_switch(&isAnalyze)->default_value(true), "run analyze")
  ("solve,s"   , boost::program_options::bool_switch(&isSolve)->default_value(false), "run solve")
  ;
  
  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
  boost::program_options::notify(vm);    
  if (vm.count("help")) {
    std::cout << desc << "\n";
    return 0;
  }

  if(isTest)
  {
    Tester tester;
    tester.Test();
  }
  
  if(isAnalyze)
  {
    Analyzer an;
    an.Analyze(filename);
  }
  
  if(isSolve)
  {
    try
    {
      CSVDataImporter dataImporter;
      Solver solver(std::make_unique<RegressionModelLn1>(dataImporter.read(filename)));
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
}
