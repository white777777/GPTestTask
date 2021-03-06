#include "analyze.h"

#include <iostream>
#include <cstdlib>
#include "solver.h"
#include "dataimporter.h"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/vector.hpp>

#include "boost_serialization_eigen.h"

#include <fstream>

#include "gnuplot-iostream.h"

#include <algorithm>
#include <numeric>
#include <string>

struct AnalyzeSet
{
  Eigen::VectorXd params;
  std::string name;
  Eigen::VectorXd delta;
  template <typename Archive>  
  void serialize(Archive &ar, const unsigned int version)
  {
    ar & params;
    ar & name;
    ar & delta;
  }
};

void load(std::vector<AnalyzeSet> & results)
{
  std::ifstream ifs("file.txt", std::ofstream::binary);
  boost::archive::text_iarchive ia(ifs);
  ia>>results;
}

void calcAndSave(TaskData taskDataOrig, std::vector<AnalyzeSet> & results)
{
  std::vector<std::string> names{"QOil1", "QOil2", "QOil4", "QWater1", "QWater2"};
  std::vector<std::unique_ptr<IRegressionModel>> models;
  
  models.push_back(std::make_unique<RegressionModelLn1>(taskDataOrig));
  models.push_back(std::make_unique<RegressionModelLn3>(taskDataOrig));
  models.push_back(std::make_unique<RegressionModelLn4>(taskDataOrig));
  TaskDataHelper::SwapOilWater(taskDataOrig);
  models.push_back(std::make_unique<RegressionModelLn1>(taskDataOrig));
  models.push_back(std::make_unique<RegressionModelLn2>(taskDataOrig));
  
  Solver::SolverParams sp;
  sp.verbose = 2;
  sp.enableNormalizer = true;
  sp.nMaxIter = 25;
  try
  {
    for(size_t i=0; i<models.size(); ++i)
    {
      std::cout<<i<<std::endl;
      Solver solver(std::move(models[i]));
      solver.SolverInit(sp);
      if(!solver.Solve())
        //throw std::logic_error("Can' solve");
        std::cout<<"not solved"<<std::endl;
      results.push_back({solver.GetResult(), names[i], solver.GetWorkingSet().yMinusF});
    }
  }
  catch(std::logic_error &e)
  {
    std::cout<<e.what()<<std::endl;
  }
  
  std::cout<<"Analyze calculation is finished"<<std::endl;
  
  {
    std::cout<<"Writing data to file"<<std::endl;//TODO: read data Q(t)
    
    std::ofstream ofs("file.txt", std::ofstream::binary);
    boost::archive::text_oarchive oa(ofs);
    oa<<results;
  }
}

typedef std::vector<double> dvec;

void Analyzer::Analyze(const std::string & filename)
{
  
  CSVDataImporter dataImporter;
  TaskData taskDataOrig = dataImporter.read(filename);
  OptimizedTaskData oTD(taskDataOrig);

  std::vector<AnalyzeSet> results;
  bool isSuccessffullRead = false;
  try
  {
    load(results);
    isSuccessffullRead = true;
  }
  catch(std::exception &e)
  {
    std::cout<<e.what()<<std::endl;
  }
  if(!isSuccessffullRead)
  {
    calcAndSave(taskDataOrig, results);
  }

  // get timeVec and maxT
  dvec time;
  double maxT = 0;
  for(auto & hole: oTD.holes)
  {
    double t = hole.sumT.back();
    maxT = maxT>t ? maxT : t;
  }
  const double tStep = 100;
  for(size_t i = 0; i<maxT/tStep; ++i)
    time.push_back(i*tStep);
  
  
  auto drawFunc = [&time](auto func, Eigen::VectorXd params,  std::string name){
    std::cout<<params.transpose()<<std::endl;
    Gnuplot gp;
    gp<<"set terminal postscript eps enhanced color font 'Helvetica,10'"<<std::endl;
    gp<<"set output '"<<name<<".ps'"<<std::endl;
    dvec val;
    for(size_t i=0; i<time.size(); ++i)
      val.push_back(func(params, time[i]));
    gp << "plot '-' tit '"<<name<<"'"<<std::endl;
    gp.send1d(boost::make_tuple(time, val));
  };

  drawFunc(Function1::CalcFT, results[0].params.tail(Function1::nParams), std::string("func_")+results[0].name);
  drawFunc(Function2::CalcFT, results[1].params.tail(Function2::nParams), std::string("func_")+results[1].name);
  drawFunc(Function4::CalcFT, results[2].params.tail(Function4::nParams), std::string("func_")+results[2].name);
  drawFunc(Function1::CalcFT, results[3].params.tail(Function1::nParams), std::string("func_")+results[3].name);
  drawFunc(Function2::CalcFT, results[4].params.tail(Function2::nParams), std::string("func_")+results[4].name);

  auto saveToCSV = [](const auto & arr, std::string name){
    std::ofstream ofs(name);
    for(size_t i = 0; i< arr.size(); ++i)
      ofs<<arr[i]<<","<<std::endl;
  };
  

  for(const AnalyzeSet & result: results)
  {
    std::ofstream ofs(result.name + "_params.txt");
    ofs<<result.params;
    saveToCSV(result.delta, std::string("delta_") + result.name +".csv");
  }
  
  
};