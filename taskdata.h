#ifndef TASKDATA_H
#define TASKDATA_H
#include <vector>
#include <string>
#include <iostream>

/// Statistical input data which describes task
//TODO: read data Q(t)
//TODO: convert data to Q'(t)
struct HoleData
{
  /// hole name
  std::string name;
  /// hours per month
  std::vector<double> ts;
  /// Oil per month for every borehole
  std::vector<double> qOils;
  /// Water per month for every borehole
  std::vector<double> qWaters;
};

struct TaskData
{
  std::vector<HoleData> holes;
};

class TaskDataHelper
{
public:
  static size_t GetTaskSize(const TaskData & taskData)
  {
    std::cout<<"Task size:"<<std::endl;
    size_t taskDataSize = 0;
    for(const HoleData & holeData: taskData.holes)
    {
      double sumQOil = 0;
      double sumT = 0;
      for(size_t i = 0; i<holeData.ts.size(); ++i)
      {
        sumQOil += holeData.qOils[i];
        sumT += holeData.ts[i];
      }
      std::cout<<"  Hole '"<< holeData.name << "' size: " << holeData.ts.size() << ". Stat: T: "<< sumT<< ", QOil: "<< sumQOil<< std::endl;
      taskDataSize += holeData.ts.size();
    }
    std::cout<<"Number holes: "<< taskData.holes.size() <<". Whole stat size: "<< taskDataSize<<std::endl;
    return taskDataSize;
  }
};

#endif // TASKDATA_H
