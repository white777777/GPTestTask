#include "taskdata.h"

#include <iostream>
#include <algorithm>

size_t TaskDataHelper::GetTaskSize(const TaskData& taskData)
{
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
    //std::cout<<"  Hole '"<< holeData.name << "' size: " << holeData.ts.size() << ". Stat: T: "<< sumT<< ", QOil: "<< sumQOil<< std::endl;
    taskDataSize += holeData.ts.size();
  }
  std::cout<<"Task size: "<< taskDataSize<<", number holes: "<< taskData.holes.size() <<std::endl;
  return taskDataSize;
}
  
void TaskDataHelper::StripTaskData(TaskData& taskData, size_t iHole, size_t nHole, size_t nQ)
{
  if(iHole>=taskData.holes.size())
  {
    taskData = TaskData();
    return;
  }
  size_t nn=taskData.holes.size() - iHole;
  nHole = nHole>nn? nn:nHole;
  for(size_t j = 0; j<nHole; ++j)
  {
    taskData.holes[j] = taskData.holes[j+iHole];
    size_t nnQ = taskData.holes[j].ts.size();
    if(nQ<nnQ)
    {
      taskData.holes[j].ts.resize(nQ);
      taskData.holes[j].qOils.resize(nQ);
    }
  }
  taskData.holes.resize(nHole);
}

void TaskDataHelper::SwapOilWater(TaskData& taskData)
{
  for(size_t iHole = 0; iHole<taskData.holes.size(); ++iHole)
  {
    std::swap(taskData.holes[iHole].qOils, taskData.holes[iHole].qWaters);
  }
}
