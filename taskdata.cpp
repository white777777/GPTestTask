#include "taskdata.h"
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
  
void TaskDataHelper::StripTaskData(TaskData& taskData, std::size_t i, std::size_t n, std::size_t nQ)
{
  if(i>=taskData.holes.size())
  {
    taskData = TaskData();
    return;
  }
  size_t nn=taskData.holes.size() - i;
  n = n>nn? nn:n;
  for(size_t j = 0; j<n; ++j)
  {
    taskData.holes[j] = taskData.holes[j+i];
    size_t nnQ = taskData.holes[j].ts.size();
    if(nQ<nnQ)
    {
      taskData.holes[j].ts.resize(nQ);
      taskData.holes[j].qOils.resize(nQ);
    }
  }
  taskData.holes.resize(n);
}
