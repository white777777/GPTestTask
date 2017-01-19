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
  static size_t GetTaskSize(const TaskData & taskData);
  
  static void StripTaskData(TaskData & taskData, size_t i, size_t n, size_t nQ);
};

#endif // TASKDATA_H
