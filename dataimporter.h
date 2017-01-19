#ifndef DATAIMPORTER_H
#define DATAIMPORTER_H
#include "taskdata.h"
#include <map>
#include <fstream>
#include <sstream>

/// Imports data from csv file
/// CSV file format:
/// time(int), name(string), working hours(int), oil tons(int), water tons(int)
/// comma separeted, without header
class CSVDataImporter
{
private:
  struct LineData
  {
    unsigned int time;
    std::string holeName;
    double  workHours;
    double  oilTons;
    double  waterTons;
  };
public:
  TaskData read(std::string filename);
};

#endif // DATAIMPORTER_H
