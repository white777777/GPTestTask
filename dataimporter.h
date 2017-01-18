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
  TaskData read(std::string filename)
  {
    TaskData data;
    std::stringstream ss;
    std::string cell;
    std::map<std::string, size_t> holeNamesToIndex;
    
    std::ifstream f(filename, std::ifstream::in);
    while(!f.eof())
    {
      // TODO: escape stringstream overhead
      std::string lineStr;
      std::getline(f, lineStr);
      if(lineStr.size()<4)
        continue;
      std::stringstream lineStream(lineStr);
      LineData lineData{0};
      std::getline(lineStream, cell, ',');ss = std::stringstream(cell);
      ss>>lineData.time;
      std::getline(lineStream, lineData.holeName, ',');
      std::getline(lineStream, cell, ',');ss = std::stringstream(cell);
      ss>>lineData.workHours;
      std::getline(lineStream, cell, ',');ss = std::stringstream(cell);
      ss>>lineData.oilTons;
      std::getline(lineStream, cell, ',');ss = std::stringstream(cell);
      ss>>lineData.waterTons;

      if(holeNamesToIndex.find(lineData.holeName) == holeNamesToIndex.end())
      {
        holeNamesToIndex[lineData.holeName] = data.holes.size();
        HoleData holeData;
        holeData.name = lineData.holeName;
        data.holes.push_back(holeData);
      }
      size_t holeIndex = holeNamesToIndex[lineData.holeName];
      HoleData& workingHole = data.holes[holeIndex];
      workingHole.ts.push_back     (lineData.workHours);
      workingHole.qOils.push_back  (lineData.oilTons);
      workingHole.qWaters.push_back(lineData.waterTons);
    }
    return data;
  }
};

#endif // DATAIMPORTER_H
