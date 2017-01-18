#ifndef TESTER_H
#define TESTER_H
#include <iostream>
#include <cstdlib>
#include "solver.h"
#include "dataimporter.h"

class Tester
{
  void testExactSolution();
  void testSolver();
public:
  void Test()
  {
    testExactSolution();
    testSolver();
  }
};


#endif // TESTER_H
