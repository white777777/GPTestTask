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
  void testRealWorld();
public:
  void Test()
  {
    testExactSolution();
    testSolver();
    testRealWorld();
  }
};


#endif // TESTER_H
