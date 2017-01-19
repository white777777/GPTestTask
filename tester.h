#ifndef TESTER_H
#define TESTER_H
#include <iostream>
#include <cstdlib>
#include "solver.h"
#include "dataimporter.h"

class Tester
{
  void testExactSolution();
  void testExactSolutionPrint();
  void testSolver();
  void testRealWorldIterative();
  void testRealWorld();
public:
  void Test();
};


#endif // TESTER_H
