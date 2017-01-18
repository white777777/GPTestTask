cmake_minimum_required(VERSION 2.6)
project(gptesttask)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g3 -Wall")

add_executable(gptesttask 
main.cpp
dataimporter.h
dataimporter.cpp
functions.h
functions.cpp
solver.h
solver.cpp
taskdata.h
taskdata.cpp
regressionmodels.h
regressionmodels.cpp
tester.h
tester.cpp
)

install(TARGETS gptesttask RUNTIME DESTINATION bin)
