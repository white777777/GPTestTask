cmake_minimum_required(VERSION 2.6)
project(gptesttask)

#find_package(EIGEN REQUIRED) 
#include_directories(${EIGEN_INCLUDE_DIR})

find_package (Threads REQUIRED)
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

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
tester.h
tester.cpp
analyze.h
analyze.cpp
)

target_link_libraries (gptesttask ${CMAKE_THREAD_LIBS_INIT})

install(TARGETS gptesttask RUNTIME DESTINATION bin)
