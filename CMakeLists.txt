cmake_minimum_required(VERSION 2.6)
project(gptesttask)

#find_package(EIGEN REQUIRED) 
#include_directories(${EIGEN_INCLUDE_DIR})
include_directories("/usr/include/eigen3")

find_package (Threads REQUIRED)
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

find_package(Boost REQUIRED COMPONENTS filesystem iostreams program_options system serialization)
include_directories(${Boost_INCLUDE_DIR})

add_executable(gptesttask 
dataimporter.h
dataimporter.cpp
functions.h
solver.h
solver.cpp
taskdata.h
taskdata.cpp
regressionmodels.h
boost_serialization_eigen.h
main.cpp
tester.h
tester.cpp
analyze.h
analyze.cpp
)

target_link_libraries (gptesttask 
  ${CMAKE_THREAD_LIBS_INIT}
  ${Boost_LIBRARIES}
s)

install(TARGETS gptesttask RUNTIME DESTINATION bin)
