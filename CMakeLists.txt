cmake_minimum_required(VERSION 2.6)
project(gptesttask)

add_executable(gptesttask main.cpp)

install(TARGETS gptesttask RUNTIME DESTINATION bin)
