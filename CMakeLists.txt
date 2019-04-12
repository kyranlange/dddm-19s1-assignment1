cmake_minimum_required(VERSION 3.13)
project(assignment1)

set(CMAKE_CXX_STANDARD 14)

#add_executable(assignment1 AAGenerator.cpp AAGenerator.h CAGenerator.cpp CAGenerator.h)
add_executable(AAGenerator AAGenerator.cpp)
add_executable(CAGenerator CAGenerator.cpp)