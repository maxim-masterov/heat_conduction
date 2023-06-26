cmake_minimum_required(VERSION 3.0)
project(heat_conduction)

set(CMAKE_CXX_STANDARD 17)

add_executable(heat_conduction main.cpp equation.cpp equation.h grid.h)

target_link_libraries(${PROJECT_NAME} -lpetsc -lmpi)