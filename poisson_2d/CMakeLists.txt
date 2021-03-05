cmake_minimum_required(VERSION 3.9)

project(poisson_2d CXX)

add_executable(poisson_2d poisson_2d.cpp main.cpp)

find_package(MPI)

target_link_libraries(poisson_2d PUBLIC MPI::MPI_CXX)