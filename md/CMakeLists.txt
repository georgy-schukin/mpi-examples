cmake_minimum_required(VERSION 3.9)

project(md CXX)

add_executable(md
        cell.cpp
        cell_block.cpp
        main.cpp
        md.cpp
        neighbor.cpp
        output.cpp
        particle.cpp
        particle_ops.cpp
        slicer.cpp
        vector3.cpp)

find_package(MPI)

target_link_libraries(md PUBLIC MPI::MPI_CXX)

