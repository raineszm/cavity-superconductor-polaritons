cmake_minimum_required(VERSION 3.8)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(bardasis_schrieffer)

# Dependencies

add_subdirectory(extern/pybind11)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(GSL REQUIRED)

add_subdirectory(src)
