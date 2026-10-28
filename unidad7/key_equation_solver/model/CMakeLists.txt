cmake_minimum_required(VERSION 3.20)
project(kes_model LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Header-only core: GF(2^m) arithmetic and the three solvers.
add_library(kes_core INTERFACE)
target_include_directories(kes_core INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Python module kesmodel.
find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
find_package(pybind11 CONFIG REQUIRED)
pybind11_add_module(kesmodel bind/kesmodel_py.cpp)
target_link_libraries(kesmodel PRIVATE kes_core)
target_compile_options(kesmodel PRIVATE -Wall -Wextra -Wpedantic)

# Python module kessc: the SystemC model, only when SystemC is available.
find_package(SystemCLanguage CONFIG QUIET)
if(SystemCLanguage_FOUND)
    pybind11_add_module(kessc bind/kessc_py.cpp)
    target_link_libraries(kessc PRIVATE kes_core SystemC::systemc)
    target_compile_options(kessc PRIVATE -Wall -Wextra -Wpedantic)
else()
    message(STATUS "SystemC not found: kessc is not built")
endif()
