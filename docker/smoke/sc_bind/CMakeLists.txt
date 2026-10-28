cmake_minimum_required(VERSION 3.20)
project(sc_bind LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(SystemCLanguage CONFIG REQUIRED)
find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(sc_bind sc_bind.cpp)
target_link_libraries(sc_bind PRIVATE SystemC::systemc)
