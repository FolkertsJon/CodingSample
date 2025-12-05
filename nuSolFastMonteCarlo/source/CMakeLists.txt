# Cmake minimum version
cmake_minimum_required(VERSION 2.8)


# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)


# set the project name and version
project(nuSolPerformance)


# Locate the ROOT package and define a number of useful targets and variables.
find_package(ROOT REQUIRED COMPONENTS RIO Net RooFit RooFitCore Gui)
message(${ROOT_INCLUDE_DIRS})
# find and link Open MP for parallel
#find_package(OpenMP)
include_directories(${ROOT_INCLUDE_DIRS})
include_directories(${CMAKE_SOURCE_DIR}/../../mySource/include/)
include_directories(${CMAKE_SOURCE_DIR}/../../eigen-3.3.9/)
message(${CMAKE_SOURCE_DIR})


# Put sources and header together	
file(GLOB sources ${CMAKE_SOURCE_DIR}/../../mySource/src/*.cpp)
file(GLOB headers ${CMAKE_SOURCE_DIR}/../../mySource/include/*.hh)


# Create the nuSolPerformance program using the library.
add_executable(nuSolPerformance nuSolPerformance.cpp ${headers} ${sources})
target_link_libraries(nuSolPerformance PUBLIC ROOT::Physics ROOT::Gpad ROOT::Tree ROOT::Hist fftw3)
# Things I've tried:
#  ${LIBPATH} ROOT::Canvas ROOT::Pad ROOT::canvas




