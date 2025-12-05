# Cmake minimum version
cmake_minimum_required(VERSION 2.8)


# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)


# set the project name and version
project(nuSolFastMonteCarlo)



# Locate the ROOT package and define a number of useful targets and variables.
find_package(ROOT REQUIRED COMPONENTS RIO Net)
message(${ROOT_INCLUDE_DIRS})
# find and link Open MP for parallel
find_package(OpenMP)
include_directories(${ROOT_INCLUDE_DIRS})
include_directories(${CMAKE_SOURCE_DIR}/../../mySource/include/)
include_directories(${CMAKE_SOURCE_DIR}/../../eigen-3.3.9/)


# Put sources and header together	
file(GLOB sources ${CMAKE_SOURCE_DIR}/../../mySource/src/*.cpp)
file(GLOB headers ${CMAKE_SOURCE_DIR}/../../mySource/include/*.hh)


# Create the nuSolFastMonteCarlo program using the library.
add_executable(nuSolFastMonteCarlo main.cpp ${headers} ${sources})
target_link_libraries(nuSolFastMonteCarlo PUBLIC OpenMP::OpenMP_CXX ROOT::Core ROOT::Physics ROOT::Tree ROOT::Gpad ROOT::Hist ROOT::RooFitCore fftw3)





