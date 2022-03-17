cmake_minimum_required(VERSION 3.19)
# set the project name
project(FluidSimulation)

set(MSVC_VERSION 1920)
set(MSVC true)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(SRC_FILES FluidSimulator.cpp Grid.cpp LevelSet.cpp MacGrid.cpp main.cpp)

# add the executable
add_executable(main ${SRC_FILES})

include_directories(../sfml-install/include/ headers C:/Users/ppap/Documents/Projects/eigen-3.4.0)

target_link_directories(main PUBLIC ../sfml-install/lib)
target_link_libraries(main PUBLIC sfml-window-d.lib sfml-system-d.lib sfml-graphics-d.lib)

