cmake_minimum_required(VERSION 3.19)
# set the project name
project(FluidSimulation)

set(MSVC_VERSION 1920)
set(MSVC true)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(SRC_FILES src/FluidSimulator.cpp src/FluidDomain.cpp src/Grid.cpp src/LevelSet.cpp src/MacGrid.cpp src/ParticleSet.cpp src/Renderer.cpp main.cpp)

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /O2" )
# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /O2" )

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /O2" )
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} /O2" )


# set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2" )
# set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /O2" )
STRING (REGEX REPLACE "/RTC(su|[1su])" "" CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")

# add the executable
add_executable(main ${SRC_FILES})

include_directories(../sfml-install/include/ headers C:/Users/ppap/Documents/Projects/eigen-3.4.0)

target_link_directories(main PUBLIC ../sfml-install/lib)
target_link_libraries(main PUBLIC sfml-window-d.lib sfml-system-d.lib sfml-graphics-d.lib)

