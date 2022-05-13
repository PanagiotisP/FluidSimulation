cmake_minimum_required(VERSION 3.19)
# set the project name
project(FluidSimulation)

set(MSVC_VERSION 1920)
set(MSVC true)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(SRC_FILES src/FluidSimulator.cpp src/FluidDomain.cpp src/LevelSet.cpp src/MacGrid.cpp src/ParticleSet.cpp src/Renderer.cpp main.cpp)
# set(SRC_FILES main.cpp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /O2 /openmp" )
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /O2 /openmp" )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj" )

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
# set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /O2" )
# set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} /O2" )
set(OPENVDB_INSTALLATION_DIR FILEPATH C:/Users/ppap/Documents/Projects/openvdb-install-debug)
elseif (CMAKE_BUILD_TYPE STREQUAL "Release")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2" )
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /O2" )
set(OPENVDB_INSTALLATION_DIR FILEPATH C:/Users/ppap/Documents/Projects/openvdb-install-release)
STRING (REGEX REPLACE "/RTC(su|[1su])" "" CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
endif()

message(STATUS "MY_VARIABLE=${CONFIGURATION}")
# message(STATUS "${CMAKE_BUILD_TYPE}")



# add the executable
add_executable(main ${SRC_FILES})

# OpenVDB stuff
cmake_minimum_required(VERSION 3.15)
include(../vcpkg/scripts/buildsystems/vcpkg.cmake)
list(APPEND CMAKE_MODULE_PATH "..${OPENVDB_INSTALLATION_DIR}/lib/cmake/OpenVDB")
find_package(OpenVDB REQUIRED)
# end

include_directories(../sfml-install/include/ headers C:/Users/ppap/Documents/Projects/eigen-3.4.0 "${OPENVDB_INSTALLATION_DIR}/include")


target_link_directories(main PUBLIC ../sfml-install-release/lib ../sfml-install/lib C:/Users/ppap/Documents/Projects "${OPENVDB_INSTALLATION_DIR}/lib")
target_link_libraries(main debug sfml-window-d.lib debug sfml-system-d.lib debug sfml-graphics-d.lib OpenVDB::openvdb optimized sfml-window.lib optimized sfml-system.lib optimized sfml-graphics.lib OpenVDB::openvdb)

