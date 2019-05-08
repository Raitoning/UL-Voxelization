project(Ex_Display3D)
cmake_minimum_required(VERSION 2.8)
aux_source_directory(. SRC_LIST)

#Required in DGtal
FIND_PACKAGE(DGtal REQUIRED)
INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})

LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})

add_definitions(-g)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

ADD_EXECUTABLE(Voxelizer src/Voxelizer.cpp src/Algorithms.cpp)
ADD_EXECUTABLE(Tests src/Tests.cpp src/Algorithms.cpp)
ADD_EXECUTABLE(Viewer src/Viewer.cpp)

TARGET_LINK_LIBRARIES(Voxelizer ${DGTAL_LIBRARIES})
TARGET_LINK_LIBRARIES(Tests ${DGTAL_LIBRARIES})
TARGET_LINK_LIBRARIES(Viewer ${DGTAL_LIBRARIES})
