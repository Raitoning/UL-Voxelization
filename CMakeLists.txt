project(Ex_Display3D)
cmake_minimum_required(VERSION 2.8)
aux_source_directory(. SRC_LIST)

#Required in DGtal
FIND_PACKAGE(DGtal REQUIRED)
INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})

add_definitions(-g)

LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})

ADD_EXECUTABLE(display3D src/display3D.cpp)
TARGET_LINK_LIBRARIES(display3D ${DGTAL_LIBRARIES})
