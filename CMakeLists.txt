cmake_minimum_required(VERSION 3.26)
set(CMAKE_CXX_STANDARD 20)
project(NLA_MTH9821)

# Add your source files and headers
file(GLOB MATRIX_SOURCES "MatrixClass/*.hpp")
add_executable(NLA_MTH9821 main.cpp ${MATRIX_SOURCES})
#add_library(NLA_MTH9821 SHARED main.cpp include/MatrixClass)

# Set up include directories
target_include_directories(NLA_MTH9821 PUBLIC include)
