cmake_minimum_required(VERSION 3.26)
project(NLA_MTH9821)

# Add your source files and headers
add_library(NLA_MTH9821 SHARED src/main.cpp include/MatrixClass)

# Set up include directories
target_include_directories(NLA_MTH9821 PUBLIC include)
