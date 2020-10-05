cmake_minimum_required(VERSION 3.10)

# set the project name and version
project(AlignmentActive)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package (Threads)

add_subdirectory(../AlignmentRepaC ../AlignmentRepaC_build)

include_directories(../rapidjson/include ../AlignmentC ../AlignmentRepaC)

add_library(AlignmentActive AlignmentActive.cpp)

add_executable(AlignmentActive_test main.cpp)

target_link_libraries(AlignmentActive_test PUBLIC AlignmentActive AlignmentRepaC AlignmentC Threads::Threads)

target_include_directories(AlignmentActive_test PUBLIC "${PROJECT_BINARY_DIR}")

install(TARGETS AlignmentActive_test DESTINATION lib)
install(FILES AlignmentActive.h DESTINATION include)
