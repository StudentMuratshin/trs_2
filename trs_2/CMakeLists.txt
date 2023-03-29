# CMakeList.txt : CMake project for trs_2, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

# Add source to this project's executable.
add_executable (trs_2.1 "task_1.cpp")
add_executable (trs_2.2 "task_2.cpp")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET trs_2.1 PROPERTY CXX_STANDARD 20)
  set_property(TARGET trs_2.2 PROPERTY CXX_STANDARD 20)
endif()

# TODO: Add tests and install targets if needed.
