cmake_minimum_required(VERSION 2.8)

project(afemhw)

set(CMAKE_BUILD_TYPE debug)
set(CMAKE_C_FLAGS_DEBUG "-Wall -Wextra -g -std=c99 -pedantic -pedantic-errors")
set(CMAKE_C_FLAGS_RELEASE "-Wall -Wextra -g -std=c99 -pedantic -pedantic-errors -O3")
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -g -std=c++11 -pedantic -pedantic-errors")
set(CMAKE_CXX_FLAGS_RELEASE "-Wall -Wextra -g -std=c++11 -pedantic -pedantic-errors -O3")

# build directories
set(CMAKE_BINARY_DIR build)
set(EXECUTABLE_OUTPUT_PATH bin)

include_directories(../src)
include_directories(../../share)
link_directories(../src)

# dependencies
add_executable(test_element test_element.cpp ../src/element.cpp ../src/element.hpp)
add_executable(test_pt test_pt.cpp ../src/point.cpp ../src/point.hpp)

target_link_libraries(test_element armadillo lapack blas)

install(TARGETS test_element test_pt DESTINATION bin)
