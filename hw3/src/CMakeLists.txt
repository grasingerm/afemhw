cmake_minimum_required(VERSION 2.8)

project(afemhw)

set(CC clang)
set(CXX clang)

set(CMAKE_BUILD_TYPE debug)
set(CMAKE_C_FLAGS_DEBUG "-Wall -Wextra -g -std=c99 -pedantic -pedantic-errors")
set(CMAKE_C_FLAGS_RELEASE "-Wall -Wextra -std=c99 -pedantic -pedantic-errors -O3")
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -g -std=c++11 -pedantic -pedantic-errors")
set(CMAKE_CXX_FLAGS_RELEASE "-Wall -Wextra -std=c++11 -pedantic -pedantic-errors -O3")

# build directories
set(CMAKE_BINARY_DIR build)
set(EXECUTABLE_OUTPUT_PATH bin)

# include and link directories
include_directories(../src)
include_directories(../../share)
link_directories(../src)

# dependencies
add_executable(hw3 main.cpp point.cpp point.hpp element.cpp element.hpp)

# link libraries
target_link_libraries(hw3 armadillo)

if (UNIX)
    target_link_libraries(hw3 m)
endif (UNIX)

install(TARGETS hw3 DESTINATION bin)
