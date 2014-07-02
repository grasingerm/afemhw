cmake_minimum_required(VERSION 2.8)

project(shared_tests)

# compiler flags
set(CMAKE_BUILD_TYPE debug)
set(CMAKE_C_FLAGS_DEBUG "-Wall -Wextra -g -std=c99 -pedantic -pedantic-errors -O1")
set(CMAKE_C_FLAGS_RELEASE "-Wall -Wextra -g -std=c99 -pedantic -pedantic-errors -O3")
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -g -std=c++11 -pedantic -pedantic-errors -O1")
set(CMAKE_CXX_FLAGS_RELEASE "-Wall -Wextra -g -std=c++11 -pedantic -pedantic-errors -O3")

# build directories
set(CMAKE_BINARY_DIR build)
set(EXECUTABLE_OUTPUT_PATH bin)

# dependencies
add_executable(test_int test_int.cpp)

install(TARGETS test_int DESTINATION bin)

