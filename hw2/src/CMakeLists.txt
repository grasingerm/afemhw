cmake_minimum_required(VERSION 2.8)

project(afemhw)

# compiler flags
set(CMAKE_BUILD_TYPE debug)
set(CMAKE_C_FLAGS_DEBUG "-Wall -Wextra -g -std=c99 -pedantic -pedantic-errors -O1")
set(CMAKE_C_FLAGS_RELEASE "-Wall -Wextra -g -std=c99 -pedantic -pedantic-errors -O3")
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -g -std=c++11 -pedantic -pedantic-errors -O1")
set(CMAKE_CXX_FLAGS_RELEASE "-Wall -Wextra -g -std=c++11 -pedantic -pedantic-errors -O3")

# build directories
set(CMAKE_BINARY_DIR build)
set(EXECUTABLE_OUTPUT_PATH bin)

if(UNIX)
    link_libraries(hw2 m)
endif(UNIX)

# dependencies
add_executable(hw2 main.cpp element.cpp)

target_link_libraries(hw2 armadillo)

install(TARGETS hw2 DESTINATION bin)

