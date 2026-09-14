cmake_minimum_required(VERSION 3.16)

project(nbody_simulator LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(nbody
    src/main.cpp
    src/gravity.cpp
    src/quadtree.cpp
    src/integrators.cpp
    src/scenarios.cpp
    src/io.cpp
    src/experiments.cpp
)

target_include_directories(nbody PRIVATE include)

add_executable(nbody_tests
    tests/tests.cpp
    src/gravity.cpp
    src/quadtree.cpp
    src/integrators.cpp
    src/scenarios.cpp
)

target_include_directories(nbody_tests PRIVATE include)

enable_testing()
add_test(NAME nbody_tests COMMAND nbody_tests)

if(MSVC)
    target_compile_options(nbody PRIVATE /W4)
    target_compile_options(nbody_tests PRIVATE /W4)
else()
    target_compile_options(nbody PRIVATE -Wall -Wextra -Wpedantic)
    target_compile_options(nbody_tests PRIVATE -Wall -Wextra -Wpedantic)
endif()
