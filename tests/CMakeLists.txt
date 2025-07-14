cmake_minimum_required(VERSION 3.22)
project(sanafe_tests LANGUAGES CXX)

include(FetchContent)

FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

set(SANAFE_SRC
    ../src/models.cpp
    ../src/pipeline.cpp
    tests.cpp
)

add_executable(sanafe_tests ${SANAFE_SRC})

target_include_directories(sanafe_tests PRIVATE ../src)

target_link_libraries(sanafe_tests PRIVATE GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(sanafe_tests)