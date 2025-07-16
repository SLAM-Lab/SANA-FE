include(FetchContent)

FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# add unit_tests directory to the gtest suite
# MAKE SURE TO SET UP CMakeLists.txt IN ANY DIRECTORY YOU WANT TO ADD TO THE TEST SUITE
# you can follow the template in tests/unit_tests/CMakeLists.txt
add_subdirectory(unit_tests)

# # don't worry about this
# file(GLOB TEST_SOURCES
#     "${CMAKE_CURRENT_SOURCE_DIR}/test_*.cpp"
# )

# add_executable(sanafe_tests ${TEST_SOURCES})

# target_link_libraries(
#     sanafe_tests
#     GTest::gtest_main
# )

# target_link_libraries(
#     sanafe_tests
#     sanafe
# )

# include(GoogleTest)
# gtest_discover_tests(sanafe_tests)