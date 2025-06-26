include(FetchContent)

FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_executable(
    sanafe_tests
    tests.cpp
)

target_link_libraries(
    sanafe_tests
    GTest::gtest_main
)

target_link_libraries(
    sanafe_tests
    sanafe
)

include(GoogleTest)
gtest_discover_tests(sanafe_tests)