file(GLOB TEST_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/test_*.cpp"
)

add_executable(unit_tests ${TEST_SOURCES})

target_link_libraries(
    unit_tests
    GTest::gtest_main
)

# links the unit tests with the main sanafe library
target_link_libraries(
    unit_tests
    sanafe
)

include(GoogleTest)
gtest_discover_tests(unit_tests PROPERTIES LABELS "unit_tests")