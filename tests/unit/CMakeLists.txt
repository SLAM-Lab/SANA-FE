file(GLOB TEST_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
)

add_executable(unit ${TEST_SOURCES})

target_link_libraries(
    unit
    GTest::gtest_main
)

# links the unit tests with the main sanafe library
target_link_libraries(
    unit
    sanafe
)

# Get the absolute path to the root project directory (SANA-FE)
get_filename_component(SANAFE_ROOT "${CMAKE_SOURCE_DIR}" ABSOLUTE)

target_compile_definitions(unit PRIVATE
    SANAFE_ROOT_PATH="${SANAFE_ROOT}"
)

include(GoogleTest)
gtest_discover_tests(unit
  PROPERTIES LABELS "dev;ci;unit"
)
