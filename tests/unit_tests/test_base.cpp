#include <gtest/gtest.h>

TEST(InitialTest, CheckTestFunctionality) {
    // Making sure that the test framework is working correctly
    EXPECT_EQ(0, 0);
}

/*
    FOR WRITING NEW TESTS WITH GTEST PLEASE:
    1. follow the naming convention test_<test name>.cpp to avoid linking issues
    2. put any helper functions in the anonymous namespace at the top of the file
    3. any files needed for the test can be put in the tests/data directory (to reference them in your tests use the relative path from the subdirectory of the test in the build directory,
       e.g. "../../../tests/data/arch/test1.yaml")
*/