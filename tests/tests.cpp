#include <gtest/gtest.h>
#include "chip.hpp"


TEST(InitialTest, CheckTestFunctionality) {
    // making sure that the test framework is working correctly
    EXPECT_EQ(0, 0);
}

TEST(InitialTest, CheckChipInitialization) {
    char *a = new char[100];
    (void)a;
}