/*
Example unit testing project from:
https://raymii.org/s/tutorials/Cpp_project_setup_with_cmake_and_unit_tests.html
*/

#include "gtest/gtest.h"


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}