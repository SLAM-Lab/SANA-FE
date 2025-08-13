#include <gtest/gtest.h>
#include <string>
#include <vector>
#include <stdexcept>
#include <filesystem>
#include "arg_parsing.hpp"

// Testing valid argument parsing
TEST(BasicInputTest, ParseValidInput) {
    std::vector<std::string> input = {"arch.yaml", "net.yaml", "100"};
    EXPECT_NO_THROW({
        RequiredProgramArgs args = parse_required_args(input, 0);
        EXPECT_EQ(args.arch_filename, "arch.yaml");
        EXPECT_EQ(args.network_filename, "net.yaml");
        EXPECT_EQ(args.timesteps_to_execute, 100);
    });
}

// Missing arguments
TEST(BasicInputTest, MissingArguments)
{
    std::vector<std::string> args = {"arch.yaml"};
    EXPECT_THROW(parse_required_args(args, 0), std::invalid_argument);
}

// Invalid timestep (non-numeric)
TEST(BasicInputTest, InvalidTimestepNonNumeric)
{
    std::vector<std::string> args = {"arch.yaml", "net.yaml", "abc"};
    EXPECT_THROW(parse_required_args(args, 0), std::invalid_argument);
}

// Invalid timestep (negative)
TEST(BasicInputTest, InvalidTimestepNegative)
{
    std::vector<std::string> args = {"arch.yaml", "net.yaml", "-10"};
    EXPECT_THROW(parse_required_args(args, 0), std::invalid_argument);
}

// Invalid timestep (zero)
TEST(BasicInputTest, InvalidTimestepZero)
{
    std::vector<std::string> args = {"arch.yaml", "net.yaml", "0"};
    EXPECT_THROW(parse_required_args(args, 0), std::invalid_argument);
}

// Trying to open a non-existent architecture file
TEST(BasicInputTest, FileDoesNotExist)
{
    std::vector<std::string> args = {"nonexistent_arch.yaml", "net.yaml", "100"};
    RequiredProgramArgs parsed_args = parse_required_args(args, 0);
    EXPECT_THROW(sanafe::load_arch(parsed_args.arch_filename), std::system_error);
}

TEST(BasicInputTest, ValidFile)
{
    std::filesystem::path path(SANAFE_ROOT_PATH);
    // std::cout << "Current path: " << path.string() << std::endl;
    std::vector<std::string> args = {path.string() + "/arch/example.yaml", path.string() + "/snn/example.yaml", "100"};
    RequiredProgramArgs parsed_args = parse_required_args(args, 0);
    EXPECT_NO_THROW(sanafe::load_arch(parsed_args.arch_filename));
}