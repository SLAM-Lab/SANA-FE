#include <gtest/gtest.h>

// A mostly hacky set of unit tests for my experimental changes. These could
//  evolve to be something more permanently useful in the future though...
//  Once I restructure the scheduler / NoC code at least

// TODO: remove this too
#define ENABLE_DEBUG_PRINTS 1

#include "print.hpp"
#include "schedule.hpp"

namespace
{

class TestNewNocInfo : public ::testing::Test
{
protected:
    class NewNocInfoExposed : public sanafe::NewNocInfo
    {
    public:
        // Constructor that forwards to parent
        NewNocInfoExposed(const sanafe::Scheduler &scheduler)
                : sanafe::NewNocInfo(scheduler)
        {
        }
    };

    sanafe::Scheduler scheduler;
    std::unique_ptr<NewNocInfoExposed> noc;

    TestNewNocInfo()
    {
        // Initialize scheduler properties
        scheduler.noc_width_in_tiles = 8;
        scheduler.noc_height_in_tiles = 4;
        scheduler.max_cores_per_tile = 4;

        // Create noc after scheduler is initialized
        noc = std::make_unique<NewNocInfoExposed>(scheduler);
    }
};

TEST_F(TestNewNocInfo, CheckXY)
{
    sanafe::Flow flow1(0, 4, 0);
    auto [x, y] = sanafe::get_path_xy_increments(flow1, *noc);
    EXPECT_EQ(x, 0);
    EXPECT_EQ(y, 1);

    sanafe::Flow flow2(4, 0, 0);
    std::tie(x, y) = sanafe::get_path_xy_increments(flow2, *noc);
    EXPECT_EQ(x, 0);
    EXPECT_EQ(y, -1);

    sanafe::Flow flow3(0, 16, 0);
    std::tie(x, y) = sanafe::get_path_xy_increments(flow3, *noc);
    EXPECT_EQ(x, 1);
    EXPECT_EQ(y, 0);

    sanafe::Flow flow4(16, 0, 0);
    std::tie(x, y) = sanafe::get_path_xy_increments(flow4, *noc);
    EXPECT_EQ(x, -1);
    EXPECT_EQ(y, 0);
}

// TODO: very rough and ready tests
TEST_F(TestNewNocInfo, GetPath)
{
    sanafe::Flow flow1(0, 15, 0);

    auto path1 = flow1.generate_path(*noc);
    INFO("path length:%zu\n", path1.size());
    EXPECT_EQ(path1.size(), 5);
    // Sending core has different link, but we don't model the receiving cores
    EXPECT_EQ(path1[0], noc->idx(0, 0, sanafe::ndirections));
    EXPECT_EQ(path1[1], noc->idx(0, 1, 0));
    EXPECT_EQ(path1[2], noc->idx(0, 2, 0));
    EXPECT_EQ(path1[3], noc->idx(0, 3, 0));
    EXPECT_EQ(path1[4], noc->idx(0, 3, sanafe::ndirections + noc->max_cores_per_tile + 3));

    sanafe::Flow flow2(1, 15, 0);
    auto path2 = flow2.generate_path(*noc);
    EXPECT_EQ(path2.size(), 5);
    EXPECT_EQ(path2[0], noc->idx(0, 0, sanafe::ndirections + 1));
    EXPECT_EQ(path2[1], noc->idx(0, 1, 0));
    EXPECT_EQ(path2[2], noc->idx(0, 2, 0));
    EXPECT_EQ(path2[3], noc->idx(0, 3, 0));
    EXPECT_EQ(path2[4], noc->idx(0, 3, sanafe::ndirections + noc->max_cores_per_tile + 3));
}

TEST_F(TestNewNocInfo, OverlappingFlows)
{
    sanafe::Flow flow1(0, 15, 0);
    sanafe::Flow flow2(1, 15, 0);
    auto [first, count] = noc->find_overlapping_links(flow1, flow2);
    INFO("first:%zu count:%zu\n", first.value(), count);
    EXPECT_EQ(first, 1);
    EXPECT_EQ(count, 4);

    sanafe::Flow flow3(0, 12, 0); // X=0, Y=0->3
    sanafe::Flow flow4(4, 8, 0); // X=0, Y=1->2
    std::tie(first, count) = noc->find_overlapping_links(flow3, flow4);
    INFO("first:%zu count:%zu\n", first.value(), count);
    EXPECT_EQ(first, 2);
    EXPECT_EQ(count, 1);
    // Check the reverse applies too
    std::tie(first, count) = noc->find_overlapping_links(flow4, flow3);
    EXPECT_EQ(first, 1);
    EXPECT_EQ(count, 1);
}

}
