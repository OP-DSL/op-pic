#include "gtest/gtest.h"
#include <filesystem>

#include "../src/opp_params.h"

using namespace opp;

TEST(Params, test_params)
{
    std::string file_path = std::__fs::filesystem::current_path().string() + "/test/system.param";

    Params p(file_path);
    
    p.write(std::cout);

    EXPECT_EQ(p.get<STRING>("inlet_mesh"), "/ext-home/zl/phd/mini-pic/examples/one_stream/coarse/inlet.dat");
    EXPECT_EQ(p.get<REAL>("dt"), REAL(1e-12));
    EXPECT_EQ(p.get<INT>("max_iter"), INT(250));
    EXPECT_EQ(p.get<BOOL>("invert_normals"), BOOL(true));

    EXPECT_EQ(p.get<STRING>("AAA"), "NOT FOUND");
    EXPECT_EQ(p.get<REAL>("BBB"), DBL_MAX);
    EXPECT_EQ(p.get<INT>("CCC"), INT_MAX);
    EXPECT_EQ(p.get<BOOL>("DDD"), BOOL(false));

    EXPECT_EQ(p.get<float>("EEE"), NULL);

    p.add<STRING>("inlet_mesh", std::string("/ext-home/"));
    p.add<REAL>("dt", REAL(1e+11));
    p.add<INT>("max_iter", 300);
    p.add<BOOL>("invert_normals", false);
    p.add<INT>("new", 500);

    EXPECT_EQ(p.get<STRING>("inlet_mesh"), "/ext-home/");
    EXPECT_EQ(p.get<REAL>("dt"), REAL(1e+11));
    EXPECT_EQ(p.get<INT>("max_iter"), INT(300));
    EXPECT_EQ(p.get<BOOL>("invert_normals"), BOOL(false));
    EXPECT_EQ(p.get<INT>("new"), INT(500));
}