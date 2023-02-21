#include "gtest/gtest.h"
#include "../src/opp_lib.h"

using namespace opp;

TEST(MeshSet, test_mesh_set_and_dat)
{
    std::shared_ptr<MeshSet> set_1 = std::make_shared<MeshSet>(std::string("set_1"), 10);
    std::shared_ptr<MeshSet> set_2 = std::make_shared<MeshSet>(std::string("set_2"), 20);

    std::shared_ptr<MeshDat<INT>> set_1_dat_1 = set_1->createDat<INT>("set_1_dat_1", 1);
    std::shared_ptr<MeshDat<REAL>> set_1_dat_2 = set_1->createDat<REAL>("set_1_dat_2", 2);
    std::shared_ptr<MeshDat<INT>> set_1_dat_3 = set_1->createDat<INT>("set_1_dat_3", 3);

    EXPECT_EQ(set_1->getName(), std::string("set_1"));	
    EXPECT_EQ(set_2->getName(), std::string("set_2"));	
    EXPECT_EQ(set_1->getSize(), 10);	
    EXPECT_EQ(set_2->getSize(), 20);	

    EXPECT_EQ(set_1_dat_1->getName(), std::string("set_1_dat_1"));	
    EXPECT_EQ(set_1_dat_2->getName(), std::string("set_1_dat_2"));	
    EXPECT_EQ(set_1_dat_3->getName(), std::string("set_1_dat_3"));	

    EXPECT_EQ(set_1_dat_1->getDimension(), 1);	
    EXPECT_EQ(set_1_dat_2->getDimension(), 2);	
    EXPECT_EQ(set_1_dat_3->getDimension(), 3);	

    EXPECT_EQ(set_1_dat_1->getSize(), 1 * sizeof(INT));	
    EXPECT_EQ(set_1_dat_2->getSize(), 2 * sizeof(REAL));	
    EXPECT_EQ(set_1_dat_3->getSize(), 3 * sizeof(INT));
}

TEST(MeshSet, test_mesh_set_and_dat_with_data)
{
    std::shared_ptr<MeshSet> set_1 = std::make_shared<MeshSet>(std::string("set_1"), 3);
    std::shared_ptr<MeshSet> set_2 = std::make_shared<MeshSet>(std::string("set_2"), 20);

    double test_double[6] = { 1.01, 2.02, 3.03, 4.04, 5.05, 6.06 };
    int test_int[10] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };

    std::shared_ptr<MeshDat<INT>> set_1_dat_1 = set_1->createDat<INT>("set_1_dat_1", 1);
    std::shared_ptr<MeshDat<REAL>> set_1_dat_2 = set_1->createDat<REAL>("set_1_dat_2", 2, test_double);
    std::shared_ptr<MeshDat<INT>> set_1_dat_3 = set_1->createDat<INT>("set_1_dat_3", 3, test_int);

    EXPECT_EQ(set_1->getName(), std::string("set_1"));	
    EXPECT_EQ(set_2->getName(), std::string("set_2"));	
    EXPECT_EQ(set_1->getSize(), 3);	
    EXPECT_EQ(set_2->getSize(), 20);	

    EXPECT_EQ(set_1_dat_1->getName(), std::string("set_1_dat_1"));	
    EXPECT_EQ(set_1_dat_2->getName(), std::string("set_1_dat_2"));	
    EXPECT_EQ(set_1_dat_3->getName(), std::string("set_1_dat_3"));	

    EXPECT_EQ(set_1_dat_1->getDimension(), 1);	
    EXPECT_EQ(set_1_dat_2->getDimension(), 2);	
    EXPECT_EQ(set_1_dat_3->getDimension(), 3);	

    EXPECT_EQ(set_1_dat_1->getSize(), 1 * sizeof(INT));	
    EXPECT_EQ(set_1_dat_2->getSize(), 2 * sizeof(REAL));	
    EXPECT_EQ(set_1_dat_3->getSize(), 3 * sizeof(INT));

    std::vector<double>& double_vec = set_1_dat_2->getData();
    std::vector<int>& int_vec = set_1_dat_3->getData();

    EXPECT_EQ(double_vec.size(), 6);
    EXPECT_EQ(int_vec.size(), 9);

    for (int i = 0; i < 6; i++)
        EXPECT_EQ(test_double[i], double_vec[i]);
    
    for (int i = 0; i < 9; i++)
        EXPECT_EQ(test_int[i], int_vec[i]);
}

TEST(ParticleSet, test_particle_set_and_dat)
{
    std::shared_ptr<MeshSet> set_1 = std::make_shared<MeshSet>(std::string("set_1"), 10);
    std::shared_ptr<ParticleSet> set_2 = std::make_shared<ParticleSet>(std::string("set_2"), set_1, 20);
    std::shared_ptr<ParticleSet> set_3 = std::make_shared<ParticleSet>(std::string("set_3"), set_1);

    std::shared_ptr<ParticleDat<INT>> set_2_dat_1 = set_2->createDat<INT>("set_2_dat_1", 1);
    std::shared_ptr<ParticleDat<REAL>> set_2_dat_2 = set_2->createDat<REAL>("set_2_dat_2", 2);
    std::shared_ptr<ParticleDat<INT>> set_2_dat_3 = set_2->createDat<INT>("set_2_dat_3", 3);

    std::shared_ptr<ParticleDat<INT>> set_3_dat_1 = set_3->createDat<INT>("set_3_dat_1", 1);
    std::shared_ptr<ParticleDat<REAL>> set_3_dat_2 = set_3->createDat<REAL>("set_3_dat_2", 2);
    std::shared_ptr<ParticleDat<INT>> set_3_dat_3 = set_3->createDat<INT>("set_3_dat_3", 3);

    EXPECT_EQ(set_1->getName(), std::string("set_1"));	
    EXPECT_EQ(set_2->getName(), std::string("set_2"));
    EXPECT_EQ(set_3->getName(), std::string("set_3"));	
    
    EXPECT_EQ(set_1->getSize(), 10);	
    EXPECT_EQ(set_2->getSize(), 20);	
    EXPECT_EQ(set_3->getSize(), 0);

    EXPECT_EQ(set_2->getCellsSet()->getName(), std::string("set_1"));	
    EXPECT_EQ(set_3->getCellsSet()->getName(), std::string("set_1"));

    EXPECT_EQ(set_2_dat_1->getName(), std::string("set_2_dat_1"));	
    EXPECT_EQ(set_2_dat_2->getName(), std::string("set_2_dat_2"));	
    EXPECT_EQ(set_2_dat_3->getName(), std::string("set_2_dat_3"));	

    EXPECT_EQ(set_3_dat_1->getName(), std::string("set_3_dat_1"));	
    EXPECT_EQ(set_3_dat_2->getName(), std::string("set_3_dat_2"));	
    EXPECT_EQ(set_3_dat_3->getName(), std::string("set_3_dat_3"));	

    EXPECT_EQ(set_2_dat_1->getDimension(), 1);	
    EXPECT_EQ(set_2_dat_2->getDimension(), 2);	
    EXPECT_EQ(set_2_dat_3->getDimension(), 3);	

    EXPECT_EQ(set_3_dat_1->getDimension(), 1);	
    EXPECT_EQ(set_3_dat_2->getDimension(), 2);	
    EXPECT_EQ(set_3_dat_3->getDimension(), 3);	

    EXPECT_EQ(set_2_dat_1->getSize(), 1 * sizeof(INT));	
    EXPECT_EQ(set_2_dat_2->getSize(), 2 * sizeof(REAL));	
    EXPECT_EQ(set_2_dat_3->getSize(), 3 * sizeof(INT));

    EXPECT_EQ(set_3_dat_1->getSize(), 1 * sizeof(INT));	
    EXPECT_EQ(set_3_dat_2->getSize(), 2 * sizeof(REAL));	
    EXPECT_EQ(set_3_dat_3->getSize(), 3 * sizeof(INT));
}

TEST(ParticleSet, test_particle_set_and_dat_with_data)
{
    std::shared_ptr<MeshSet> set_0 = std::make_shared<MeshSet>(std::string("set_0"), 10);
    std::shared_ptr<ParticleSet> set_1 = std::make_shared<ParticleSet>(std::string("set_1"), set_0, 3);

    EXPECT_EQ(set_1->getCellsSet()->getName(), std::string("set_0"));

    double test_double[6] = { 1.01, 2.02, 3.03, 4.04, 5.05, 6.06 };
    int test_int[10] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
    int test_cell_idx[4] = { 10, 9, 8, 7 };

    std::shared_ptr<ParticleDat<INT>> set_1_dat_1 = set_1->createDat<INT>("set_1_dat_1", 1, test_cell_idx, OPP_DAT_CELL_INDEX);
    std::shared_ptr<ParticleDat<REAL>> set_1_dat_2 = set_1->createDat<REAL>("set_1_dat_2", 2, test_double);
    std::shared_ptr<ParticleDat<INT>> set_1_dat_3 = set_1->createDat<INT>("set_1_dat_3", 3, test_int);

    EXPECT_EQ(set_0->getName(), std::string("set_0"));	
    EXPECT_EQ(set_1->getName(), std::string("set_1"));	
    EXPECT_EQ(set_0->getSize(), 10);	
    EXPECT_EQ(set_1->getSize(), 3);	

    EXPECT_EQ(set_1_dat_1->getName(), std::string("set_1_dat_1"));	
    EXPECT_EQ(set_1_dat_2->getName(), std::string("set_1_dat_2"));	
    EXPECT_EQ(set_1_dat_3->getName(), std::string("set_1_dat_3"));	

    EXPECT_EQ(set_1_dat_1->getDimension(), 1);	
    EXPECT_EQ(set_1_dat_2->getDimension(), 2);	
    EXPECT_EQ(set_1_dat_3->getDimension(), 3);	

    EXPECT_EQ(set_1_dat_1->getSize(), 1 * sizeof(INT));	
    EXPECT_EQ(set_1_dat_2->getSize(), 2 * sizeof(REAL));	
    EXPECT_EQ(set_1_dat_3->getSize(), 3 * sizeof(INT));

    std::vector<int>& icell_idx_vec = set_1_dat_1->getData();
    std::vector<double>& double_vec = set_1_dat_2->getData();
    std::vector<int>& int_vec = set_1_dat_3->getData();

    EXPECT_EQ(icell_idx_vec.size(), 3);
    EXPECT_EQ(double_vec.size(), 6);
    EXPECT_EQ(int_vec.size(), 9);

    for (int i = 0; i < 3; i++)
        EXPECT_EQ(test_cell_idx[i], icell_idx_vec[i]);

    for (int i = 0; i < 6; i++)
        EXPECT_EQ(test_double[i], double_vec[i]);
    
    for (int i = 0; i < 9; i++)
        EXPECT_EQ(test_int[i], int_vec[i]);
}

TEST(MeshSet, test_mesh_set_and_map)
{
    std::shared_ptr<MeshSet> set_1 = std::make_shared<MeshSet>(std::string("set_1"), 2);
    std::shared_ptr<MeshSet> set_2 = std::make_shared<MeshSet>(std::string("set_2"), 3);
    std::shared_ptr<MeshSet> set_3 = std::make_shared<MeshSet>(std::string("set_3"), 4);

    int test_int1[10] = { 11, 12, 13, 14, 15, 16, 17, 18, 19 };
    int test_int2[10] = { 21, 22, 23, 24, 25, 26, 27, 28, 29 };
    int test_int3[10] = { 31, 32, 33, 34, 35, 36, 37, 38, 39 };

    std::shared_ptr<Map> set_1_map_2 = set_1->createMapTo(set_2, "map_1_2", 2, test_int1);
    std::shared_ptr<Map> set_1_map_3 = set_1->createMapTo(set_3, "map_1_3", 3, test_int2);
    std::shared_ptr<Map> set_2_map_3 = set_2->createMapTo(set_3, "map_2_3", 1, test_int3);

    EXPECT_EQ(set_1->getName(), std::string("set_1"));	
    EXPECT_EQ(set_2->getName(), std::string("set_2"));
    EXPECT_EQ(set_3->getName(), std::string("set_3"));	

    EXPECT_EQ(set_1->getSize(), 2);	
    EXPECT_EQ(set_2->getSize(), 3);	
    EXPECT_EQ(set_3->getSize(), 4);

    EXPECT_EQ(set_1_map_2->getName(), std::string("map_1_2"));	
    EXPECT_EQ(set_1_map_3->getName(), std::string("map_1_3"));	
    EXPECT_EQ(set_2_map_3->getName(), std::string("map_2_3"));	

    EXPECT_EQ(set_1_map_2->getFrom()->getName(), std::string("set_1"));	
    EXPECT_EQ(set_1_map_3->getFrom()->getName(), std::string("set_1"));	
    EXPECT_EQ(set_2_map_3->getFrom()->getName(), std::string("set_2"));	

    EXPECT_EQ(set_1_map_2->getTo()->getName(), std::string("set_2"));	
    EXPECT_EQ(set_1_map_3->getTo()->getName(), std::string("set_3"));	
    EXPECT_EQ(set_2_map_3->getTo()->getName(), std::string("set_3"));	

    EXPECT_EQ(set_1_map_2->getDimension(), 2);	
    EXPECT_EQ(set_1_map_3->getDimension(), 3);	
    EXPECT_EQ(set_2_map_3->getDimension(), 1);	

    const std::vector<INT>& int_vec1 = set_1_map_2->getMapData();
    const std::vector<INT>& int_vec2 = set_1_map_3->getMapData();
    const std::vector<INT>& int_vec3 = set_2_map_3->getMapData();

    EXPECT_EQ(int_vec1.size(), 4);
    EXPECT_EQ(int_vec2.size(), 6);
    EXPECT_EQ(int_vec3.size(), 3);

    for (int i = 0; i < 4; i++)
        EXPECT_EQ(test_int1[i], int_vec1[i]);
    
    for (int i = 0; i < 6; i++)
        EXPECT_EQ(test_int2[i], int_vec2[i]);
    
    for (int i = 0; i < 3; i++)
        EXPECT_EQ(test_int3[i], int_vec3[i]);

    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 3; j++)
        EXPECT_EQ(test_int2[i * 3 + j], set_1_map_3->get(i, j));
    
    for (int i = 0; i < 3; i++)
        EXPECT_EQ(test_int3[i], set_2_map_3->get(i));
}

TEST(MeshSet, test_mesh_set_and_map_covert_to_zero_indexing)
{
    std::shared_ptr<MeshSet> set_1 = std::make_shared<MeshSet>(std::string("set_1"), 2);
    std::shared_ptr<MeshSet> set_2 = std::make_shared<MeshSet>(std::string("set_2"), 3);

    int test_int1[10] = { 11, 12, 13, 14, 15, 16, 17, 18, 19 };

    SIM->OPP_maps_base_index = 1;

    std::shared_ptr<Map> set_1_map_2 = set_1->createMapTo(set_2, "map_1_2", 2, test_int1);

    EXPECT_EQ(set_1->getName(), std::string("set_1"));	

    EXPECT_EQ(set_1->getSize(), 2);	
    EXPECT_EQ(set_2->getSize(), 3);

    EXPECT_EQ(set_1_map_2->getName(), std::string("map_1_2"));	
    EXPECT_EQ(set_1_map_2->getFrom()->getName(), std::string("set_1"));	
    EXPECT_EQ(set_1_map_2->getTo()->getName(), std::string("set_2"));	
    EXPECT_EQ(set_1_map_2->getDimension(), 2);	

    const std::vector<INT>& int_vec1 = set_1_map_2->getMapData();

    EXPECT_EQ(int_vec1.size(), 4);

    for (int i = 0; i < 4; i++)
        EXPECT_EQ((test_int1[i] - 1), int_vec1[i]);
    
    SIM->releaseInstance();
}

TEST(MeshSet, test_mesh_set_dat_and_arg)
{
    std::shared_ptr<MeshSet> set_1 = std::make_shared<MeshSet>(std::string("set_1"), 2);
    std::shared_ptr<MeshSet> set_2 = std::make_shared<MeshSet>(std::string("set_2"), 3);

    INT test_int1[10] = { 11, 12, 13, 14, 15, 16, 17, 18, 19 };
    REAL test_doub2[10] = { 21.01, 22.01, 23.01, 24.01, 25.01, 26.01, 27.01, 28.01, 29.01 };
    INT test_int3[10] = { 31, 32, 33, 34, 35, 36, 37, 38, 39 };

    std::shared_ptr<Map> set_1_map_2 = set_1->createMapTo(set_2, "map_1_2", 2, test_int1);

    std::shared_ptr<MeshDat<REAL>> set_1_dat_1 = set_1->createDat<REAL>("set_1_dat_1", 2, test_doub2);
    std::shared_ptr<MeshDat<INT>> set_1_dat_2 = set_1->createDat<INT>("set_1_dat_2", 3, test_int3);

    Arg arg11 = set_1_dat_1->createArg(OPP_WRITE);
    Arg arg12 = set_1_dat_2->createArg(5, set_1_map_2, OPP_INC, true);
    Arg arg13 = set_1_dat_1->createArg(6, set_1_map_2, 4, std::string("double"), OPP_MAX, true);

    EXPECT_EQ(arg11.type, std::string("double"));	
    EXPECT_EQ(arg12.type, std::string("int"));
    EXPECT_EQ(arg13.type, std::string("double"));	

    EXPECT_EQ(arg11.dim, 2);	
    EXPECT_EQ(arg12.dim, 3);
    EXPECT_EQ(arg13.dim, 4);	

    EXPECT_EQ(arg11.idx, 0);	
    EXPECT_EQ(arg12.idx, 5);
    EXPECT_EQ(arg13.idx, 6);

    EXPECT_EQ(arg11.acc, OPP_WRITE);	
    EXPECT_EQ(arg12.acc, OPP_INC);
    EXPECT_EQ(arg13.acc, OPP_MAX);	

    EXPECT_EQ(arg11.argtype, OPP_ARG_DAT);	
    EXPECT_EQ(arg12.argtype, OPP_ARG_DAT);
    EXPECT_EQ(arg13.argtype, OPP_ARG_DAT);	

    EXPECT_EQ(arg11.map_with_cell_index, false);	
    EXPECT_EQ(arg12.map_with_cell_index, true);
    EXPECT_EQ(arg13.map_with_cell_index, true);	

    EXPECT_EQ((std::dynamic_pointer_cast<MeshDat<REAL>>(arg11.dat))->getName(), "set_1_dat_1");	
    EXPECT_EQ((std::dynamic_pointer_cast<MeshDat<INT>>(arg12.dat))->getName(), "set_1_dat_2");
    EXPECT_EQ((std::dynamic_pointer_cast<MeshDat<REAL>>(arg13.dat))->getName(), "set_1_dat_1");	

    EXPECT_EQ(arg11.map, nullptr);	
    EXPECT_EQ(arg12.map->getName(), "map_1_2");
    EXPECT_EQ(arg13.map->getName(), "map_1_2");

    std::vector<REAL>& double_vec1 = arg11.getRData();
    std::vector<INT>& int_vec1 = arg12.getIData();
    std::vector<REAL>& double_vec2 = arg13.getRData();

    EXPECT_EQ(double_vec1.size(), 4);
    EXPECT_EQ(int_vec1.size(), 6);
    EXPECT_EQ(double_vec2.size(), 4);

    for (int i = 0; i < 4; i++)
    {
        EXPECT_EQ(test_doub2[i], double_vec1[i]);
        EXPECT_EQ(test_doub2[i], double_vec2[i]);
    }
    for (int i = 0; i < 6; i++)
    {
        EXPECT_EQ((test_int3[i]), int_vec1[i]);
    }

    const std::vector<INT>& int_map1 = arg12.getMapData();
    const std::vector<INT>& int_map2 = arg13.getMapData();

    EXPECT_EQ(int_map1.size(), 4);
    EXPECT_EQ(int_map2.size(), 4);

    for (int i = 0; i < 4; i++)
    {
        EXPECT_EQ(test_int1[i], int_map1[i]);
        EXPECT_EQ(test_int1[i], int_map2[i]);
    }
}

TEST(ParticleSet, test_particle_set_dat_and_arg)
{
    std::shared_ptr<MeshSet> set_0 = std::make_shared<MeshSet>(std::string("set_0"), 2);
    std::shared_ptr<MeshSet> set_1 = std::make_shared<MeshSet>(std::string("set_1"), 3);
    std::shared_ptr<ParticleSet> set_2 = std::make_shared<ParticleSet>(std::string("set_2"), set_1, 3);

    INT test_int1[10] = { 11, 12, 13, 14, 15, 16, 17, 18, 19 };
    REAL test_doub2[10] = { 21.01, 22.01, 23.01, 24.01, 25.01, 26.01, 27.01, 28.01, 29.01 };
    INT test_int3[10] = { 31, 32, 33, 34, 35, 36, 37, 38, 39 };

    std::shared_ptr<Map> set_1_map_0 = set_1->createMapTo(set_0, "map_1_0", 2, test_int1);

    std::shared_ptr<ParticleDat<REAL>> set_2_dat_1 = set_2->createDat<REAL>("set_2_dat_1", 2, test_doub2);
    std::shared_ptr<ParticleDat<INT>> set_2_dat_2 = set_2->createDat<INT>("set_2_dat_2", 3, test_int3);

    Arg arg11 = set_2_dat_1->createArg(OPP_WRITE);
    Arg arg12 = set_2_dat_2->createArg(5, set_1_map_0, OPP_INC, true);                          // Incorrect mapping, only for testing
    Arg arg13 = set_2_dat_1->createArg(6, set_1_map_0, 4, std::string("double"), OPP_MAX, true);   // Incorrect mapping, only for testing

    EXPECT_EQ(arg11.type, std::string("double"));	
    EXPECT_EQ(arg12.type, std::string("int"));
    EXPECT_EQ(arg13.type, std::string("double"));	

    EXPECT_EQ(arg11.dim, 2);	
    EXPECT_EQ(arg12.dim, 3);
    EXPECT_EQ(arg13.dim, 4);	

    EXPECT_EQ(arg11.idx, 0);	
    EXPECT_EQ(arg12.idx, 5);
    EXPECT_EQ(arg13.idx, 6);

    EXPECT_EQ(arg11.acc, OPP_WRITE);	
    EXPECT_EQ(arg12.acc, OPP_INC);
    EXPECT_EQ(arg13.acc, OPP_MAX);	

    EXPECT_EQ(arg11.argtype, OPP_ARG_DAT);	
    EXPECT_EQ(arg12.argtype, OPP_ARG_DAT);
    EXPECT_EQ(arg13.argtype, OPP_ARG_DAT);	

    EXPECT_EQ(arg11.map_with_cell_index, false);	
    EXPECT_EQ(arg12.map_with_cell_index, true);
    EXPECT_EQ(arg13.map_with_cell_index, true);	

    EXPECT_EQ((std::dynamic_pointer_cast<ParticleDat<REAL>>(arg11.dat))->getName(), "set_2_dat_1");	
    EXPECT_EQ((std::dynamic_pointer_cast<ParticleDat<INT>>(arg12.dat))->getName(), "set_2_dat_2");
    EXPECT_EQ((std::dynamic_pointer_cast<ParticleDat<REAL>>(arg13.dat))->getName(), "set_2_dat_1");	

    EXPECT_EQ(arg11.map, nullptr);	
    EXPECT_EQ(arg12.map->getName(), "map_1_0");
    EXPECT_EQ(arg13.map->getName(), "map_1_0");
}

TEST(MeshSet, test_mesh_set_map_and_arg)
{
    std::shared_ptr<MeshSet> set_1 = std::make_shared<MeshSet>(std::string("set_1"), 4);
    std::shared_ptr<MeshSet> set_2 = std::make_shared<MeshSet>(std::string("set_2"), 3);

    INT test_int1[10] = { 11, 12, 13, 14, 15, 16, 17, 18, 19 };

    std::shared_ptr<Map> set_1_map_2 = set_1->createMapTo(set_2, "map_1_2", 2, test_int1);

    Arg arg11 = set_1_map_2->createArg(OPP_WRITE);
    Arg arg12 = set_1_map_2->createArg(OPP_INC, true);

    EXPECT_EQ(arg11.type, std::string("int"));	
    EXPECT_EQ(arg12.type, std::string("int"));

    EXPECT_EQ(arg11.dim, 2);	
    EXPECT_EQ(arg12.dim, 2);

    EXPECT_EQ(arg11.idx, 0);	
    EXPECT_EQ(arg12.idx, 0);

    EXPECT_EQ(arg11.acc, OPP_WRITE);	
    EXPECT_EQ(arg12.acc, OPP_INC);

    EXPECT_EQ(arg11.argtype, OPP_ARG_MAP);	
    EXPECT_EQ(arg12.argtype, OPP_ARG_MAP);

    EXPECT_EQ(arg11.map_with_cell_index, false);	
    EXPECT_EQ(arg12.map_with_cell_index, true);

    EXPECT_EQ(arg11.dat, nullptr);	
    EXPECT_EQ(arg12.dat, nullptr);

    EXPECT_EQ(arg11.map->getName(), "map_1_2");	
    EXPECT_EQ(arg12.map->getName(), "map_1_2");

    const std::vector<INT>& int_map1 = arg11.getMapData();
    const std::vector<INT>& int_map2 = arg12.getMapData();

    EXPECT_EQ(int_map1.size(), 8);
    EXPECT_EQ(int_map2.size(), 8);

    for (int i = 0; i < 8; i++)
    {
        EXPECT_EQ(test_int1[i], int_map1[i]);
        EXPECT_EQ(test_int1[i], int_map2[i]);
    }
}



// TEST(MeshSet, XXXXX)
// {
//     const int size = 200000000;

//     std::vector<int> vec;
//     vec.reserve(size);

//     for (int i = 0; i < size; i++)
//         vec.push_back(i);

//     std::vector<int>::iterator it = vec.begin();

//     EXPECT_EQ(vec[995], it[995]);

//     for (int i = 0; i < size; i+=1)
//     {
//         EXPECT_EQ(vec[i], it[i]);
//     }
// }
