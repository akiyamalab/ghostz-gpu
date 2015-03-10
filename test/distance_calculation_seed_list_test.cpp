#include <gtest/gtest.h>
#include <string>
#include <stdint.h>
#include "../src/distance_calculation_seed_list.h"
#include "../src/host_seeds_memory.h"

using namespace std;

TEST(DistanceCalculationSeedList, SetUpAndAdd) {
	HostSeedsMemory host_memory;
	size_t host_memory_size = 2;
	host_memory.Init(0, host_memory_size);
	DistanceCalculationSeedList seed_list;
	seed_list.SetUp(host_memory_size, &host_memory);
	seed_list.StartAddingSeed(10);
	seed_list.AddSeed(0, 0);
	seed_list.AddSeed(1, 1);
	seed_list.FinishedAddingSeed();
	EXPECT_EQ(1, seed_list.GetNumberHashPositionDataListIds());
	EXPECT_EQ(10, seed_list.GetHashPositionDataListIds()[0]);

	EXPECT_EQ(0, host_memory.GetQueryConcatenatedPositions()[0]);
	EXPECT_EQ(1, host_memory.GetQueryConcatenatedPositions()[1]);

	host_memory.Release();
	seed_list.SetUp(host_memory_size, &host_memory);
	EXPECT_EQ(0, seed_list.GetNumberHashPositionDataListIds());
}

TEST(DistanceCalculationSeedList, SetUpAndAddOver) {
	HostSeedsMemory host_memory;
	size_t host_memory_size = 3;
	host_memory.Init(0, host_memory_size);
	DistanceCalculationSeedList seed_list;
	seed_list.SetUp(host_memory_size, &host_memory);
	seed_list.StartAddingSeed(10);
	seed_list.AddSeed(0, 0);
	seed_list.AddSeed(1, 1);
	seed_list.FinishedAddingSeed();
	seed_list.StartAddingSeed(11);
	seed_list.AddSeed(2, 2);
	seed_list.AddSeed(3, 3);
	seed_list.FinishedAddingSeed();

	EXPECT_EQ(2, seed_list.GetNumberHashPositionDataListIds());
	EXPECT_EQ(10, seed_list.GetHashPositionDataListIds()[0]);
	EXPECT_EQ(11, seed_list.GetHashPositionDataListIds()[1]);
	EXPECT_EQ(0, seed_list.GetNumberCalculatedSeeds());
	EXPECT_EQ(0, seed_list.GetSeedOffset(0));
	EXPECT_EQ(2, seed_list.GetSeedOffset(1));
	EXPECT_EQ(4, seed_list.GetSeedOffset(2));
	EXPECT_EQ(0, host_memory.GetQueryConcatenatedPositions()[0]);
	EXPECT_EQ(1, host_memory.GetQueryConcatenatedPositions()[1]);
	EXPECT_EQ(2, host_memory.GetQueryConcatenatedPositions()[2]);

	host_memory.Release();
	seed_list.SetUp(host_memory_size, &host_memory);
	EXPECT_EQ(1, seed_list.GetNumberHashPositionDataListIds());
	EXPECT_EQ(11, seed_list.GetHashPositionDataListIds()[0]);
	EXPECT_EQ(1, seed_list.GetNumberCalculatedSeeds());
	EXPECT_EQ(3, host_memory.GetQueryConcatenatedPositions()[0]);
}
