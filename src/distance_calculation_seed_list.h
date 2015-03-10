/*
 * distance_calculation_seed_list.h
 *
 *  Created on: Aug 27, 2014
 *      Author: shu
 */

#ifndef DISTANCE_CALCULATION_SEED_LIST_H_
#define DISTANCE_CALCULATION_SEED_LIST_H_

#include <stdint.h>
#include <limits.h>
#include <vector>
#include "distance_calculator_gpu.h"
#include "host_seeds_memory.h"

class DistanceCalculationSeedList {
public:
	DistanceCalculationSeedList();
	virtual ~DistanceCalculationSeedList();

	void SetUp(size_t threhold_size, HostSeedsMemory* buffer);
	HostSeedsMemory* MoveHostMemory() {
		HostSeedsMemory *seeds_memory = seeds_memory_;
		//seeds_memory_query_ids_ = NULL;
		seeds_memory_query_concatenated_positions_ = NULL;
		seeds_memory_database_positions_ = NULL;
		seeds_memory_ = NULL;
		return seeds_memory;
	}
	void ReleaseHostMemory() {
		//seeds_memory_query_ids_ = NULL;
		seeds_memory_query_concatenated_positions_ = NULL;
		seeds_memory_database_positions_ = NULL;
		seeds_memory_->Release();
		seeds_memory_ = NULL;
	}
	void Clear() {
		ReleaseHostMemory();
		size_ = 0;
		threhold_size_ = 0;
		number_calculated_seeds_ = 0;
	}
	bool FilledBuffer();
	void StartAddingSeed(uint32_t id);
	void FinishedAddingSeed();
	void AddSeed(uint32_t query_concatenated_position,
			uint32_t database_position);
	bool IsCompletedHashPositionDataList(uint32_t id);
	size_t GetNumberCalculatedSeeds() {
		return number_calculated_seeds_;
	}
	size_t GetNumberHashPositionDataListIds();
	size_t GetNumberSeeds() {
		return size_;
	}

	size_t GetNumberSeedsThreshold() {
		return threhold_size_;
	}
	size_t GetSeedOffset(uint32_t hash_position_data_list_id);
	uint32_t* GetHashPositionDataListIds();
	DistanceCalculatorGpu::Distance* GetDistances();
	size_t GetDistancesLength() {
		return distnaces_.size();
	}
#if 0
	std::vector<uint32_t> host_query_ids_;
	std::vector<uint32_t> host_query_concatenated_positions_;
	std::vector<uint32_t> host_database_positions_;
#endif

private:
	size_t size_;
	size_t threhold_size_;
	size_t number_calculated_seeds_;
	HostSeedsMemory *seeds_memory_;
	uint32_t *seeds_memory_query_concatenated_positions_;
	uint32_t *seeds_memory_database_positions_;
	std::vector<uint32_t> hash_position_data_list_ids_;
	std::vector<uint32_t> seed_offsets_;
	std::vector<uint32_t> over_query_concatenated_positions_;
	std::vector<uint32_t> over_database_positions_;
	std::vector<DistanceCalculatorGpu::Distance> distnaces_;
};

inline bool DistanceCalculationSeedList::FilledBuffer() {
	return size_ >= threhold_size_;
}

inline void DistanceCalculationSeedList::StartAddingSeed(uint32_t id) {
	assert(!FilledBuffer());
	hash_position_data_list_ids_.push_back(id);
	if (seed_offsets_.empty()) {
		seed_offsets_.push_back(0);
	}
	assert(hash_position_data_list_ids_.size() == seed_offsets_.size());
}

inline void DistanceCalculationSeedList::FinishedAddingSeed() {
	seed_offsets_.push_back(size_);
	assert(hash_position_data_list_ids_.size() + 1 == seed_offsets_.size());
}

inline void DistanceCalculationSeedList::AddSeed(
		uint32_t query_concatenated_position, uint32_t database_position) {
	if (size_ < threhold_size_) {
		seeds_memory_query_concatenated_positions_[size_] =
				query_concatenated_position;
		seeds_memory_database_positions_[size_] = database_position;
#if 0
		host_query_ids_.push_back(query_id);
		host_query_concatenated_positions_.push_back(
				query_concatenated_position);
		host_database_positions_.push_back(database_position);
#endif
	} else {
		over_query_concatenated_positions_.push_back(
				query_concatenated_position);
		over_database_positions_.push_back(database_position);
	}
	++size_;
}

inline bool DistanceCalculationSeedList::IsCompletedHashPositionDataList(
		uint32_t id) {
	return GetSeedOffset(id + 1) <= threhold_size_;
}

inline size_t DistanceCalculationSeedList::GetNumberHashPositionDataListIds() {
	return hash_position_data_list_ids_.size();
}

inline size_t DistanceCalculationSeedList::GetSeedOffset(
		uint32_t hash_position_data_list_ids_i) {
	assert(hash_position_data_list_ids_i < seed_offsets_.size());
	return seed_offsets_[hash_position_data_list_ids_i];
}

inline uint32_t* DistanceCalculationSeedList::GetHashPositionDataListIds() {
	return &hash_position_data_list_ids_[0];
}

inline DistanceCalculatorGpu::Distance* DistanceCalculationSeedList::GetDistances() {
	return &distnaces_[0];
}

#endif /* DISTANCE_CALCULATION_SEED_LIST_H_ */
