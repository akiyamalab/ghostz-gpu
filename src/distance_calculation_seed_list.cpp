/*
 * distance_calculation_seed_cpp
 *
 *  Created on: Aug 27, 2014
 *      Author: shu
 */

#include <iostream>
#include <limits.h>
#include <assert.h>
#include "distance_calculation_seed_list.h"

using namespace std;

DistanceCalculationSeedList::DistanceCalculationSeedList() :
		size_(0), threhold_size_(0), number_calculated_seeds_(0), seeds_memory_(
				NULL),
		//seeds_memory_query_ids_(NULL),
		seeds_memory_query_concatenated_positions_(NULL), seeds_memory_database_positions_(
				NULL) {

}

DistanceCalculationSeedList::~DistanceCalculationSeedList() {

}

void DistanceCalculationSeedList::SetUp(size_t threhold_size,
		HostSeedsMemory* buffer) {
	seeds_memory_ = buffer;
	seeds_memory_->Allocate();
	size_ = 0;
	size_t prev_threshold = threhold_size_;
	threhold_size_ = threhold_size;
	//seeds_memory_query_ids_ = seeds_memory_->GetQueryIds();
	seeds_memory_query_concatenated_positions_ =
			seeds_memory_->GetQueryConcatenatedPositions();
	seeds_memory_database_positions_ = seeds_memory_->GetDatabasePositions();
#if 0
	host_query_ids_.clear();
	host_query_concatenated_positions_.clear();
	host_database_positions_.clear();
#endif
	if (distnaces_.size() < threhold_size) {
		distnaces_.resize(threhold_size);
	}
	if (over_database_positions_.size() > 0) {
		size_t copied_size = std::min(threhold_size_,
				over_database_positions_.size());
		uint32_t last_hash_position_data_list_id =
				hash_position_data_list_ids_.back();
		uint32_t last_seed_offset =
				seed_offsets_[hash_position_data_list_ids_.size() - 1];
		size_t number_previous_calculated_seeds = prev_threshold
				- last_seed_offset;
		if (hash_position_data_list_ids_.size() == 1) {
			number_calculated_seeds_ = number_calculated_seeds_
					+ number_previous_calculated_seeds;
		} else {
			number_calculated_seeds_ = number_previous_calculated_seeds;
		}
		hash_position_data_list_ids_.clear();
		seed_offsets_.clear();
		StartAddingSeed(last_hash_position_data_list_id);
		size_ = over_database_positions_.size();
		memcpy(seeds_memory_query_concatenated_positions_,
				&over_query_concatenated_positions_[0],
				sizeof(over_query_concatenated_positions_[0]) * copied_size);
		memcpy(seeds_memory_database_positions_, &over_database_positions_[0],
				sizeof(over_database_positions_[0]) * copied_size);
#if 0
		for (size_t i = 0; i < copied_size; ++i) {
			host_query_ids_.push_back(over_query_ids_[i]);
			host_query_concatenated_positions_.push_back(
					over_query_concatenated_positions_[i]);
			host_database_positions_.push_back(over_database_positions_[i]);
		}
#endif
		FinishedAddingSeed();
		if (copied_size != over_database_positions_.size()) {
			size_t remain_size = over_database_positions_.size() - copied_size;
			for (size_t i = 0; i < remain_size; ++i) {
				//over_query_ids_[i] = over_query_ids_[copied_size + i];
				over_query_concatenated_positions_[i] =
						over_query_concatenated_positions_[copied_size + i];
				over_database_positions_[i] =
						over_database_positions_[copied_size + i];
			}
			//over_query_ids_.resize(remain_size);
			over_query_concatenated_positions_.resize(remain_size);
			over_database_positions_.resize(remain_size);
		} else {
			//over_query_ids_.clear();
			over_query_concatenated_positions_.clear();
			over_database_positions_.clear();
		}
	} else {
		number_calculated_seeds_ = 0;
		hash_position_data_list_ids_.clear();
		seed_offsets_.clear();
	}
}

