/*
 * ungapped_extension_with_trigger_seed_list.h
 *
 *  Created on: Sep 5, 2014
 *      Author: shu
 */

#ifndef UNGAPPED_EXTENSION_WITH_TRIGGER_SEED_LIST_H_
#define UNGAPPED_EXTENSION_WITH_TRIGGER_SEED_LIST_H_

#include <cstdlib>
#include <stdint.h>
#include <vector>
#include <iostream>
#include "host_seeds_memory.h"

class UngappedExtensionWithTriggerSeedList {
public:
	UngappedExtensionWithTriggerSeedList();
	virtual ~UngappedExtensionWithTriggerSeedList();
	void SetUp(size_t threhold_size, HostSeedsMemory* host_memory);

	HostSeedsMemory* MoveHostMemory() {
		assert(seeds_memory_ != NULL);
		HostSeedsMemory *seeds_memory = seeds_memory_;
		seeds_memory_query_seed_ids_ = NULL;
		seeds_memory_query_seed_starts_ = NULL;
		seeds_memory_database_positions_ = NULL;
		seeds_memory_ = NULL;
		return seeds_memory;
	}

	void ReleaseHostMemory() {
		assert(seeds_memory_ != NULL);
		seeds_memory_query_seed_ids_ = NULL;
		seeds_memory_query_seed_starts_ = NULL;
		seeds_memory_database_positions_ = NULL;
		seeds_memory_->Release();
		seeds_memory_ = NULL;
	}
	void Clear() {
		ReleaseHostMemory();
		ClearOverSeeds();
		length_ = 0;
		seeds_length_threhold_ = 0;
		query_seeds_length_ = 0;
	}
	size_t GetSeedsLength() {
		return length_;
	}

	size_t GetQuerySeedsLength() {
		return query_seeds_length_;
	}

	size_t GetSeedsLengthThreshold() {
		return seeds_length_threhold_;
	}

	bool Filled();
	void AddSeed(uint32_t query_seed_id, uint32_t database_position);
	void AddSeed(uint32_t query_seed_id, const uint32_t* database_positions,
			uint32_t database_positions_length);
	void ClearOverSeeds();

private:
	size_t length_;
	size_t seeds_length_threhold_;
	size_t query_seeds_length_;
	HostSeedsMemory *seeds_memory_;
	uint32_t *seeds_memory_query_seed_ids_;
	uint32_t *seeds_memory_query_seed_starts_;
	uint32_t *seeds_memory_database_positions_;
	std::vector<uint32_t> over_query_seed_ids_;
	std::vector<uint32_t> over_query_seed_starts_;
	std::vector<uint32_t> over_database_positions_;
};

inline bool UngappedExtensionWithTriggerSeedList::Filled() {
	return length_ >= seeds_length_threhold_;
}

inline void UngappedExtensionWithTriggerSeedList::AddSeed(
		uint32_t query_seed_id, uint32_t database_position) {
	if (length_ < seeds_length_threhold_) {
		seeds_memory_query_seed_ids_[query_seeds_length_] = query_seed_id;
		seeds_memory_query_seed_starts_[query_seeds_length_] = length_;
		seeds_memory_database_positions_[length_] = database_position;
		++query_seeds_length_;
	} else {
		over_query_seed_ids_.push_back(query_seed_id);
		over_query_seed_starts_.push_back(length_);
		over_database_positions_.push_back(database_position);
	}
	++length_;
}

inline void UngappedExtensionWithTriggerSeedList::AddSeed(
		uint32_t query_seed_id, const uint32_t* database_positions,
		uint32_t database_positions_length) {
	size_t buffer_set_length = 0;
	if (length_ < seeds_length_threhold_) {
		buffer_set_length = std::min(database_positions_length,
				(uint32_t) (seeds_length_threhold_ - length_));
		seeds_memory_query_seed_ids_[query_seeds_length_] = query_seed_id;
		seeds_memory_query_seed_starts_[query_seeds_length_] = length_;
		memcpy(&seeds_memory_database_positions_[length_], database_positions,
				sizeof(uint32_t) * buffer_set_length);
		++query_seeds_length_;
	}
	if (buffer_set_length != database_positions_length) {
		over_query_seed_ids_.push_back(query_seed_id);
		over_query_seed_starts_.push_back(length_ + buffer_set_length);
		for (uint32_t i = buffer_set_length; i < database_positions_length;
				++i) {
			over_database_positions_.push_back(database_positions[i]);
		}
	}
	length_ += database_positions_length;
}

inline void UngappedExtensionWithTriggerSeedList::ClearOverSeeds() {
	over_query_seed_ids_.clear();
	over_query_seed_starts_.clear();
	over_database_positions_.clear();
}

#endif /* UNGAPPED_EXTENSION_WITH_TRIGGER_SEED_LIST_H_ */
