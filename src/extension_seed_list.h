/*
 * extension_seed_list.h
 *
 *  Created on: Aug 27, 2014
 *      Author: shu
 */

#ifndef EXTENSION_SEED_LIST_H_
#define EXTENSION_SEED_LIST_H_

class Queries;

#include <cstdlib>
#include <stdint.h>
#include <vector>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/tuple/tuple.hpp>
#include "aligner_common.h"

class ExtensionSeedList {
public:
	enum State {
		kRun, kEnd, kIdle,
	};

	struct SeedData {
		int score;
		uint32_t query_id;
		uint32_t foward_seed_id;
		uint32_t reverse_seed_id;
		AlignerCommon::Coordinate hit;
	};

	struct CpuSeed {
		bool reverse;
		int score;
		uint32_t query_id;
		AlignerCommon::Coordinate coordinate;
	};

	ExtensionSeedList();
	virtual ~ExtensionSeedList();

	State GetState() {
		return state_;
	}

	void StartRun() {
		state_ = kRun;
	}
	void FinishRun() {
		state_ = kEnd;
	}

	void SetIdle() {
		state_ = kIdle;
	}

	void SetUp(size_t threhold_size, HostSeedsMemory* host_memory);
	HostSeedsMemory * GetSeedsMemory() {
		return seeds_memory_;
	}

	void ReleaseHostMemory() {
		seeds_memory_query_ids_ = NULL;
		seeds_memory_query_concatenated_positions_ = NULL;
		seeds_memory_database_positions_ = NULL;
		seeds_memory_temp_seed_ids_ = NULL;
		seeds_memory_->Release();
		seeds_memory_ = NULL;
	}

	size_t GetForwardGpuSeedsOffset() {
		return forward_offset_;
	}

	size_t GetReverseGpuSeedsOffset() {
		return reverse_offset_;
	}

	size_t GetForwardGpuSeedsLength() {
		return forward_length_;
	}
	size_t GetReverseGpuSeedsLength() {
		return reverse_length_;
	}

	size_t GetSeedsLength() {
		return seed_data_list_.size();
	}

	int PushExtensionSeeds(uint32_t query_id, uint32_t query_offset,
			std::vector<SeedSearcherCommon::Hit> &hits);

	void ConvertGappedExtensionSeeds(Queries &queries,
			std::vector<uint32_t>& query_offsets);

	void SetComputedFlagForBackSeeds(size_t length) {
		forward_length_ -= length;
		reverse_length_ -= length;
	}

	SeedData* GetSeedDataList();

	void Clear();
private:
	State state_;
	size_t forward_length_;
	size_t reverse_length_;
	size_t forward_offset_;
	size_t reverse_offset_;
	size_t seeds_memory_size_;
	size_t threhold_size_;
	HostSeedsMemory *seeds_memory_;
	uint32_t *seeds_memory_query_ids_;
	uint32_t *seeds_memory_query_concatenated_positions_;
	uint32_t *seeds_memory_database_positions_;
	uint32_t *seeds_memory_temp_seed_ids_;
	std::vector<SeedData> seed_data_list_;
	std::vector<uint32_t> temp_lengths_;
	//std::vector<uint32_t> temp_seed_ids_;
	std::vector<uint32_t> temp_query_ids_;
	std::vector<uint32_t> temp_query_concatenated_positions_;
	std::vector<uint32_t> temp_database_positions_;
};

#endif /* EXTENSION_SEED_LIST_H_ */
