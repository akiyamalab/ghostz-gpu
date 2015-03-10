/*
 * gapped_extension_seed_list.cpp
 *
 *  Created on: Aug 27, 2014
 *      Author: shu
 */

#include <assert.h>
#include <thrust/sort.h>
#include "extension_seed_list.h"

using namespace std;

ExtensionSeedList::ExtensionSeedList() :
		state_(kIdle), forward_length_(0), reverse_length_(0), forward_offset_(
				0), reverse_offset_(0), seeds_memory_size_(0), threhold_size_(
				0), seeds_memory_(NULL), seeds_memory_query_ids_(NULL), seeds_memory_query_concatenated_positions_(
				NULL), seeds_memory_database_positions_(NULL) {
}

ExtensionSeedList::~ExtensionSeedList() {
}

void ExtensionSeedList::SetUp(size_t threhold_size,
		HostSeedsMemory* host_memory) {
	assert(state_ == kIdle);
	forward_length_ = 0;
	reverse_length_ = 0;
	forward_offset_ = threhold_size;
	reverse_offset_ = 0;
	threhold_size_ = threhold_size;
	seeds_memory_size_ = 2 * threhold_size;
	seeds_memory_ = host_memory;
	seeds_memory_->Allocate();
	seeds_memory_query_ids_ = seeds_memory_->GetQueryIds();
	seeds_memory_query_concatenated_positions_ =
			seeds_memory_->GetQueryConcatenatedPositions();
	seeds_memory_database_positions_ = seeds_memory_->GetDatabasePositions();
	seeds_memory_temp_seed_ids_ = seeds_memory_->GetQuerySeedIds();
	seed_data_list_.clear();
	temp_lengths_.resize(seeds_memory_size_);
	temp_query_ids_.resize(seeds_memory_size_);
	temp_query_concatenated_positions_.resize(seeds_memory_size_);
	temp_database_positions_.resize(seeds_memory_size_);
}

int ExtensionSeedList::PushExtensionSeeds(uint32_t query_id,
		uint32_t query_offset, std::vector<SeedSearcherCommon::Hit> &hits) {
	if (seed_data_list_.size() + hits.size() < threhold_size_) {
		SeedData seed_data;
		seed_data.score = 0;
		for (size_t i = 0; i < hits.size(); ++i) {
			SeedSearcherCommon::Hit &hit = hits[i];
			seed_data.query_id = query_id;
			seed_data.hit.query_position = hit.query_sequence_position;
			seed_data.hit.database_position = hit.database_sequence_position;
			size_t reverse_p = GetReverseGpuSeedsOffset()
					+ GetReverseGpuSeedsLength();
			seeds_memory_query_ids_[reverse_p] = query_id;
			seeds_memory_query_concatenated_positions_[reverse_p] =
					hit.query_sequence_position + query_offset - 1;
			seeds_memory_database_positions_[reverse_p] =
					hit.database_sequence_position - 1;
			++reverse_length_;
			size_t foward_p = GetForwardGpuSeedsOffset()
					+ GetForwardGpuSeedsLength();
			seeds_memory_query_ids_[foward_p] = query_id;
			seeds_memory_query_concatenated_positions_[foward_p] =
					hit.query_sequence_position + query_offset;
			seeds_memory_database_positions_[foward_p] =
					hit.database_sequence_position;
			++forward_length_;
			seed_data.reverse_seed_id = reverse_p;
			seed_data.foward_seed_id = foward_p;
			seed_data_list_.push_back(seed_data);
		}
		return 0;
	} else {
		return 1;
	}
}

void ExtensionSeedList::ConvertGappedExtensionSeeds(Queries &queries,
		std::vector<uint32_t>& query_offsets) {
	if (seed_data_list_.empty()) {
		return;
	}
	uint32_t query_id = seed_data_list_[0].query_id;
	Query *query = queries.GetQuery(query_id);
	uint32_t query_offset = query_offsets[query_id];
	size_t number_seeds = GetSeedsLength();
	int *result_values = seeds_memory_->GetResultValues();
	forward_length_ = number_seeds;
	reverse_length_ = number_seeds;
	for (size_t seed_i = 0; seed_i < number_seeds; ++seed_i) {
		ExtensionSeedList::SeedData &seed_data = seed_data_list_[seed_i];
		if (query_id != seed_data.query_id) {
			query_id = seed_data.query_id;
			query = queries.GetQuery(query_id);
			query_offset = query_offsets[query_id];
		}

		size_t reverse_p = seed_data.reverse_seed_id;
		size_t foward_p = seed_data.foward_seed_id;
		seed_data.score += result_values[seed_data.reverse_seed_id];
		seed_data.score += result_values[seed_data.foward_seed_id];
		seeds_memory_temp_seed_ids_[reverse_p] = seed_i;
		seeds_memory_temp_seed_ids_[foward_p] = seed_i;
		temp_lengths_[reverse_p] =
				query->GetDistanceFromStartDelimiter(
						seeds_memory_query_concatenated_positions_[seed_data.reverse_seed_id]
								- query_offset);
		temp_lengths_[foward_p] =
				query->GetDistanceFromEndDelimiter(
						seeds_memory_query_concatenated_positions_[seed_data.foward_seed_id]
								- query_offset);
	}

	/*
	 memcpy(&temp_query_concatenated_positions_[GetReverseGpuSeedsOffset()],
	 &seeds_memory_query_concatenated_positions_[GetReverseGpuSeedsOffset()],
	 sizeof(temp_query_concatenated_positions_[0])
	 * GetReverseGpuSeedsLength());

	 memcpy(&temp_query_concatenated_positions_[GetForwardGpuSeedsOffset()],
	 &seeds_memory_query_concatenated_positions_[GetForwardGpuSeedsOffset()],
	 sizeof(temp_query_concatenated_positions_[0])
	 * GetForwardGpuSeedsLength());

	 memcpy(&temp_database_positions_[GetReverseGpuSeedsOffset()],
	 &seeds_memory_database_positions_[GetReverseGpuSeedsOffset()],
	 sizeof(temp_database_positions_[0]) * GetReverseGpuSeedsLength());

	 memcpy(&temp_database_positions_[GetForwardGpuSeedsOffset()],
	 &seeds_memory_database_positions_[GetForwardGpuSeedsOffset()],
	 sizeof(temp_database_positions_[0]) * GetForwardGpuSeedsLength());
	 */

	thrust::sort_by_key(temp_lengths_.begin() + GetReverseGpuSeedsOffset(),
			temp_lengths_.begin() + GetReverseGpuSeedsOffset()
					+ GetReverseGpuSeedsLength(),
			seeds_memory_temp_seed_ids_ + GetReverseGpuSeedsOffset());

	thrust::sort_by_key(temp_lengths_.begin() + GetForwardGpuSeedsOffset(),
			temp_lengths_.begin() + GetForwardGpuSeedsOffset()
					+ GetForwardGpuSeedsLength(),
			seeds_memory_temp_seed_ids_ + GetForwardGpuSeedsOffset());

	size_t reverse_end = GetReverseGpuSeedsOffset()
			+ GetReverseGpuSeedsLength();
	for (size_t i = GetReverseGpuSeedsOffset(); i < reverse_end; ++i) {
		uint32_t seed_id = seeds_memory_temp_seed_ids_[i];
		/*
		 seeds_memory_query_concatenated_positions_[i] =
		 temp_query_concatenated_positions_[seed_data_list_[seed_id].reverse_seed_id]
		 - 1;
		 seeds_memory_database_positions_[i] =
		 temp_database_positions_[seed_data_list_[seed_id].reverse_seed_id]
		 - 1;
		 */
		seed_data_list_[seed_id].reverse_seed_id = i;
#if 0
		if (seed_data_list_[seed_id].query_id == 9995) {
			cout << "length : " << temp_lengths_[i] << " : q "
			<< seed_data_list_[seed_id].query_id << " , "
			<< seeds_memory_query_concatenated_positions_[i] << ", d "
			<< seeds_memory_database_positions_[i] << endl;
		}
#endif
	}

	size_t foward_end = GetForwardGpuSeedsOffset() + GetForwardGpuSeedsLength();
	for (size_t i = GetForwardGpuSeedsOffset(); i < foward_end; ++i) {
		uint32_t seed_id = seeds_memory_temp_seed_ids_[i];
		/*
		 seeds_memory_query_concatenated_positions_[i] =
		 temp_query_concatenated_positions_[seed_data_list_[seed_id].foward_seed_id]
		 + 1;
		 seeds_memory_database_positions_[i] =
		 temp_database_positions_[seed_data_list_[seed_id].foward_seed_id]
		 + 1;
		 */
		seed_data_list_[seed_id].foward_seed_id = i;
#if 0
		if (seed_data_list_[seed_id].query_id == 9995) {
			cout << "length : " << temp_lengths_[i] << " : q "
			<< seed_data_list_[seed_id].query_id << " , "
			<< seeds_memory_query_concatenated_positions_[i] << ", d "
			<< seeds_memory_database_positions_[i] << endl;
		}
#endif
	}

	return;
}

ExtensionSeedList::SeedData * ExtensionSeedList::GetSeedDataList() {
	return &seed_data_list_[0];
}

void ExtensionSeedList::Clear() {
	state_ = kIdle;
	forward_length_ = 0;
	reverse_length_ = 0;
	seed_data_list_.clear();
}
