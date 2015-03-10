/*
 * ungapped_extension_with_trigger_seed_list.cpp
 *
 *  Created on: Sep 5, 2014
 *      Author: shu
 */

#include "ungapped_extension_with_trigger_seed_list.h"
#include <assert.h>

UngappedExtensionWithTriggerSeedList::UngappedExtensionWithTriggerSeedList() :
		length_(0), seeds_length_threhold_(0), query_seeds_length_(0), seeds_memory_(
				NULL), seeds_memory_query_seed_ids_(NULL), seeds_memory_query_seed_starts_(
				NULL),
		//seeds_memory_query_ids_(NULL), seeds_memory_query_concatenated_positions_(NULL),
		seeds_memory_database_positions_(NULL) {
	// TODO Auto-generated constructor stub

}

UngappedExtensionWithTriggerSeedList::~UngappedExtensionWithTriggerSeedList() {
	// TODO Auto-generated destructor stub
}

void UngappedExtensionWithTriggerSeedList::SetUp(size_t threhold_size,
		HostSeedsMemory* buffer) {
	seeds_memory_ = buffer;
	seeds_memory_->Allocate();
	length_ = 0;
	size_t prev_seeds_length_threshold = seeds_length_threhold_;
	seeds_length_threhold_ = threhold_size;
	query_seeds_length_ = 0;
	seeds_memory_query_seed_ids_ = seeds_memory_->GetQuerySeedIds();
	seeds_memory_query_seed_starts_ = seeds_memory_->GetQuerySeedStarts();
	seeds_memory_database_positions_ = seeds_memory_->GetDatabasePositions();
	if (over_database_positions_.size() > 0) {
		length_ = over_database_positions_.size();
		size_t copied_size = std::min(threhold_size,
				over_database_positions_.size());
		seeds_memory_query_seed_ids_[0] = over_query_seed_ids_[0];
		seeds_memory_query_seed_starts_[0] = 0;
		size_t next_copied_query_seed_id = 1;
		for (;
				next_copied_query_seed_id < over_query_seed_ids_.size()
						&& (over_query_seed_starts_[next_copied_query_seed_id]
								- prev_seeds_length_threshold) < copied_size;
				++next_copied_query_seed_id) {
			seeds_memory_query_seed_ids_[next_copied_query_seed_id] =
					over_query_seed_ids_[next_copied_query_seed_id];
			seeds_memory_query_seed_starts_[next_copied_query_seed_id] =
					over_query_seed_starts_[next_copied_query_seed_id]
							- prev_seeds_length_threshold;
		}
		query_seeds_length_ = next_copied_query_seed_id;
		memcpy(seeds_memory_database_positions_, &over_database_positions_[0],
				sizeof(over_database_positions_[0]) * copied_size);
		if (copied_size != over_database_positions_.size()) {
			size_t remain_size = over_database_positions_.size() - copied_size;
			size_t query_seed_start_offset = prev_seeds_length_threshold
					+ copied_size;

			if (next_copied_query_seed_id >= over_query_seed_ids_.size()
					|| 0
							< (over_query_seed_starts_[next_copied_query_seed_id]
									- prev_seeds_length_threshold)) {
				--next_copied_query_seed_id;
				over_query_seed_starts_[next_copied_query_seed_id] =
						query_seed_start_offset;
			}
			size_t new_query_seed_size = over_query_seed_ids_.size()
					- next_copied_query_seed_id;
			for (size_t i = 0;
					next_copied_query_seed_id < over_query_seed_ids_.size();
					++i, ++next_copied_query_seed_id) {

				over_query_seed_ids_[i] =
						over_query_seed_ids_[next_copied_query_seed_id];
				over_query_seed_starts_[i] =
						over_query_seed_starts_[next_copied_query_seed_id]
								- query_seed_start_offset;
			}

			for (size_t i = 0; i < remain_size; ++i) {
				over_database_positions_[i] =
						over_database_positions_[copied_size + i];
			}
			over_query_seed_ids_.resize(new_query_seed_size);
			over_query_seed_starts_.resize(new_query_seed_size);
			over_database_positions_.resize(remain_size);
		} else {
			ClearOverSeeds();
		}
	}
}
