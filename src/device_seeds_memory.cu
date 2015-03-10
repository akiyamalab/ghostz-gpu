/*
 * device_seed_buffer.cpp
 *
 *  Created on: Aug 5, 2014
 *      Author: shu
 */

#include <cub/cub.cuh>
#include "device_seeds_memory.h"

DeviceSeedsMemory::DeviceSeedsMemory() :
		size_(0), temp_storage_bytes(0) {

}

DeviceSeedsMemory::~DeviceSeedsMemory() {

}

int DeviceSeedsMemory::Init(size_t size) {
	size_ = size;
	query_seed_ids_.resize(size);
	query_seed_lengths_.resize(size);
	query_ids_.resize(size);
	query_concatenated_positions_.resize(size);
	database_positions_.resize(size);
	result_values_.resize(size);
	void *d_temp_storage = NULL;
	temp_storage_bytes = 0;
	cub::DeviceScan::InclusiveScan(d_temp_storage, temp_storage_bytes,
			thrust::raw_pointer_cast(query_seed_ids_.data()),
			thrust::raw_pointer_cast(query_seed_lengths_.data()), cub::Max(),
			size);
	temp_storage_.resize(temp_storage_bytes);
	return 0;
}

