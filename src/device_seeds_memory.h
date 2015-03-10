/*
 * device_seed_buffer.h
 *
 *  Created on: Aug 5, 2014
 *      Author: shu
 */

#ifndef DEVICE_SEED_BUFFER_H_
#define DEVICE_SEED_BUFFER_H_

#include <stdint.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/system/cuda/experimental/pinned_allocator.h>
#include <thrust/copy.h>
#include "cuda_common.h"

class DeviceSeedsMemory {
public:
	DeviceSeedsMemory();
	virtual ~DeviceSeedsMemory();

	int Init(size_t size);
	size_t GetSize();
	size_t GetTempStorageSize() {
		return temp_storage_bytes;
	}

	void* GetTempStorage() {
		return (void*) thrust::raw_pointer_cast(temp_storage_.data());
	}

	uint32_t* GetQuerySeedIds() {
		return thrust::raw_pointer_cast(query_seed_ids_.data());
	}

	uint32_t* GetQuerySeedStarts() {
		return thrust::raw_pointer_cast(query_seed_lengths_.data());
	}

	uint32_t* GetQueryIds() {
		return thrust::raw_pointer_cast(query_ids_.data());
	}
	uint32_t* GetQueryConcatenatedPositions() {
		return thrust::raw_pointer_cast(query_concatenated_positions_.data());
	}
	uint32_t* GetDatabasePositions() {
		return thrust::raw_pointer_cast(database_positions_.data());
	}
	int* GetResultValues() {
		return thrust::raw_pointer_cast(result_values_.data());
	}
private:
	size_t size_;
	size_t temp_storage_bytes;
	thrust::device_vector<char> temp_storage_;
	thrust::device_vector<uint32_t> query_seed_ids_;
	thrust::device_vector<uint32_t> query_seed_lengths_;
	thrust::device_vector<uint32_t> query_ids_;
	thrust::device_vector<uint32_t> query_concatenated_positions_;
	thrust::device_vector<uint32_t> database_positions_;
	thrust::device_vector<int> result_values_;
}
;

#endif /* DEVICE_SEED_BUFFER_H_ */
