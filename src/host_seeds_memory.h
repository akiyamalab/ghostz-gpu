/*
 * host_seed_buffer.h
 *
 *  Created on: Aug 5, 2014
 *      Author: shu
 */

#ifndef HOST_SEED_BUFFER_H_
#define HOST_SEED_BUFFER_H_

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/system/cuda/experimental/pinned_allocator.h>
#include <thrust/copy.h>
#include "cuda_common.h"
#include <assert.h>

class HostSeedsMemory {
public:
	HostSeedsMemory();
	HostSeedsMemory(size_t id, size_t size);
	virtual ~HostSeedsMemory();

	int Init(size_t id, size_t size);

	size_t GetId() {
		return id_;
	}
	size_t GetSize() const;

	bool IsAllocated() {
		return allocated_;
	}
	void Allocate() {
		assert(allocated_ == false);
		allocated_ = true;
	}

	void Release() {
		assert(allocated_ == true);
		allocated_ = false;
	}

	uint32_t* GetQuerySeedIds() {
		return &query_seed_ids_[0];
	}

	uint32_t* GetQuerySeedStarts() {
		return &query_seed_lengths_[0];
	}

	uint32_t* GetQueryIds() {
		return &query_ids_[0];
	}
	uint32_t* GetQueryConcatenatedPositions() {
		return &query_concatenated_positions_[0];
	}
	uint32_t* GetDatabasePositions() {
		return &database_positions_[0];
	}
	int* GetResultValues() {
		return &result_values_[0];
	}

private:
	bool allocated_;
	size_t id_;
	size_t size_;
	thrust::host_vector<uint32_t,
			thrust::cuda::experimental::pinned_allocator<uint32_t> > query_seed_ids_;
	thrust::host_vector<uint32_t,
			thrust::cuda::experimental::pinned_allocator<uint32_t> > query_seed_lengths_;
	thrust::host_vector<uint32_t,
			thrust::cuda::experimental::pinned_allocator<uint32_t> > query_ids_;
	thrust::host_vector<uint32_t,
			thrust::cuda::experimental::pinned_allocator<uint32_t> > query_concatenated_positions_;
	thrust::host_vector<uint32_t,
			thrust::cuda::experimental::pinned_allocator<uint32_t> > database_positions_;
	thrust::host_vector<int, thrust::cuda::experimental::pinned_allocator<int> > result_values_;
};

#endif /* HOST_SEED_BUFFER_H_ */
