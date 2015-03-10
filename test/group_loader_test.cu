#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <stdint.h>
#include <thrust/sort.h>
#include <assert.h>
#include <gtest/gtest.h>
#include "../src/group_loader.h"
#include "../src/cuda_common.h"

__global__ void GroupLoadTestKernel(uint64_t* d_in, uint32_t* d_in_pos,
		uint64_t* d_out, int N) {
	const size_t length = 4;
	extern __shared__ char s[];
	uint64_t* s_in = (uint64_t*) s;
	uint64_t thread_data[length];
	GroupLoader<uint64_t*, length> group_loader;
	int group_idx = group_loader.GetGroupId();
	int idx_in_group = group_loader.GetIdInGroup();
	int number_members = group_loader.GetNumberMembers();
	int group_work_stride = group_loader.GetNumberGroups();
	int group_work_end = (N + number_members - 1) / number_members;
	for (int group_work_i = group_idx; group_work_i < group_work_end;
			group_work_i += group_work_stride) {
		int pos_idx_begin = group_work_i * number_members;
		group_loader.Load(d_in, &d_in_pos[pos_idx_begin],
				min(N - pos_idx_begin, number_members), s_in, thread_data);
		int subwork_start = pos_idx_begin;
		int output_pos = subwork_start + idx_in_group;
		uint64_t sum = 0;
		if (output_pos < N) {
			for (int i = 0; i < length; ++i) {
				sum += thread_data[i];
			}
			d_out[group_work_i * number_members + idx_in_group] = sum;
		}

	}
}

using namespace std;

class GroupLoaderTest: public ::testing::Test {
protected:
	virtual void SetUp() {
	}

	virtual void TearDown() {

	}
};

TEST_F(GroupLoaderTest, SharedMemorySize) {
	const int length = 4;
	const int number_threads = 128;
	size_t shared_memory_size =
			GroupLoader<uint64_t *, length>::GetTotalSharedMemorySize(
					number_threads);
	EXPECT_EQ(sizeof(uint64_t) * length * number_threads, shared_memory_size);
}

TEST_F(GroupLoaderTest, GroupLoad) {
	uint64_t *d_in = NULL;
	uint64_t *d_out = NULL;
	uint32_t *d_pos = NULL;
	const int load_length = 4;
	uint32_t data_size = 10;
	uint32_t work_size = data_size - load_length;
	vector<uint64_t> idata(data_size), random_pos(data_size),
			odata(data_size);
	vector<uint32_t> pos(data_size);

	for (int i = 0; i < data_size; i++) {
		idata[i] = (unsigned int) i;
		pos[i] = i;
		random_pos[i] = rand();
	}

	CUDA_CHECK_RETURN(cudaMalloc((void** ) &d_in, sizeof(uint64_t) * data_size));
	CUDA_CHECK_RETURN(
			cudaMalloc((void** ) &d_pos, sizeof(uint32_t) * data_size));
	CUDA_CHECK_RETURN(
			cudaMalloc((void** ) &d_out, sizeof(uint64_t) * data_size));

	thrust::sort_by_key(random_pos.begin(), random_pos.begin() + work_size,
			pos.begin());
	CUDA_CHECK_RETURN(
			cudaMemcpy(d_pos, &pos[0], sizeof(uint32_t) * data_size,
					cudaMemcpyHostToDevice));
	CUDA_CHECK_RETURN(
			cudaMemcpy(d_in, &idata[0], sizeof(uint64_t) * data_size,
					cudaMemcpyHostToDevice));
	int threads = 128;
	int blocks = min((work_size + threads - 1) / threads, 128);
	//device_copy_ref_kernel<<<blocks, threads>>>(d_in, d_in_pos, d_out, N);
	GroupLoadTestKernel<<<blocks, threads, GroupLoader<uint64_t *, load_length>::GetTotalSharedMemorySize(
			threads)>>>(d_in,d_pos, d_out, work_size);

	CUDA_CHECK_RETURN(cudaThreadSynchronize());	// Wait for the GPU launched work to complete
	CUDA_CHECK_RETURN(cudaGetLastError());
	CUDA_CHECK_RETURN(
			cudaMemcpy(&odata[0], d_out, sizeof(uint64_t) * data_size,
					cudaMemcpyDeviceToHost));

	for (int i = 0; i < work_size; ++i) {
		uint64_t sum = 0;
		uint64_t p = pos[i];
		for (int j = 0; j < load_length; ++j) {
			sum += idata[p + j];
		}
		EXPECT_EQ(sum, odata[i]);
	}
	CUDA_CHECK_RETURN(cudaFree((void* ) d_in));
	CUDA_CHECK_RETURN(cudaFree((void* ) d_pos));
	CUDA_CHECK_RETURN(cudaFree((void* ) d_out));
}
