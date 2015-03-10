/*
 * distance_calculator_gpu.cpp
 *
 *  Created on: Jul 28, 2014
 *      Author: shu
 */

#include "distance_calculator_gpu.h"
#include "cuda_common.h"
#include "group_loader.h"
#include "packed_alphabet_code.h"
#include <assert.h>

namespace distance_calculator_gpu_kernel {
#if 0
const int debug_seed_i = 4;
const int debug_q_p = 119;
const int debug_db_p = 461379
+ cuda_common::kMaxLoadLength
* packed_alphabet_code::kNumberOfCodesInBlock;
#endif
static const size_t kNumberBlocks = 512;
static const size_t kNumberThreads = 128;
static const size_t kLoadLength = 4;

__device__ DistanceCalculatorGpu::Distance CalculateDisntace(
		const packed_alphabet_code::PackedAlphabetCode sequence0_code_block,
		const packed_alphabet_code::PackedAlphabetCode sequence1_code_block,
		const AlphabetCoder::Code* /* __restrict__ */redueced_code_map,
		const uint32_t length) {
	DistanceCalculatorGpu::Distance distance = 0;
	for (uint32_t i = 0; i < length; ++i) {
		const uint32_t s0_c = packed_alphabet_code::GetAlphabetCode<
				cuda_common::kFoward>(sequence0_code_block, i);
		const uint32_t s1_c = packed_alphabet_code::GetAlphabetCode<
				cuda_common::kFoward>(sequence1_code_block, i);
#if 0
		if ((blockDim.x * blockIdx.x + threadIdx.x) == debug_seed_i) {
			printf("%d, %d\n", s0_c, s1_c);
		}
#endif
#if 0
		if ((blockDim.x * blockIdx.x + threadIdx.x) == debug_seed_i) {
			printf("%d, %d\n", redueced_code_map[s0_c],
					redueced_code_map[s1_c]);
		}
#endif
		distance +=
				(__ldg(&redueced_code_map[s0_c])
						!= __ldg(&redueced_code_map[s1_c])) ? 1 : 0;
	}
	return distance;
}

__device__ DistanceCalculatorGpu::Distance CalculateDisntace(
		const packed_alphabet_code::PackedAlphabetCode sequence0_code_block,
		const packed_alphabet_code::PackedAlphabetCode sequence1_code_block,
		const AlphabetCoder::Code* /* __restrict__ */redueced_code_map) {
	DistanceCalculatorGpu::Distance distance = 0;

#pragma unroll
	for (uint32_t i = 0; i < packed_alphabet_code::kNumberOfCodesInBlock; ++i) {
		const uint32_t s0_c = packed_alphabet_code::GetAlphabetCode<
				cuda_common::kFoward>(sequence0_code_block, i);
		const uint32_t s1_c = packed_alphabet_code::GetAlphabetCode<
				cuda_common::kFoward>(sequence1_code_block, i);
#if 0
		if ((blockDim.x * blockIdx.x + threadIdx.x) == debug_seed_i) {
			printf("%d, %d\n", s0_c, s1_c);
		}
#endif
#if 0
		if ((blockDim.x * blockIdx.x + threadIdx.x) == debug_seed_i) {
			printf("%d, %d\n", redueced_code_map[s0_c],
					redueced_code_map[s1_c]);
		}
#endif
		distance +=
				(redueced_code_map[s0_c] != redueced_code_map[s1_c]) ? 1 : 0;
	}
	return distance;
}

__device__ DistanceCalculatorGpu::Distance CalculateDisntaceDevice(
		const packed_alphabet_code::PackedAlphabetCode* sequence0_code_block_cache_mem,
		const uint32_t s0_start_position,
		const packed_alphabet_code::PackedAlphabetCode* sequence1_code_block_cache_mem,
		const uint32_t s1_start_position,
		const AlphabetCoder::Code* /* __restrict__ */redueced_code_map,
		const uint32_t subsequence_length) {

	const packed_alphabet_code::PackedAlphabetCode * sequence0_code_block_cache =
			sequence0_code_block_cache_mem;
	const packed_alphabet_code::PackedAlphabetCode * sequence1_code_block_cache =
			sequence1_code_block_cache_mem;
	int s0_next_bolock_position = packed_alphabet_code::GetBlockPosition(
			s0_start_position);
	int s0_offset = s0_start_position
			- s0_next_bolock_position
					* packed_alphabet_code::kNumberOfCodesInBlock;
	int s0_next_bolock_chache_position = 0;
	uint32_t s0_block_sift;
	uint32_t s0_next_block_sift;
	packed_alphabet_code::PackedAlphabetCode s0_temp_code_block;
	packed_alphabet_code::InitPackedAlphabetCode<cuda_common::kFoward>(
			sequence0_code_block_cache, s0_offset,
			&s0_next_bolock_chache_position, &s0_block_sift,
			&s0_next_block_sift, &s0_temp_code_block);

	int s1_next_bolock_position = packed_alphabet_code::GetBlockPosition(
			s1_start_position);
	int s1_offset = s1_start_position
			- s1_next_bolock_position
					* packed_alphabet_code::kNumberOfCodesInBlock;
	int s1_next_bolock_chache_position = 0;
	uint32_t s1_block_sift;
	uint32_t s1_next_block_sift;
	packed_alphabet_code::PackedAlphabetCode s1_temp_code_block;
	packed_alphabet_code::InitPackedAlphabetCode<cuda_common::kFoward>(
			sequence1_code_block_cache, s1_offset,
			&s1_next_bolock_chache_position, &s1_block_sift,
			&s1_next_block_sift, &s1_temp_code_block);
#if 0
	if ((debug_db_p - 1) == s1_start_position || debug_db_p == s1_start_position
			|| (blockDim.x * blockIdx.x + threadIdx.x) == debug_seed_i) {
		printf("temp code block %ld, %ld\n", s0_temp_code_block,
				s1_temp_code_block);
	}
#endif

	DistanceCalculatorGpu::Distance distance = 0;
	uint32_t offset = 0;
	for (;
			offset + packed_alphabet_code::kNumberOfCodesInBlock
					<= subsequence_length; offset +=
					packed_alphabet_code::kNumberOfCodesInBlock) {
		packed_alphabet_code::PackedAlphabetCode s0_code_block =
				packed_alphabet_code::GetPackedAlphabetCode<cuda_common::kFoward>(
						sequence0_code_block_cache, s0_block_sift,
						s0_next_block_sift, &s0_next_bolock_chache_position,
						&s0_temp_code_block);

		packed_alphabet_code::PackedAlphabetCode s1_code_block =
				packed_alphabet_code::GetPackedAlphabetCode<cuda_common::kFoward>(
						sequence1_code_block_cache, s1_block_sift,
						s1_next_block_sift, &s1_next_bolock_chache_position,
						&s1_temp_code_block);

#if 0
		if ((debug_db_p - 1) == s1_start_position
				|| debug_db_p == s1_start_position
				|| (blockDim.x * blockIdx.x + threadIdx.x) == debug_seed_i) {
			printf("code block %ld, %ld\n", s0_code_block, s1_code_block);
		}
#endif
		distance += CalculateDisntace(s0_code_block, s1_code_block,
				redueced_code_map);
	}

	uint32_t remain_length = subsequence_length - offset;
	if (remain_length > 0) {
		packed_alphabet_code::PackedAlphabetCode s0_code_block =
				packed_alphabet_code::GetPackedAlphabetCode<cuda_common::kFoward>(
						sequence0_code_block_cache, s0_block_sift,
						s0_next_block_sift, &s0_next_bolock_chache_position,
						&s0_temp_code_block);

		packed_alphabet_code::PackedAlphabetCode s1_code_block =
				packed_alphabet_code::GetPackedAlphabetCode<cuda_common::kFoward>(
						sequence1_code_block_cache, s1_block_sift,
						s1_next_block_sift, &s1_next_bolock_chache_position,
						&s1_temp_code_block);

#if 0
		if ((debug_db_p - 1) == s1_start_position
				|| debug_db_p == s1_start_position
				|| (blockDim.x * blockIdx.x + threadIdx.x) == debug_seed_i) {
			printf("code block %ld, %ld\n", s0_code_block, s1_code_block);
		}
#endif
		distance += CalculateDisntace(s0_code_block, s1_code_block,
				redueced_code_map, remain_length);
	}
	return distance;
}

__global__ void CalculateDisntacesKernel(
		const packed_alphabet_code::PackedAlphabetCode* sequence0_code_block,
		const packed_alphabet_code::PackedAlphabetCode* sequence1_code_block,
		const uint32_t number_distances, const uint32_t * sequence0_positions,
		const uint32_t* sequence1_positions,
		const AlphabetCoder::Code* /* __restrict__ */redueced_code_map,
		const uint32_t subsequence_length,
		DistanceCalculatorGpu::Distance *distances) {
	assert(
			subsequence_length + packed_alphabet_code::kNumberOfCodesInBlock
					<= cuda_common::kMaxLoadLength
							* packed_alphabet_code::kNumberOfCodesInBlock);
	extern __shared__ char shared_mem[];
	uint32_t* s_sequence_positions = (uint32_t *) shared_mem;
	uint64_t* s_group_load_temp =
			(packed_alphabet_code::PackedAlphabetCode *) &s_sequence_positions[blockDim.x];
	GroupLoader<const packed_alphabet_code::PackedAlphabetCode *, kLoadLength> group_loader;
	int group_idx = group_loader.GetGroupId();
	int idx_in_group = group_loader.GetIdInGroup();
	int number_members = group_loader.GetNumberMembers();
	int group_work_skip = group_loader.GetNumberGroups();
	int number_group_in_block = blockDim.x / number_members;
	int number_group_works = (number_distances + number_members - 1)
			/ number_members;
	int remain_group_in_block = number_group_works % number_group_in_block;
	int group_work_end = number_group_works
			+ (remain_group_in_block != 0 ?
					(number_group_in_block - remain_group_in_block) : 0);
	uint32_t* s_group_sequence_positions = &s_sequence_positions[(threadIdx.x
			/ number_members) * number_members];
	const uint32_t foward_direction = subsequence_length / 2;
	for (int group_work_i = group_idx; group_work_i < group_work_end;
			group_work_i += group_work_skip) {

		const int group_extension_begin = group_work_i * number_members;
		const int thread_work_id = group_extension_begin + idx_in_group;
		const int extension_id =
				thread_work_id < number_distances ? thread_work_id : 0;
		packed_alphabet_code::PackedAlphabetCode s0_cache[kLoadLength];
		packed_alphabet_code::PackedAlphabetCode s1_cache[kLoadLength];

		const uint32_t sequence0_position = sequence0_positions[extension_id]
				+ cuda_common::kMaxLoadLength
						* packed_alphabet_code::kNumberOfCodesInBlock
				- foward_direction;

		const uint32_t sequence1_position = sequence1_positions[extension_id]
				+ cuda_common::kMaxLoadLength
						* packed_alphabet_code::kNumberOfCodesInBlock
				- foward_direction;

		uint32_t sequence0_brock_position =
				packed_alphabet_code::GetBlockPosition(sequence0_position);
		s_group_sequence_positions[idx_in_group] = sequence0_brock_position;
		//__syncthreads();
		group_loader.Load(sequence0_code_block, s_group_sequence_positions,
				min((int) number_distances - group_extension_begin,
						number_members), s_group_load_temp, s0_cache);
#if 0
		// debug /////////////////////////////////////
		if (debug_db_p == sequence1_position
				|| debug_seed_i == blockIdx.x * blockDim.x + threadIdx.x) {
			printf("position : %d\n", sequence0_positions[extension_id]);
			printf("block position : %d\n",
					s_group_sequence_positions[idx_in_group]);
			for (int i = 0; i < kLoadLength; ++i) {
				printf("%d : %ld\n", i, s0_cache[i]);
			}
		}
		/////////////////////////////////////////////
#endif
		uint32_t sequence1_brock_position =
				packed_alphabet_code::GetBlockPosition(sequence1_position);
		s_group_sequence_positions[idx_in_group] = sequence1_brock_position;
		//__syncthreads();
		group_loader.Load(sequence1_code_block, s_group_sequence_positions,
				min((int) number_distances - group_extension_begin,
						number_members), s_group_load_temp, s1_cache);
#if 0
		// debug /////////////////////////////////////
		if (debug_db_p == sequence1_position
				|| debug_seed_i == blockIdx.x * blockDim.x + threadIdx.x) {
			printf("position : %d\n", sequence1_positions[extension_id]);
			printf("block position : %d\n",
					s_group_sequence_positions[idx_in_group]);
			for (int i = 0; i < kLoadLength; ++i) {
				printf("%d : %ld\n", i, s1_cache[i]);
			}
		}
		/////////////////////////////////////////////
#endif

		if (thread_work_id < number_distances) {
			distances[extension_id] = CalculateDisntaceDevice(s0_cache,
					sequence0_position, s1_cache, sequence1_position,
					redueced_code_map, subsequence_length);
		}
	}
	return;
}

}

DistanceCalculatorGpu::DistanceCalculatorGpu() :
		subsequence_length_(0), d_concatenated_query_sequence_(NULL), d_database_sequence_(
				NULL), d_redueced_code_map_(NULL) {
}

DistanceCalculatorGpu::~DistanceCalculatorGpu() {

}

int DistanceCalculatorGpu::SetQueries(
		const packed_alphabet_code::PackedAlphabetCode *d_query_sequence) {
	d_concatenated_query_sequence_ = d_query_sequence;
	return 0;
}
int DistanceCalculatorGpu::SetDatabase(
		const packed_alphabet_code::PackedAlphabetCode *d_database_sequence) {
	d_database_sequence_ = d_database_sequence;
	return 0;
}
int DistanceCalculatorGpu::SetReducedCodeMap(
		const AlphabetCoder::Code *d_redueced_code_map) {
	d_redueced_code_map_ = d_redueced_code_map;
	return 0;
}

int DistanceCalculatorGpu::SetSubsequenceLength(uint32_t length) {
	subsequence_length_ = length;
	return 0;
}

int DistanceCalculatorGpu::CalculateDistancesAsync(size_t size,
		uint32_t* query_concatenated_centers, uint32_t* database_centers,
		Distance* distances, uint32_t* d_query_concatenated_centers,
		uint32_t* d_database_centers, Distance* d_distances,
		cudaStream_t &stream) {

	CUDA_CHECK_RETURN(
			cudaMemcpyAsync(d_query_concatenated_centers,
					query_concatenated_centers,
					sizeof(query_concatenated_centers[0]) * size,
					cudaMemcpyDefault, stream));

	CUDA_CHECK_RETURN(
			cudaMemcpyAsync(d_database_centers, database_centers,
					sizeof(database_centers[0]) * size, cudaMemcpyDefault,
					stream));

	size_t shared_memory_size =
			sizeof(uint32_t) * distance_calculator_gpu_kernel::kNumberThreads
					+ GroupLoader<
							const packed_alphabet_code::PackedAlphabetCode *,
							distance_calculator_gpu_kernel::kLoadLength>::GetTotalSharedMemorySize(
							distance_calculator_gpu_kernel::kNumberThreads);
	distance_calculator_gpu_kernel::CalculateDisntacesKernel<<<
	distance_calculator_gpu_kernel::kNumberBlocks,
	distance_calculator_gpu_kernel::kNumberThreads, shared_memory_size,
	stream>>>(d_concatenated_query_sequence_, d_database_sequence_,
			size, d_query_concatenated_centers, d_database_centers,
			d_redueced_code_map_, subsequence_length_, d_distances);

	CUDA_CHECK_RETURN(
			cudaMemcpyAsync(distances, d_distances, sizeof(distances[0]) * size,
					cudaMemcpyDefault, stream));

	return 0;
}
