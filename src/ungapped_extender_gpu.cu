/*
 * ungapped_extender_gpu.cpp
 *
 *  Created on: Jul 21, 2014
 *      Author: shu
 */

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/system/cuda/experimental/pinned_allocator.h>
#include <thrust/copy.h>
#include <vector>
#include <assert.h>
#include <cstdlib>
#include <cub/cub.cuh>

#include "ungapped_extender_gpu.h"
#include "group_loader.h"
#include "packed_alphabet_code.h"
#include "score_matrix.h"
#include "cuda_common.h"

using namespace std;

namespace ungapped_extention_with_trigger_gpu_kernel {
#if 0
const int debug_q_p = 119;
const int debug_db_p = 8963
+ cuda_common::kMaxLoadLength
* packed_alphabet_code::kNumberOfCodesInBlock;
#endif
static const size_t kNumberBlocks = 128;
static const size_t kNumberThreads = 128;
static const size_t kLoadLength = 4;

template<cuda_common::Direction TDirection>
__device__ int UngappedExtendOneSideMultiAlphabetCode(
		packed_alphabet_code::PackedAlphabetCode sequence0_code_block,
		packed_alphabet_code::PackedAlphabetCode sequence1_code_block,
		const uint32_t sequence_delimiter,
		const int* /* __restrict__  */score_matrix,
		const uint32_t number_letters, const int cutoff, const int trigger,
		int *score_ptr, int *best_score_ptr) {

	int score = *score_ptr;
	int best_score = *best_score_ptr;
	int threshold = best_score - cutoff;
	uint32_t stop_flag = 0;

#pragma unroll
	for (uint32_t i = 0; i < packed_alphabet_code::kNumberOfCodesInBlock; ++i) {
		const uint32_t s0_c = packed_alphabet_code::GetAlphabetCode<TDirection>(
				sequence0_code_block, i);
		const uint32_t s1_c = packed_alphabet_code::GetAlphabetCode<TDirection>(
				sequence1_code_block, i);
#if 0
		if ((blockDim.x * blockIdx.x + threadIdx.x) == 0) {
			printf("%d, %d\n", s0_c, s1_c);
		}
#endif
		stop_flag = stop_flag | s0_c == sequence_delimiter
				| s1_c == sequence_delimiter;
		score +=
				stop_flag ?
						0 : __ldg(&score_matrix[s0_c * number_letters + s1_c]);
		best_score = max(best_score, score);
		threshold = best_score - cutoff;
		stop_flag = stop_flag | score <= threshold | best_score > trigger;
	}

	*score_ptr = score;
	*best_score_ptr = best_score;

	return stop_flag;
}

template<cuda_common::Direction TDirection>
__device__ int UngappedExtendOneSideDevice(
		const packed_alphabet_code::PackedAlphabetCode * sequence0_code_block,
		const uint32_t s0_start_position,
		const packed_alphabet_code::PackedAlphabetCode * sequence1_code_block,
		const uint32_t s1_start_position,
		const packed_alphabet_code::PackedAlphabetCode * sequence0_code_block_cache_mem,
		const packed_alphabet_code::PackedAlphabetCode * sequence1_code_block_cache_mem,
		const int cache_length, const AlphabetCoder::Code sequence_delimiter,
		const int* /* __restrict__  */score_matrix,
		const uint32_t number_letters, const int cutoff, const int trigger,
		int current_best_score) {

	int score = current_best_score;
	int best_score = current_best_score;

	const packed_alphabet_code::PackedAlphabetCode * sequence0_code_block_cache =
			TDirection == cuda_common::kReverse ?
					sequence0_code_block_cache_mem + cache_length - 1 :
					sequence0_code_block_cache_mem;
	const packed_alphabet_code::PackedAlphabetCode * sequence1_code_block_cache =
			TDirection == cuda_common::kReverse ?
					sequence1_code_block_cache_mem + cache_length - 1 :
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
	packed_alphabet_code::InitPackedAlphabetCode<TDirection>(
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
	packed_alphabet_code::InitPackedAlphabetCode<TDirection>(
			sequence1_code_block_cache, s1_offset,
			&s1_next_bolock_chache_position, &s1_block_sift,
			&s1_next_block_sift, &s1_temp_code_block);
#if 0
	if ((debug_db_p - 1) == s1_start_position || debug_db_p == s1_start_position
			|| 0 == blockIdx.x * blockDim.x + threadIdx.x) {
		printf("temp code block %ld, %ld\n", s0_temp_code_block,
				s1_temp_code_block);
	}
#endif

	int ret = 1;
	int cache_end = (cache_length)
			* (TDirection == cuda_common::kReverse ? -1 : 1);
	while (ret && (s0_next_bolock_chache_position != cache_end)) {
		packed_alphabet_code::PackedAlphabetCode s0_code_block =
				packed_alphabet_code::GetPackedAlphabetCode<TDirection>(
						sequence0_code_block_cache, s0_block_sift,
						s0_next_block_sift, &s0_next_bolock_chache_position,
						&s0_temp_code_block);

		packed_alphabet_code::PackedAlphabetCode s1_code_block =
				packed_alphabet_code::GetPackedAlphabetCode<TDirection>(
						sequence1_code_block_cache, s1_block_sift,
						s1_next_block_sift, &s1_next_bolock_chache_position,
						&s1_temp_code_block);

#if 0
		if ((debug_db_p - 1) == s1_start_position
				|| debug_db_p == s1_start_position
				|| 0 == blockIdx.x * blockDim.x + threadIdx.x) {
			printf("code block %ld, %ld\n", s0_code_block, s1_code_block);
		}
#endif
		ret = !UngappedExtendOneSideMultiAlphabetCode<TDirection>(s0_code_block,
				s1_code_block, sequence_delimiter, score_matrix, number_letters,
				cutoff, trigger, &score, &best_score);
	}

	s0_next_bolock_position += s0_next_bolock_chache_position;
	s1_next_bolock_position += s1_next_bolock_chache_position;
	while (ret) {
		packed_alphabet_code::PackedAlphabetCode s0_code_block =
				packed_alphabet_code::GetPackedAlphabetCode<TDirection>(
						sequence0_code_block, s0_block_sift, s0_next_block_sift,
						&s0_next_bolock_position, &s0_temp_code_block);

		packed_alphabet_code::PackedAlphabetCode s1_code_block =
				packed_alphabet_code::GetPackedAlphabetCode<TDirection>(
						sequence1_code_block, s1_block_sift, s1_next_block_sift,
						&s1_next_bolock_position, &s1_temp_code_block);

#if 0
		if ((debug_db_p - 1) == s1_start_position
				|| debug_db_p == s1_start_position
				|| 0 == blockIdx.x * blockDim.x + threadIdx.x) {
			printf("code_block %ld, %ld\n", s0_code_block, s1_code_block);
		}
#endif
		ret = !UngappedExtendOneSideMultiAlphabetCode<TDirection>(s0_code_block,
				s1_code_block, sequence_delimiter, score_matrix, number_letters,
				cutoff, trigger, &score, &best_score);
	}

	return best_score;
}

__global__ void
__launch_bounds__(ungapped_extention_with_trigger_gpu_kernel::kNumberThreads, 1) UngappedExtendKernel(
		const packed_alphabet_code::PackedAlphabetCode* sequence0_code_block,
		const packed_alphabet_code::PackedAlphabetCode* sequence1_code_block,
		const uint32_t number_extensions,
		const uint32_t* /* __restrict__  */query_ids,
		const uint32_t * sequence0_positions,
		const uint32_t* sequence1_positions,
		const AlphabetCoder::Code sequence_delimiter,
		const int* /* __restrict__  */score_matrix,
		const uint32_t number_letters, const int* /* __restrict__  */cutoffs,
		const int* /* __restrict__  */triggers,
		char* /* __restrict__  */next_step_flags) {

	extern __shared__ char shared_mem[];
	uint32_t* s_sequence_positions = (uint32_t *) shared_mem;
	uint64_t* s_group_loader =
			(packed_alphabet_code::PackedAlphabetCode *) &s_sequence_positions[blockDim.x];
	GroupLoader<const packed_alphabet_code::PackedAlphabetCode *, kLoadLength> group_loader;
	int group_idx = group_loader.GetGroupId();
	int idx_in_group = group_loader.GetIdInGroup();
	int number_members = group_loader.GetNumberMembers();
	int group_work_skip = group_loader.GetNumberGroups();
	int number_group_in_block = blockDim.x / number_members;
	int number_group_works = (number_extensions + number_members - 1)
			/ number_members;
	int remain_group_in_block = number_group_works % number_group_in_block;
	int group_work_end = number_group_works
			+ (remain_group_in_block != 0 ?
					(number_group_in_block - remain_group_in_block) : 0);
	uint32_t* s_group_sequence_positions = &s_sequence_positions[(threadIdx.x
			/ number_members) * number_members];

	for (int group_work_i = group_idx; group_work_i < group_work_end;
			group_work_i += group_work_skip) {

		const int group_extension_begin = group_work_i * number_members;
		const int thread_work_id = group_extension_begin + idx_in_group;
		const int extension_id =
				thread_work_id < number_extensions ? thread_work_id : 0;
		packed_alphabet_code::PackedAlphabetCode s0_cache[kLoadLength];
		packed_alphabet_code::PackedAlphabetCode s1_cache[kLoadLength];

		const uint32_t sequence0_position = sequence0_positions[extension_id]
				+ cuda_common::kMaxLoadLength
						* packed_alphabet_code::kNumberOfCodesInBlock;

		const uint32_t sequence1_position = sequence1_positions[extension_id]
				+ cuda_common::kMaxLoadLength
						* packed_alphabet_code::kNumberOfCodesInBlock;

		// reverse
		uint32_t sequence0_brock_position =
				packed_alphabet_code::GetBlockPosition(sequence0_position - 1)
						- kLoadLength + 1;
		s_group_sequence_positions[idx_in_group] = sequence0_brock_position;
		//__syncthreads();
		group_loader.Load(sequence0_code_block, s_group_sequence_positions,
				min((int) number_extensions - group_extension_begin,
						number_members), s_group_loader, s0_cache);
#if 0
		// debug /////////////////////////////////////
		if (debug_db_p == sequence1_position
				|| 0 == blockIdx.x * blockDim.x + threadIdx.x) {
			printf("block position : %d\n",
					s_group_sequence_positions[idx_in_group]);
			for (int i = 0; i < kLoadLength; ++i) {
				printf("%d : %ld\n", i, s0_cache[i]);
			}
		}
		/////////////////////////////////////////////
#endif
		uint32_t sequence1_brock_position =
				packed_alphabet_code::GetBlockPosition(sequence1_position - 1)
						- kLoadLength + 1;
		s_group_sequence_positions[idx_in_group] = sequence1_brock_position;
		//__syncthreads();
		group_loader.Load(sequence1_code_block, s_group_sequence_positions,
				min((int) number_extensions - group_extension_begin,
						number_members), s_group_loader, s1_cache);
#if 0
		// debug /////////////////////////////////////
		if (debug_db_p == sequence1_position
				|| 0 == blockIdx.x * blockDim.x + threadIdx.x) {
			printf("block position : %d\n",
					s_group_sequence_positions[idx_in_group]);
			for (int i = 0; i < kLoadLength; ++i) {
				printf("%d : %ld\n", i, s1_cache[i]);
			}
		}
		/////////////////////////////////////////////
#endif

		uint32_t query_id = 0;
		int cutoff = 0;
		int trigger = 0;
		int current_best_score = 0;
		if (thread_work_id < number_extensions) {
			query_id = query_ids[extension_id];
			cutoff = cutoffs[query_id];
			trigger = triggers[query_id];
			current_best_score = UngappedExtendOneSideDevice<
					cuda_common::kReverse>(sequence0_code_block,
					sequence0_position - 1, sequence1_code_block,
					sequence1_position - 1, s0_cache, s1_cache, kLoadLength,
					sequence_delimiter, score_matrix, number_letters, cutoff,
					trigger, current_best_score);

			// foward
			uint32_t sequence0_brock_position =
					packed_alphabet_code::GetBlockPosition(sequence0_position);
			s_group_sequence_positions[idx_in_group] = sequence0_brock_position;
			//__syncthreads();
			group_loader.Load(sequence0_code_block, s_group_sequence_positions,
					min((int) number_extensions - group_extension_begin,
							number_members), s_group_loader, s0_cache);

			uint32_t sequence1_brock_position =
					packed_alphabet_code::GetBlockPosition(sequence1_position);
			s_group_sequence_positions[idx_in_group] = sequence1_brock_position;
			//__syncthreads();
			group_loader.Load(sequence1_code_block, s_group_sequence_positions,
					min((int) number_extensions - group_extension_begin,
							number_members), s_group_loader, s1_cache);

			if (current_best_score <= trigger) {
				current_best_score = UngappedExtendOneSideDevice<
						cuda_common::kFoward>(sequence0_code_block,
						sequence0_position, sequence1_code_block,
						sequence1_position, s0_cache, s1_cache, kLoadLength,
						sequence_delimiter, score_matrix, number_letters,
						cutoff, trigger, current_best_score);
			}
			next_step_flags[extension_id] =
					current_best_score > trigger ? 1 : 0;
		}
	}
	return;
}

__global__ void SetQuerySeedIdsIAtStart(const uint32_t query_seed_ids_length,
		const uint32_t* query_seed_starts, uint32_t set_query_seed_ids_length,
		uint32_t* set_query_seed_ids) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int skip = blockDim.x * gridDim.x;
	for (int i = idx; i < query_seed_ids_length; i += skip) {
		uint32_t start = query_seed_starts[i];
		set_query_seed_ids[start] = i;
	}

#if 0
	if (idx == 0) {
		for (int i = 30735 - 10; i < 30735 + 10; ++i) {
			printf("%d : %d\n", i, set_query_seed_ids[i]);
		}
	}
#endif
}

__global__ void SetQuerySeeds(const uint32_t query_seeds_length,
		const uint32_t* query_seed_query_ids,
		const uint32_t* query_seed_query_positions,
		const uint32_t* query_seed_ids, const uint32_t* query_seed_id_i_list,
		uint32_t* query_ids, uint32_t* query_positions) {
	const size_t vector_length = 2;
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int skip = blockDim.x * gridDim.x;

#if 0
	if (idx == 0) {
		for (int i = 30735 - 10; i < 30735 + 10; ++i) {
			printf("%d : %d : %d\n", i, query_seed_id_i_list[i],
					query_seed_ids[query_seed_id_i_list[i]]);
		}
	}
#endif
	for (int i = idx; i < query_seeds_length / vector_length; i += skip) {
		uint2 temp_vec;
		reinterpret_cast<uint64_t*>(&temp_vec)[0] =
				reinterpret_cast<const uint64_t*>(query_seed_id_i_list)[i];
		temp_vec.x = query_seed_ids[temp_vec.x];
		temp_vec.y = query_seed_ids[temp_vec.y];
#if 0
		if (query_ids[i] != query_seed_query_ids[id]) {
			printf("i %d\n", i);
			printf("id %d\n", id);
			printf("e %d ,a %d\n", query_ids[i], query_seed_query_ids[id]);
		}
		if (query_positions[i] != query_seed_query_positions[id]) {
			printf("i %d\n", i);
			printf("id %d\n", id);
			printf("e %d ,a %d\n", query_positions[i],
					query_seed_query_positions[id]);
		}
		assert(query_ids[i] == query_seed_query_ids[id]);
		assert(query_positions[i] == query_seed_query_positions[id]);
#endif
		uint2 temp_out = make_uint2(query_seed_query_ids[temp_vec.x],
				query_seed_query_ids[temp_vec.y]);
		reinterpret_cast<uint64_t*>(query_ids)[i] =
				reinterpret_cast<uint64_t*>(&temp_out)[0];
		//query_ids[i] = __ldg(&query_seed_query_ids[id]);

		temp_out = make_uint2(query_seed_query_positions[temp_vec.x],
				query_seed_query_positions[temp_vec.y]);
		reinterpret_cast<uint64_t*>(query_positions)[i] =
				reinterpret_cast<uint64_t*>(&temp_out)[0];
		//query_positions[i] = __ldg(&query_seed_query_positions[id]);
	}

	int i = idx + (query_seeds_length / vector_length) * vector_length;
	if (i < query_seeds_length) {
		uint32_t id = query_seed_ids[query_seed_id_i_list[i]];
		query_ids[i] = query_seed_query_ids[id];
		query_positions[i] = query_seed_query_positions[id];
	}
}

}

namespace ungapped_extention_gpu_kernel {
#if 0
const int debug_q_p = 119;
const int debug_db_p = 8963
+ cuda_common::kMaxLoadLength
* packed_alphabet_code::kNumberOfCodesInBlock;
#endif
static const size_t kNumberBlocks = 512;
static const size_t kNumberThreads = 128;
static const size_t kLoadLength = 4;

template<cuda_common::Direction TDirection>
__device__ int UngappedExtendOneSideMultiAlphabetCode(
		packed_alphabet_code::PackedAlphabetCode sequence0_code_block,
		packed_alphabet_code::PackedAlphabetCode sequence1_code_block,
		const uint32_t sequence_delimiter,
		const int* /* __restrict__  */score_matrix,
		const uint32_t number_letters, const int cutoff, const int offset,
		int *score_ptr, int *best_score_ptr,
		int *max_score_sequence_position_ptr) {

	int score = *score_ptr;
	int best_score = *best_score_ptr;
	int threshold = best_score - cutoff;
	int best_i = -1;
	uint32_t stop_flag = 0;

#pragma unroll
	for (uint32_t i = 0; i < packed_alphabet_code::kNumberOfCodesInBlock; ++i) {
		const uint32_t s0_c = packed_alphabet_code::GetAlphabetCode<TDirection>(
				sequence0_code_block, i);
		const uint32_t s1_c = packed_alphabet_code::GetAlphabetCode<TDirection>(
				sequence1_code_block, i);
#if 0
		if ((blockDim.x * blockIdx.x + threadIdx.x) == 0) {
			printf("%d, %d\n", s0_c, s1_c);
		}
#endif
		stop_flag = stop_flag | s0_c == sequence_delimiter
				| s1_c == sequence_delimiter;
		score += stop_flag ? 0 : score_matrix[s0_c * number_letters + s1_c];
		if (score > best_score) {
			best_i = i;
			best_score = score;
			threshold = best_score - cutoff;
		}
		stop_flag = stop_flag | score <= threshold;
	}

	*score_ptr = score;
	*best_score_ptr = best_score;
	if (best_i >= 0) {
		*max_score_sequence_position_ptr = offset + best_i;
	}

	return stop_flag;
}

template<cuda_common::Direction TDirection>
__device__ int UngappedExtendOneSideDevice(
		const packed_alphabet_code::PackedAlphabetCode * sequence0_code_block,
		const uint32_t s0_start_position,
		const packed_alphabet_code::PackedAlphabetCode * sequence1_code_block,
		const uint32_t s1_start_position,
		const packed_alphabet_code::PackedAlphabetCode * sequence0_code_block_cache_mem,
		const packed_alphabet_code::PackedAlphabetCode * sequence1_code_block_cache_mem,
		const int cache_length, const AlphabetCoder::Code sequence_delimiter,
		const int* /* __restrict__  */score_matrix,
		const uint32_t number_letters, const int cutoff, int *best_score_ptr,
		int *best_score_position_ptr) {

	int score = *best_score_ptr;
	int best_score = *best_score_ptr;
	int best_score_position = -1;
	const packed_alphabet_code::PackedAlphabetCode * sequence0_code_block_cache =
			TDirection == cuda_common::kReverse ?
					sequence0_code_block_cache_mem + cache_length - 1 :
					sequence0_code_block_cache_mem;
	const packed_alphabet_code::PackedAlphabetCode * sequence1_code_block_cache =
			TDirection == cuda_common::kReverse ?
					sequence1_code_block_cache_mem + cache_length - 1 :
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
	packed_alphabet_code::InitPackedAlphabetCode<TDirection>(
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
	packed_alphabet_code::InitPackedAlphabetCode<TDirection>(
			sequence1_code_block_cache, s1_offset,
			&s1_next_bolock_chache_position, &s1_block_sift,
			&s1_next_block_sift, &s1_temp_code_block);
#if 0
	if ((debug_db_p - 1) == s1_start_position || debug_db_p == s1_start_position
			|| 0 == blockIdx.x * blockDim.x + threadIdx.x) {
		printf("temp code block %ld, %ld\n", s0_temp_code_block,
				s1_temp_code_block);
	}
#endif

	int ret = 1;
	int cache_end = (cache_length - 1)
			* (TDirection == cuda_common::kReverse ? -1 : 1);
	int sequence_offset = 0;
	while (ret && (s0_next_bolock_chache_position != cache_end)) {
		packed_alphabet_code::PackedAlphabetCode s0_code_block =
				packed_alphabet_code::GetPackedAlphabetCode<TDirection>(
						sequence0_code_block_cache, s0_block_sift,
						s0_next_block_sift, &s0_next_bolock_chache_position,
						&s0_temp_code_block);

		packed_alphabet_code::PackedAlphabetCode s1_code_block =
				packed_alphabet_code::GetPackedAlphabetCode<TDirection>(
						sequence1_code_block_cache, s1_block_sift,
						s1_next_block_sift, &s1_next_bolock_chache_position,
						&s1_temp_code_block);

#if 0
		if ((debug_db_p - 1) == s1_start_position
				|| debug_db_p == s1_start_position
				|| 0 == blockIdx.x * blockDim.x + threadIdx.x) {
			printf("code block %ld, %ld\n", s0_code_block, s1_code_block);
		}
#endif
		ret = !UngappedExtendOneSideMultiAlphabetCode<TDirection>(s0_code_block,
				s1_code_block, sequence_delimiter, score_matrix, number_letters,
				cutoff, sequence_offset, &score, &best_score,
				&best_score_position);
		sequence_offset += packed_alphabet_code::kNumberOfCodesInBlock;
	}

	s0_next_bolock_position += s0_next_bolock_chache_position;
	s1_next_bolock_position += s1_next_bolock_chache_position;
	while (ret) {
		packed_alphabet_code::PackedAlphabetCode s0_code_block =
				packed_alphabet_code::GetPackedAlphabetCode<TDirection>(
						sequence0_code_block, s0_block_sift, s0_next_block_sift,
						&s0_next_bolock_position, &s0_temp_code_block);

		packed_alphabet_code::PackedAlphabetCode s1_code_block =
				packed_alphabet_code::GetPackedAlphabetCode<TDirection>(
						sequence1_code_block, s1_block_sift, s1_next_block_sift,
						&s1_next_bolock_position, &s1_temp_code_block);

#if 0
		if ((debug_db_p - 1) == s1_start_position
				|| debug_db_p == s1_start_position
				|| 0 == blockIdx.x * blockDim.x + threadIdx.x) {
			printf("code_block %ld, %ld\n", s0_code_block, s1_code_block);
		}
#endif
		ret = !UngappedExtendOneSideMultiAlphabetCode<TDirection>(s0_code_block,
				s1_code_block, sequence_delimiter, score_matrix, number_letters,
				cutoff, sequence_offset, &score, &best_score,
				&best_score_position);
		sequence_offset += packed_alphabet_code::kNumberOfCodesInBlock;
	}
	*best_score_ptr = best_score;
	*best_score_position_ptr = best_score_position;
	return 0;
}

__global__ void
__launch_bounds__(ungapped_extention_with_trigger_gpu_kernel::kNumberThreads, 1) ExtendOneSideKernel(
		const packed_alphabet_code::PackedAlphabetCode* sequence0_code_block,
		const packed_alphabet_code::PackedAlphabetCode* sequence1_code_block,
		const uint32_t number_extensions, const bool reverse,
		const uint32_t* /* __restrict__  */query_ids,
		const AlphabetCoder::Code sequence_delimiter,
		const int* /* __restrict__  */score_matrix,
		const uint32_t number_letters, const int* /* __restrict__  */cutoffs,
		uint32_t * sequence0_positions, uint32_t* sequence1_positions,
		int* /* __restrict__  */best_scores) {

	extern __shared__ char shared_mem[];
	uint32_t* s_sequence_positions = (uint32_t *) shared_mem;
	uint64_t* s_group_loader =
			(packed_alphabet_code::PackedAlphabetCode *) &s_sequence_positions[blockDim.x];
	GroupLoader<const packed_alphabet_code::PackedAlphabetCode *, kLoadLength> group_loader;
	int group_idx = group_loader.GetGroupId();
	int idx_in_group = group_loader.GetIdInGroup();
	int number_members = group_loader.GetNumberMembers();
	int group_work_skip = group_loader.GetNumberGroups();
	int number_group_in_block = blockDim.x / number_members;
	int number_group_works = (number_extensions + number_members - 1)
			/ number_members;
	int remain_group_in_block = number_group_works % number_group_in_block;
	int group_work_end = number_group_works
			+ (remain_group_in_block != 0 ?
					(number_group_in_block - remain_group_in_block) : 0);
	uint32_t* s_group_sequence_positions = &s_sequence_positions[(threadIdx.x
			/ number_members) * number_members];

	if (reverse) {
		// reverse
		for (int group_work_i = group_idx; group_work_i < group_work_end;
				group_work_i += group_work_skip) {

			const int group_extension_begin = group_work_i * number_members;
			const int thread_work_id = group_extension_begin + idx_in_group;
			const int extension_id =
					thread_work_id < number_extensions ? thread_work_id : 0;
			packed_alphabet_code::PackedAlphabetCode s0_cache[kLoadLength];
			packed_alphabet_code::PackedAlphabetCode s1_cache[kLoadLength];

			const uint32_t sequence0_position =
					sequence0_positions[extension_id]
							+ cuda_common::kMaxLoadLength
									* packed_alphabet_code::kNumberOfCodesInBlock;

			const uint32_t sequence1_position =
					sequence1_positions[extension_id]
							+ cuda_common::kMaxLoadLength
									* packed_alphabet_code::kNumberOfCodesInBlock;

			uint32_t sequence0_brock_position =
					packed_alphabet_code::GetBlockPosition(sequence0_position)
							- kLoadLength + 1;
			s_group_sequence_positions[idx_in_group] = sequence0_brock_position;
			//__syncthreads();
			group_loader.Load(sequence0_code_block, s_group_sequence_positions,
					min((int) number_extensions - group_extension_begin,
							number_members), s_group_loader, s0_cache);
#if 0
			// debug /////////////////////////////////////
			if (debug_db_p == sequence1_position
					|| 0 == blockIdx.x * blockDim.x + threadIdx.x) {
				printf("block position : %d\n",
						s_group_sequence_positions[idx_in_group]);
				for (int i = 0; i < kLoadLength; ++i) {
					printf("%d : %ld\n", i, s0_cache[i]);
				}
			}
			/////////////////////////////////////////////
#endif
			uint32_t sequence1_brock_position =
					packed_alphabet_code::GetBlockPosition(sequence1_position)
							- kLoadLength + 1;
			s_group_sequence_positions[idx_in_group] = sequence1_brock_position;
			//__syncthreads();
			group_loader.Load(sequence1_code_block, s_group_sequence_positions,
					min((int) number_extensions - group_extension_begin,
							number_members), s_group_loader, s1_cache);
#if 0
			// debug /////////////////////////////////////
			if (debug_db_p == sequence1_position
					|| 0 == blockIdx.x * blockDim.x + threadIdx.x) {
				printf("block position : %d\n",
						s_group_sequence_positions[idx_in_group]);
				for (int i = 0; i < kLoadLength; ++i) {
					printf("%d : %ld\n", i, s1_cache[i]);
				}
			}
			/////////////////////////////////////////////
#endif

			uint32_t query_id = 0;
			int cutoff = 0;
			if (thread_work_id < number_extensions) {
				query_id = query_ids[extension_id];
				cutoff = cutoffs[query_id];
				int best_score = 0;
				int best_score_position = 0;
				UngappedExtendOneSideDevice<cuda_common::kReverse>(
						sequence0_code_block, sequence0_position,
						sequence1_code_block, sequence1_position, s0_cache,
						s1_cache, kLoadLength, sequence_delimiter, score_matrix,
						number_letters, cutoff, &best_score,
						&best_score_position);

				best_scores[extension_id] = best_score;
				sequence0_positions[extension_id] = sequence0_position
						- best_score_position
						- cuda_common::kMaxLoadLength
								* packed_alphabet_code::kNumberOfCodesInBlock;
				sequence1_positions[extension_id] = sequence1_position
						- best_score_position
						- cuda_common::kMaxLoadLength
								* packed_alphabet_code::kNumberOfCodesInBlock;
			}
		}
	} else {
		// foward
		for (int group_work_i = group_idx; group_work_i < group_work_end;
				group_work_i += group_work_skip) {
			const int group_extension_begin = group_work_i * number_members;
			const int thread_work_id = group_extension_begin + idx_in_group;
			const int extension_id =
					thread_work_id < number_extensions ? thread_work_id : 0;
			packed_alphabet_code::PackedAlphabetCode s0_cache[kLoadLength];
			packed_alphabet_code::PackedAlphabetCode s1_cache[kLoadLength];

			const uint32_t sequence0_position =
					sequence0_positions[extension_id]
							+ cuda_common::kMaxLoadLength
									* packed_alphabet_code::kNumberOfCodesInBlock;

			const uint32_t sequence1_position =
					sequence1_positions[extension_id]
							+ cuda_common::kMaxLoadLength
									* packed_alphabet_code::kNumberOfCodesInBlock;

			uint32_t sequence0_brock_position =
					packed_alphabet_code::GetBlockPosition(sequence0_position);
			s_group_sequence_positions[idx_in_group] = sequence0_brock_position;
			//__syncthreads();
			group_loader.Load(sequence0_code_block, s_group_sequence_positions,
					min((int) number_extensions - group_extension_begin,
							number_members), s_group_loader, s0_cache);

			uint32_t sequence1_brock_position =
					packed_alphabet_code::GetBlockPosition(sequence1_position);
			s_group_sequence_positions[idx_in_group] = sequence1_brock_position;
			//__syncthreads();
			group_loader.Load(sequence1_code_block, s_group_sequence_positions,
					min((int) number_extensions - group_extension_begin,
							number_members), s_group_loader, s1_cache);
			uint32_t query_id = 0;
			int cutoff = 0;
			if (thread_work_id < number_extensions) {
				query_id = query_ids[extension_id];
				cutoff = cutoffs[query_id];
				int best_score = 0;
				int best_score_position = 0;
				UngappedExtendOneSideDevice<cuda_common::kFoward>(
						sequence0_code_block, sequence0_position,
						sequence1_code_block, sequence1_position, s0_cache,
						s1_cache, kLoadLength, sequence_delimiter, score_matrix,
						number_letters, cutoff, &best_score,
						&best_score_position);

				best_scores[extension_id] = best_score;
				sequence0_positions[extension_id] = sequence0_position
						+ best_score_position
						- cuda_common::kMaxLoadLength
								* packed_alphabet_code::kNumberOfCodesInBlock;
				sequence1_positions[extension_id] = sequence1_position
						+ best_score_position
						- cuda_common::kMaxLoadLength
								* packed_alphabet_code::kNumberOfCodesInBlock;

			}

		}
	}
	return;
}

}

UngappedExtenderGpu::UngappedExtenderGpu() :
		number_letters_(0), sequence_delimiter_(0), d_database_sequence_(NULL), d_concatenated_query_sequence_(
				NULL), d_score_matrix_(NULL), d_ungapped_extension_cutoffs_(
				NULL), d_gapped_extension_triggers_(NULL) {
}

UngappedExtenderGpu::~UngappedExtenderGpu() {

}

int UngappedExtenderGpu::SetQuerySeedDataList(const uint32_t* d_query_ids,
		const uint32_t* d_query_positions) {
	d_query_seed_query_ids_ = d_query_ids;
	d_query_seed_query_positions_ = d_query_positions;
	return 0;
}

int UngappedExtenderGpu::SetQueries(AlphabetCoder::Code sequence_delimiter,
		const packed_alphabet_code::PackedAlphabetCode *d_concatenated_query_sequence,
		const int *d_ungapped_extension_cutoffs,
		const int *d_gapped_extension_triggers) {
	sequence_delimiter_ = sequence_delimiter;
	d_concatenated_query_sequence_ = d_concatenated_query_sequence;
	d_ungapped_extension_cutoffs_ = d_ungapped_extension_cutoffs;
	d_gapped_extension_triggers_ = d_gapped_extension_triggers;
	return 0;
}
int UngappedExtenderGpu::SetDatabase(
		const packed_alphabet_code::PackedAlphabetCode *d_database_sequence) {
	d_database_sequence_ = d_database_sequence;
	return 0;
}
int UngappedExtenderGpu::SetScoreMatrix(const int *d_score_matrix,
		uint32_t number_letters) {
	number_letters_ = number_letters;
	d_score_matrix_ = d_score_matrix;
	return 0;
}

int UngappedExtenderGpu::SetQuerySeeds(size_t query_seeds_size, size_t size,
		uint32_t* query_seed_ids, uint32_t* query_seed_starts,
		size_t temp_storage_bytes, void* d_temp_storage,
		uint32_t* d_query_seed_ids, uint32_t* d_query_seed_starts,
		uint32_t* d_query_ids, uint32_t* d_query_concatenated_positions,
		uint32_t* temp_array, cudaStream_t &stream) const {
	CUDA_CHECK_RETURN(
			cudaMemcpyAsync(d_query_seed_starts, query_seed_starts,
					sizeof(query_seed_starts[0]) * query_seeds_size,
					cudaMemcpyDefault, stream));

	uint32_t* d_set_query_seed_ids_i_start = d_query_ids;
	CUDA_CHECK_RETURN(
			cudaMemsetAsync(d_set_query_seed_ids_i_start, 0,
					sizeof(uint32_t) * size, stream));

	ungapped_extention_with_trigger_gpu_kernel::SetQuerySeedIdsIAtStart<<<256, 64, 0, stream>>>(
			query_seeds_size, d_query_seed_starts, size, d_set_query_seed_ids_i_start);

	uint32_t* d_set_query_seed_id_i_list = temp_array;
	CUDA_CHECK_RETURN(
			cub::DeviceScan::InclusiveScan(d_temp_storage, temp_storage_bytes,
					d_set_query_seed_ids_i_start, d_set_query_seed_id_i_list,
					cub::Max(), size, stream));

	CUDA_CHECK_RETURN(
			cudaMemcpyAsync(d_query_seed_ids, query_seed_ids,
					sizeof(query_seed_ids[0]) * query_seeds_size,
					cudaMemcpyDefault, stream));

	ungapped_extention_with_trigger_gpu_kernel::SetQuerySeeds<<<256, 64, 0, stream>>>(
			size, d_query_seed_query_ids_, d_query_seed_query_positions_, d_query_seed_ids, d_set_query_seed_id_i_list, d_query_ids, d_query_concatenated_positions);

	return 0;
}

int UngappedExtenderGpu::ExtendWithTriggerAsync(size_t size,
		uint32_t* query_ids, uint32_t* query_concatenated_positions,
		uint32_t* database_positions, char* flags, uint32_t* d_query_ids,
		uint32_t* d_query_concatenated_positions,
		uint32_t* d_database_positions, char* d_flags, int* d_temp_array,
		cudaStream_t &stream) const {

	if (query_ids != NULL) {
		CUDA_CHECK_RETURN(
				cudaMemcpyAsync(d_query_ids, query_ids,
						sizeof(query_ids[0]) * size, cudaMemcpyDefault, stream));
	}
	if (query_concatenated_positions != NULL) {
		CUDA_CHECK_RETURN(
				cudaMemcpyAsync(d_query_concatenated_positions,
						query_concatenated_positions,
						sizeof(query_concatenated_positions[0]) * size,
						cudaMemcpyDefault, stream));
	}

	CUDA_CHECK_RETURN(
			cudaMemcpyAsync(d_database_positions, database_positions,
					sizeof(database_positions[0]) * size, cudaMemcpyDefault,
					stream));

	size_t shared_memory_size =
			sizeof(uint32_t)
					* ungapped_extention_with_trigger_gpu_kernel::kNumberThreads
					+ GroupLoader<
							const packed_alphabet_code::PackedAlphabetCode *,
							ungapped_extention_with_trigger_gpu_kernel::kLoadLength>::GetTotalSharedMemorySize(
							ungapped_extention_with_trigger_gpu_kernel::kNumberThreads);
	ungapped_extention_with_trigger_gpu_kernel::UngappedExtendKernel<<<
	ungapped_extention_with_trigger_gpu_kernel::kNumberBlocks,
	ungapped_extention_with_trigger_gpu_kernel::kNumberThreads, shared_memory_size,
	stream>>>(d_concatenated_query_sequence_, d_database_sequence_,
			size, d_query_ids, d_query_concatenated_positions,
			d_database_positions, sequence_delimiter_, d_score_matrix_,
			number_letters_ + 1, d_ungapped_extension_cutoffs_,
			d_gapped_extension_triggers_, d_flags);

	CUDA_CHECK_RETURN(
			cudaMemcpyAsync(flags, d_flags, sizeof(flags[0]) * size,
					cudaMemcpyDefault, stream));
	return 0;
}

int UngappedExtenderGpu::ExtendOneSideAsync(size_t size, bool reverse,
		uint32_t* query_ids, uint32_t* query_concatenated_positions,
		uint32_t* database_positions, int* scores, uint32_t* d_query_ids,
		uint32_t* d_query_concatenated_positions,
		uint32_t* d_database_positions, int* d_scores,
		cudaStream_t &stream) const {

	CUDA_CHECK_RETURN(
			cudaMemcpyAsync(d_query_ids, query_ids, sizeof(query_ids[0]) * size,
					cudaMemcpyDefault, stream));

	CUDA_CHECK_RETURN(
			cudaMemcpyAsync(d_query_concatenated_positions,
					query_concatenated_positions,
					sizeof(query_concatenated_positions[0]) * size,
					cudaMemcpyDefault, stream));

	CUDA_CHECK_RETURN(
			cudaMemcpyAsync(d_database_positions, database_positions,
					sizeof(database_positions[0]) * size, cudaMemcpyDefault,
					stream));

	size_t shared_memory_size =
			sizeof(uint32_t) * ungapped_extention_gpu_kernel::kNumberThreads
					+ GroupLoader<
							const packed_alphabet_code::PackedAlphabetCode *,
							ungapped_extention_gpu_kernel::kLoadLength>::GetTotalSharedMemorySize(
							ungapped_extention_gpu_kernel::kNumberThreads);
	ungapped_extention_gpu_kernel::ExtendOneSideKernel<<<
	ungapped_extention_gpu_kernel::kNumberBlocks,
	ungapped_extention_gpu_kernel::kNumberThreads, shared_memory_size,
	stream>>>(d_concatenated_query_sequence_, d_database_sequence_,
			size, reverse, d_query_ids, sequence_delimiter_, d_score_matrix_,
			number_letters_ + 1, d_ungapped_extension_cutoffs_,
			d_query_concatenated_positions, d_database_positions, d_scores);

	CUDA_CHECK_RETURN(
			cudaMemcpyAsync(query_concatenated_positions,
					d_query_concatenated_positions,
					sizeof(query_concatenated_positions[0]) * size,
					cudaMemcpyDefault, stream));
	CUDA_CHECK_RETURN(
			cudaMemcpyAsync(database_positions, d_database_positions,
					sizeof(database_positions[0]) * size, cudaMemcpyDefault,
					stream));
	CUDA_CHECK_RETURN(
			cudaMemcpyAsync(scores, d_scores, sizeof(scores[0]) * size,
					cudaMemcpyDefault, stream));

	return 0;
}

