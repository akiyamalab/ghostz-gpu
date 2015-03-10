/*
 * gapped_extender_gpu.cpp
 *
 *  Created on: Aug 11, 2014
 *      Author: shu
 */

#include "gapped_extender_gpu.h"

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/system/cuda/experimental/pinned_allocator.h>
#include <thrust/copy.h>
#include <vector>
#include <assert.h>
#include <iostream>

#include "packed_alphabet_code.h"
#include "group_loader.h"
#include "score_matrix.h"
#include "cuda_common.h"

using namespace std;

namespace gapped_extender_gpu_kernel {
#if 0
const int debug_q_p = 119;
const int debug_db_p = 19959
+ cuda_common::kMaxLoadLength
* packed_alphabet_code::kNumberOfCodesInBlock;
const int debug_thread_id = 0;
#endif
static const size_t kNumberBlocks = 128;
static const size_t kNumberThreads = 128;
static const size_t kLoadLength = 4;
static const size_t kDpRowLength = packed_alphabet_code::kNumberOfCodesInBlock
		/ 3 + 1;

template<cuda_common::Direction TDirection>
__device__ int InitSequence(
		const packed_alphabet_code::PackedAlphabetCode* sequence_cache_mem,
		const packed_alphabet_code::PackedAlphabetCode* sequence_code_block,
		const int sequence_position, const uint32_t sequence_delimiter,
		AlphabetCoder::Code* cached_sequence) {

	const packed_alphabet_code::PackedAlphabetCode * sequence_code_block_cache =
			TDirection == cuda_common::kReverse ?
					sequence_cache_mem + kLoadLength - 1 : sequence_cache_mem;
	int next_bolock_position = packed_alphabet_code::GetBlockPosition(
			sequence_position);
	int offset = sequence_position
			- next_bolock_position
					* packed_alphabet_code::kNumberOfCodesInBlock;
	int next_bolock_chache_position = 0;
	uint32_t block_sift;
	uint32_t next_block_sift;
	packed_alphabet_code::PackedAlphabetCode temp_code_block;
	packed_alphabet_code::InitPackedAlphabetCode<TDirection>(
			sequence_code_block_cache, offset, &next_bolock_chache_position,
			&block_sift, &next_block_sift, &temp_code_block);
#if 0
	if ((debug_q_p - 1) == sequence_position || debug_q_p == sequence_position
			|| blockIdx.x * blockDim.x + threadIdx.x == 0) {
		printf("temp code block %ld\n", temp_code_block);
	}
#endif

	int cache_end = ((int) kLoadLength)
			* (TDirection == cuda_common::kReverse ? -1 : 1);
	uint32_t stop_flag = 0;
	int p = 0;
	int increment = TDirection == cuda_common::kReverse ? -1 : 1;
	while (!stop_flag && (next_bolock_chache_position != cache_end)) {
		packed_alphabet_code::PackedAlphabetCode code_block =
				packed_alphabet_code::GetPackedAlphabetCode<TDirection>(
						sequence_code_block_cache, block_sift, next_block_sift,
						&next_bolock_chache_position, &temp_code_block);
#if 0
		if ((debug_q_p - 1) == sequence_position
				|| debug_q_p == sequence_position
				|| blockIdx.x * blockDim.x + threadIdx.x == 0) {
			printf("code block %ld\n", code_block);
		}
#endif
#pragma unroll
		for (uint32_t i = 0; i < packed_alphabet_code::kNumberOfCodesInBlock;
				++i) {
			const uint32_t c =
					packed_alphabet_code::GetAlphabetCode<TDirection>(
							code_block, i);
#if 0
			if ((blockDim.x * blockIdx.x + threadIdx.x) == 0) {
				printf("%d : %d\n", i, c);
			}
#endif
			stop_flag = stop_flag | c == sequence_delimiter;
			cached_sequence[p] = c;
			p += increment;
		}
	}
	next_bolock_position += next_bolock_chache_position;
	while (!stop_flag) {
		packed_alphabet_code::PackedAlphabetCode code_block =
				packed_alphabet_code::GetPackedAlphabetCode<TDirection>(
						sequence_code_block, block_sift, next_block_sift,
						&next_bolock_position, &temp_code_block);
#pragma unroll
		for (uint32_t i = 0; i < packed_alphabet_code::kNumberOfCodesInBlock;
				++i) {
			const uint32_t c =
					packed_alphabet_code::GetAlphabetCode<TDirection>(
							code_block, i);
			stop_flag = stop_flag | c == sequence_delimiter;
			cached_sequence[p] = c;
			p += increment;
		}
	}

	return 0;
}

__device__ int InitDpColumn(const AlphabetCoder::Code* sequence,
		const AlphabetCoder::Code sequence_delimiter, int gap_init,
		int gap_extention, int default_cutoff, int increment,
		GappedExtenderGpu::DpColumnCell* dp_column) {
	int score = -gap_init;
	GappedExtenderGpu::DpColumnCell new_cell;
	new_cell.h = 0;
	new_cell.e = -gap_init;
	reinterpret_cast<uint64_t *>(dp_column)[0] =
			reinterpret_cast<uint64_t *>(&new_cell)[0];
	int sequence_position = 0;
	int array_index = 0;
	int cutoff = default_cutoff;
	if (cutoff < gap_init) {
		cutoff = gap_init;
	}
	for (array_index = 1; sequence[sequence_position] != sequence_delimiter;
			++array_index, sequence_position += increment) {
		if (score < -cutoff) {
			break;
		}
		new_cell.h = score;
		new_cell.e = score - gap_init;
		reinterpret_cast<uint64_t *>(dp_column)[array_index] =
				reinterpret_cast<uint64_t *>(&new_cell)[0];
		score -= gap_extention;
	}

	return array_index;
}

template<cuda_common::Direction TDirection>
__device__ uint32_t InitDpRow(
		const packed_alphabet_code::PackedAlphabetCode seq1_code_block,
		const uint32_t seq1_offset,
		const AlphabetCoder::Code sequence_delimiter, const int gap_init,
		const int gap_extention, const int default_start_dp_h,
		const int default_start_dp_e, const int threshold,
		AlphabetCoder::Code *seq_1_subsequence, int *dp_h, int *dp_f,
		int *last_dp_e) {

	int start_dp_h = default_start_dp_h;
	int start_dp_e = default_start_dp_e;
	if (threshold > default_start_dp_h) {
		start_dp_h = GappedExtenderGpu::kInitScore;
		start_dp_e = GappedExtenderGpu::kInitScore;
	}
	dp_h[0] = start_dp_h;
	dp_f[0] = start_dp_h - gap_init;
	int dp_e = start_dp_e;
	const int s_dp_row_skip = blockDim.x;
	int s_dp_row_i = s_dp_row_skip;
#pragma unroll
	for (uint32_t i = 1; i < kDpRowLength; ++i, s_dp_row_i += s_dp_row_skip) {
		if (threshold > dp_e) {
			dp_e = GappedExtenderGpu::kInitScore;
		}
		dp_h[s_dp_row_i] = dp_e;
		dp_f[s_dp_row_i] = dp_e - gap_init;
		dp_e -= gap_extention;
	}
	*last_dp_e = dp_e;

	const uint32_t seq1_subsequence_i_end = seq1_offset + kDpRowLength - 1;
	uint32_t seq1_subsequence_i = seq1_offset;
	uint32_t seq_1_subsequence_length = 0;
	for (; seq1_subsequence_i < seq1_subsequence_i_end; ++seq1_subsequence_i) {
		const AlphabetCoder::Code c = packed_alphabet_code::GetAlphabetCode<
				TDirection>(seq1_code_block, seq1_subsequence_i);
		if (c == sequence_delimiter) {
			break;
		}
		seq_1_subsequence[seq_1_subsequence_length] = c;
		++seq_1_subsequence_length;
	}

#if 0
	if (blockIdx.x * blockDim.x + threadIdx.x == debug_thread_id) {
		printf("start dp_h %d, start dp_e %d\n", default_start_dp_h,
				default_start_dp_e);
		printf("sizeof(int) = %d\n", sizeof(int));
		printf("     ");
		printf("     ");
		for (int i = 0; i < seq_1_subsequence_length; ++i) {
			printf("%5d", seq_1_subsequence[i]);
		}
		printf("\n");
		printf("     ");
		for (int i = 0; i < kDpRowLength; ++i) {
			printf("%5d", dp_h[i * blockDim.x]);
		}
		printf("\n");
	}
#endif
	return seq_1_subsequence_length;
}

__device__ int UpdateDpRow(const AlphabetCoder::Code *seq_1_subsequence,
		const int sequence0_position, const int sequence1_offset,
		const uint32_t seq_1_subsequence_length, int increment, int cutoff,
		const int *score_matrix_row, const int gap_init,
		const int gap_extention, int* dp_row_h, int* dp_row_f, int start_dp_e,
		int* max_score_sequence0_position, int* max_score_sequence1_position,
		int* max_score_ptr, int* last_dp_e) {

	int max_score_sequence1_position_in_row = -1;
	int max_score = *max_score_ptr;
	int dp_prev_h = dp_row_h[0];
	int dp_e = start_dp_e;
	const int dp_row_length = seq_1_subsequence_length + 1;
	int start_update_cell = dp_row_length;
	const int s_dp_row_skip = blockDim.x;
	int s_dp_row_i = s_dp_row_skip;
	int dp_f = GappedExtenderGpu::kInitScore;

	for (int i = 1; i < dp_row_length; ++i, s_dp_row_i += s_dp_row_skip) {
		int score = dp_prev_h
				+ __ldg(&score_matrix_row[seq_1_subsequence[i - 1]]);
		dp_prev_h = dp_row_h[s_dp_row_i];
		dp_f = dp_row_f[s_dp_row_i];
		score = max(score, dp_f);
		score = max(score, dp_e);
		if (score > max_score) {
			max_score = score;
			max_score_sequence1_position_in_row = i;
		}

		if (max_score - score > cutoff) {
			score = GappedExtenderGpu::kInitScore;
			dp_f = GappedExtenderGpu::kInitScore;
			dp_e = GappedExtenderGpu::kInitScore;
		} else {
			start_update_cell = min(start_update_cell, i);
		}

		dp_row_h[s_dp_row_i] = score;
		dp_row_f[s_dp_row_i] = max(score - gap_init, dp_f - gap_extention);
		dp_e = max(score - gap_init, dp_e - gap_extention);
	}

	*last_dp_e = dp_e;
	if (max_score_sequence1_position_in_row >= 0) {
		*max_score_ptr = max_score;
		*max_score_sequence0_position = sequence0_position;
		*max_score_sequence1_position = sequence1_offset
				+ (max_score_sequence1_position_in_row - 1) * increment;
	}

	return start_update_cell == dp_row_length;
}

template<cuda_common::Direction TDirection>
__device__ int ExtendOneSideScoreOnlyDevice(
		const AlphabetCoder::Code* sequence0,
		const packed_alphabet_code::PackedAlphabetCode* sequence1_code_block,
		const uint32_t sequence1_offset,
		const packed_alphabet_code::PackedAlphabetCode* sequence1_cache_mem,
		const bool reverse, const AlphabetCoder::Code sequence_delimiter,
		const int* score_matrix, const uint32_t number_letters,
		const int gap_open, const int gap_extention, const int cutoff,
		char shared_mem[], uint32_t* best_sequence0_position,
		uint32_t* best_sequence1_position, int* best_score) {

	const packed_alphabet_code::PackedAlphabetCode * sequence1_code_block_cache =
			reverse ?
					sequence1_cache_mem + kLoadLength - 1 : sequence1_cache_mem;
	int next_bolock_position = packed_alphabet_code::GetBlockPosition(
			sequence1_offset);
	int sequence1_code_block_offset = sequence1_offset
			- next_bolock_position
					* packed_alphabet_code::kNumberOfCodesInBlock;
	int next_bolock_chache_position = 0;
	uint32_t block_sift;
	uint32_t next_block_sift;
	packed_alphabet_code::PackedAlphabetCode temp_code_block;
	packed_alphabet_code::InitPackedAlphabetCode<TDirection>(
			sequence1_code_block_cache, sequence1_code_block_offset,
			&next_bolock_chache_position, &block_sift, &next_block_sift,
			&temp_code_block);

	int increment = reverse ? -1 : 1;
	GappedExtenderGpu::DpColumnCell dp_column[GappedExtenderGpu::kMaxSequence0Length];

	int sequence0_position_start = 0;
	int dp_column_end = 0;
	int gap_init = gap_open + gap_extention;
	dp_column_end = InitDpColumn(sequence0, sequence_delimiter, gap_init,
			gap_extention, cutoff, increment, dp_column);

	int cache_end = ((int) kLoadLength)
			* (TDirection == cuda_common::kReverse ? -1 : 1);
	uint32_t stop_flag = 0;
	int max_score = 0;
	int max_score_sequence0_position = -increment;
	int max_score_sequence1_position = -increment;
	int sequence1_position = 0;

	AlphabetCoder::Code seq_1_subsequence[kDpRowLength];
	int *dp_row_mem = (int*) shared_mem;
	int *dp_row_h = &dp_row_mem[threadIdx.x];
	int *dp_row_f = &dp_row_mem[blockDim.x * kDpRowLength + threadIdx.x];
	uint32_t seq_1_subsequence_offset =
			packed_alphabet_code::kNumberOfCodesInBlock;
	packed_alphabet_code::PackedAlphabetCode code_block = 0;

	while (!stop_flag && (next_bolock_chache_position != cache_end)) {
		if (seq_1_subsequence_offset
				>= packed_alphabet_code::kNumberOfCodesInBlock) {
			code_block =
					packed_alphabet_code::GetPackedAlphabetCode<TDirection>(
							sequence1_code_block_cache, block_sift,
							next_block_sift, &next_bolock_chache_position,
							&temp_code_block);
			seq_1_subsequence_offset = 0;
		}
#if 0
		if (blockIdx.x * blockDim.x + threadIdx.x == debug_thread_id) {
			printf("code block %ld \n", code_block);
		}

		if (blockIdx.x * blockDim.x + threadIdx.x == debug_thread_id) {
			printf("column %d:%2d\n", 0, dp_column_end);
		}
#endif
		int last_dp_e = 0;
		GappedExtenderGpu::DpColumnCell dp_column_cell;
		reinterpret_cast<uint64_t *>(&dp_column_cell)[0] =
				reinterpret_cast<uint64_t *>(dp_column)[0];
		uint32_t seq_1_subsequence_length = InitDpRow<TDirection>(code_block,
				seq_1_subsequence_offset, sequence_delimiter, gap_init,
				gap_extention, dp_column_cell.h, dp_column_cell.e,
				max_score - cutoff, seq_1_subsequence, dp_row_h, dp_row_f,
				&last_dp_e);
		seq_1_subsequence_offset += seq_1_subsequence_length;
		int last_s_dp_row_i = seq_1_subsequence_length * blockDim.x;
		dp_column_cell.h = dp_row_h[last_s_dp_row_i];
		dp_column_cell.e = last_dp_e;
		reinterpret_cast<uint64_t *>(dp_column)[0] =
				reinterpret_cast<uint64_t *>(&dp_column_cell)[0];
		int start_drop_flag = 1;
		int dp_column_stored_i = 1;
		int last_updated_dp_column_stored_i = -1;
		int ret = 0;
		int sequence0_position = sequence0_position_start;
		AlphabetCoder::Code s0_c;
		for (int column_i = 1; column_i < dp_column_end | !ret;
				++column_i, sequence0_position += increment) {
			s0_c = sequence0[sequence0_position];
			if (s0_c == sequence_delimiter) {
				break;
			}
			reinterpret_cast<uint64_t *>(&dp_column_cell)[0] =
					reinterpret_cast<uint64_t *>(dp_column)[column_i];
			int start_dp_e =
					column_i < dp_column_end ?
							dp_column_cell.e : GappedExtenderGpu::kInitScore;
			const int *score_matrix_row = &score_matrix[s0_c * number_letters];
#if 0
			if (blockIdx.x * blockDim.x + threadIdx.x == debug_thread_id) {
				printf("%2d:%2d", column_i, s0_c);
			}
#endif
			last_dp_e = 0;

			ret = UpdateDpRow(seq_1_subsequence, sequence0_position,
					sequence1_position, seq_1_subsequence_length, increment,
					cutoff, score_matrix_row, gap_init, gap_extention, dp_row_h,
					dp_row_f, start_dp_e, &max_score_sequence0_position,
					&max_score_sequence1_position, &max_score, &last_dp_e);

			if (start_drop_flag
					& dp_row_h[last_s_dp_row_i]
							== GappedExtenderGpu::kInitScore) {
				dp_column_stored_i = 0;
				sequence0_position_start += increment;
			} else {
				start_drop_flag = 0;
			}
			int temp_score =
					column_i < dp_column_end ?
							dp_column_cell.h : GappedExtenderGpu::kInitScore;

			dp_column_cell.h =
					ret ? GappedExtenderGpu::kInitScore : dp_row_h[last_s_dp_row_i];
			dp_column_cell.e = ret ? GappedExtenderGpu::kInitScore : last_dp_e;
			reinterpret_cast<uint64_t *>(dp_column)[dp_column_stored_i] =
					reinterpret_cast<uint64_t *>(&dp_column_cell)[0];
			++dp_column_stored_i;
			last_updated_dp_column_stored_i =
					ret ? last_updated_dp_column_stored_i : dp_column_stored_i;
			dp_row_h[0] = temp_score;
#if 0
			if (blockIdx.x * blockDim.x + threadIdx.x == debug_thread_id) {
				for (int i = 0; i <= seq_1_subsequence_length; ++i) {
					printf("%5d", dp_row_h[i * blockDim.x]);
				}
				printf(" (ret = %d, dp_column_h[column_i] = %d)", ret,
						dp_column_h[column_i]);
				printf("\n");
			}
#endif
		}
		if (seq_1_subsequence_length < (kDpRowLength - 1)
				|| 1 >= dp_column_stored_i) {
			stop_flag = true;
			break;
		}
		dp_column_end = last_updated_dp_column_stored_i + 1;
		sequence1_position += seq_1_subsequence_length * increment;
	}

	next_bolock_position += next_bolock_chache_position;
	while (!stop_flag) {
		if (seq_1_subsequence_offset
				>= packed_alphabet_code::kNumberOfCodesInBlock) {
			code_block =
					packed_alphabet_code::GetPackedAlphabetCode<TDirection>(
							sequence1_code_block, block_sift, next_block_sift,
							&next_bolock_position, &temp_code_block);
			seq_1_subsequence_offset = 0;
		}
#if 0
		if (blockIdx.x * blockDim.x + threadIdx.x == debug_thread_id) {
			printf("code block %ld \n", code_block);
		}

		if (blockIdx.x * blockDim.x + threadIdx.x == debug_thread_id) {
			printf("column %d:%2d\n", 0, dp_column_end);
		}
#endif
		int last_dp_e = 0;
		GappedExtenderGpu::DpColumnCell dp_column_cell;
		reinterpret_cast<uint64_t *>(&dp_column_cell)[0] =
				reinterpret_cast<uint64_t *>(dp_column)[0];
		uint32_t seq_1_subsequence_length = InitDpRow<TDirection>(code_block,
				seq_1_subsequence_offset, sequence_delimiter, gap_init,
				gap_extention, dp_column_cell.h, dp_column_cell.e,
				max_score - cutoff, seq_1_subsequence, dp_row_h, dp_row_f,
				&last_dp_e);
		seq_1_subsequence_offset += seq_1_subsequence_length;
		int last_s_dp_row_i = seq_1_subsequence_length * blockDim.x;
		dp_column_cell.h = dp_row_h[last_s_dp_row_i];
		dp_column_cell.e = last_dp_e;
		reinterpret_cast<uint64_t *>(dp_column)[0] =
				reinterpret_cast<uint64_t *>(&dp_column_cell)[0];
		int start_drop_flag = 1;
		int dp_column_stored_i = 1;
		int last_updated_dp_column_stored_i = -1;
		int ret = 0;
		int sequence0_position = sequence0_position_start;
		AlphabetCoder::Code s0_c;
		for (int column_i = 1; column_i < dp_column_end | !ret;
				++column_i, sequence0_position += increment) {
			s0_c = sequence0[sequence0_position];
			if (s0_c == sequence_delimiter) {
				break;
			}
			reinterpret_cast<uint64_t *>(&dp_column_cell)[0] =
					reinterpret_cast<uint64_t *>(dp_column)[column_i];
			int start_dp_e =
					column_i < dp_column_end ?
							dp_column_cell.e : GappedExtenderGpu::kInitScore;
			const int *score_matrix_row = &score_matrix[s0_c * number_letters];
#if 0
			if (blockIdx.x * blockDim.x + threadIdx.x == debug_thread_id) {
				printf("%2d:%2d", column_i, s0_c);
			}
#endif
			last_dp_e = 0;

			ret = UpdateDpRow(seq_1_subsequence, sequence0_position,
					sequence1_position, seq_1_subsequence_length, increment,
					cutoff, score_matrix_row, gap_init, gap_extention, dp_row_h,
					dp_row_f, start_dp_e, &max_score_sequence0_position,
					&max_score_sequence1_position, &max_score, &last_dp_e);

			if (start_drop_flag
					& dp_row_h[last_s_dp_row_i]
							== GappedExtenderGpu::kInitScore) {
				dp_column_stored_i = 0;
				sequence0_position_start += increment;
			} else {
				start_drop_flag = 0;
			}
			int temp_score =
					column_i < dp_column_end ?
							dp_column_cell.h : GappedExtenderGpu::kInitScore;

			dp_column_cell.h =
					ret ? GappedExtenderGpu::kInitScore : dp_row_h[last_s_dp_row_i];
			dp_column_cell.e = ret ? GappedExtenderGpu::kInitScore : last_dp_e;
			reinterpret_cast<uint64_t *>(dp_column)[dp_column_stored_i] =
					reinterpret_cast<uint64_t *>(&dp_column_cell)[0];
			++dp_column_stored_i;
			last_updated_dp_column_stored_i =
					ret ? last_updated_dp_column_stored_i : dp_column_stored_i;
			dp_row_h[0] = temp_score;
#if 0
			if (blockIdx.x * blockDim.x + threadIdx.x == debug_thread_id) {
				for (int i = 0; i <= seq_1_subsequence_length; ++i) {
					printf("%5d", dp_row_h[i * blockDim.x]);
				}
				printf(" (ret = %d, dp_column_h[column_i] = %d)", ret,
						dp_column_h[column_i]);
				printf("\n");
			}
#endif
		}
		if (seq_1_subsequence_length < (kDpRowLength - 1)
				|| 1 >= dp_column_stored_i) {
			stop_flag = true;
			break;
		}
		dp_column_end = last_updated_dp_column_stored_i + 1;
		sequence1_position += seq_1_subsequence_length * increment;
	}

	*best_score = max_score;
	*best_sequence0_position = max_score_sequence0_position;
	*best_sequence1_position = max_score_sequence1_position;
	return 0;
}

__global__ void
__launch_bounds__(gapped_extender_gpu_kernel::kNumberThreads, 10) ExtendOneSideScoreOnlyKernel(
		const packed_alphabet_code::PackedAlphabetCode* sequence0_code_block,
		const packed_alphabet_code::PackedAlphabetCode* sequence1_code_block,
		const uint32_t number_extensions, const bool reverse,
		const AlphabetCoder::Code sequence_delimiter, const int* score_matrix,
		const uint32_t number_letters, const int gap_open,
		const int gap_extention, const int cutoff,
		uint32_t* sequence0_positions, uint32_t* sequence1_positions,
		int* best_scores) {

	extern __shared__ char shared_mem[];
	uint32_t* s_sequence_positions = (uint32_t *) shared_mem;
	packed_alphabet_code::PackedAlphabetCode* s_group_loader =
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
		packed_alphabet_code::PackedAlphabetCode sequence_cache[kLoadLength];

		const uint32_t sequence0_position = sequence0_positions[extension_id]
				+ cuda_common::kMaxLoadLength
						* packed_alphabet_code::kNumberOfCodesInBlock;
		const uint32_t sequence0_brock_position =
				packed_alphabet_code::GetBlockPosition(sequence0_position)
						+ (reverse ? -kLoadLength + 1 : 0);
		s_group_sequence_positions[idx_in_group] = sequence0_brock_position;
		//__syncthreads();
		group_loader.Load(sequence0_code_block, s_group_sequence_positions,
				min((int) number_extensions - group_extension_begin,
						number_members), s_group_loader, sequence_cache);
#if 0
		// debug /////////////////////////////////////
		if (debug_db_p == sequence0_position
				|| blockIdx.x * blockDim.x + threadIdx.x == debug_thread_id) {
			printf("block position : %d\n",
					s_group_sequence_positions[idx_in_group]);
			for (int i = 0; i < kLoadLength; ++i) {
				printf("%d : %ld\n", i, sequence_cache[i]);
			}
		}
		/////////////////////////////////////////////
#endif

		AlphabetCoder::Code cached_sequence0_mem[GappedExtenderGpu::kMaxSequence0Length];
		AlphabetCoder::Code *cached_sequence0 =
				reverse ?
						&cached_sequence0_mem[GappedExtenderGpu::kMaxSequence0Length
								- 1] :
						&cached_sequence0_mem[0];

		if (thread_work_id < number_extensions) {
			if (reverse) {
				InitSequence<cuda_common::kReverse>(sequence_cache,
						sequence0_code_block, sequence0_position,
						sequence_delimiter, cached_sequence0);
			} else {
				InitSequence<cuda_common::kFoward>(sequence_cache,
						sequence0_code_block, sequence0_position,
						sequence_delimiter, cached_sequence0);
			}
		}
		const uint32_t sequence1_position = sequence1_positions[extension_id]
				+ cuda_common::kMaxLoadLength
						* packed_alphabet_code::kNumberOfCodesInBlock;
		uint32_t sequence1_brock_position =
				packed_alphabet_code::GetBlockPosition(sequence1_position)
						+ (reverse ? -kLoadLength + 1 : 0);
		s_group_sequence_positions[idx_in_group] = sequence1_brock_position;
//__syncthreads();
		group_loader.Load(sequence1_code_block, s_group_sequence_positions,
				min((int) number_extensions - group_extension_begin,
						number_members), s_group_loader, sequence_cache);
#if 0
// debug /////////////////////////////////////
		if (debug_db_p == sequence1_position
				|| blockIdx.x * blockDim.x + threadIdx.x == debug_thread_id) {
			printf("thread id : %d\n", blockIdx.x * blockDim.x + threadIdx.x);
			printf("sequence1_position : %d\n", sequence1_position);
			printf("block position : %d\n",
					s_group_sequence_positions[idx_in_group]);
			for (int i = 0; i < kLoadLength; ++i) {
				printf("%d : %ld\n", i, sequence_cache[i]);
			}
		}
/////////////////////////////////////////////
#endif
		__syncthreads();
		if (thread_work_id < number_extensions) {
			uint32_t best_sequence_0_p = 0;
			uint32_t best_sequence_1_p = 0;
			int best_score = 0;
			if (reverse) {
				ExtendOneSideScoreOnlyDevice<cuda_common::kReverse>(
						cached_sequence0, sequence1_code_block,
						sequence1_position, sequence_cache, reverse,
						(uint32_t) sequence_delimiter, score_matrix,
						number_letters, gap_open, gap_extention, cutoff,
						shared_mem, &best_sequence_0_p, &best_sequence_1_p,
						&best_score);

				sequence0_positions[extension_id] = sequence0_position
						+ best_sequence_0_p
						- cuda_common::kMaxLoadLength
								* packed_alphabet_code::kNumberOfCodesInBlock;
				sequence1_positions[extension_id] = sequence1_position
						+ best_sequence_1_p
						- cuda_common::kMaxLoadLength
								* packed_alphabet_code::kNumberOfCodesInBlock;

				best_scores[extension_id] = best_score;
#if 0
// debug /////////////////////////////////////
				if (debug_db_p == sequence1_position
						|| blockIdx.x * blockDim.x + threadIdx.x
						== debug_thread_id) {
					printf("reverse gapped extend score : %d\n", best_score);
				}
/////////////////////////////////////////////
#endif
			} else {
				ExtendOneSideScoreOnlyDevice<cuda_common::kFoward>(
						cached_sequence0, sequence1_code_block,
						sequence1_position, sequence_cache, reverse,
						(uint32_t) sequence_delimiter, score_matrix,
						number_letters, gap_open, gap_extention, cutoff,
						shared_mem, &best_sequence_0_p, &best_sequence_1_p,
						&best_score);

				sequence0_positions[extension_id] = sequence0_position
						+ best_sequence_0_p
						- cuda_common::kMaxLoadLength
								* packed_alphabet_code::kNumberOfCodesInBlock;
				sequence1_positions[extension_id] = sequence1_position
						+ best_sequence_1_p
						- cuda_common::kMaxLoadLength
								* packed_alphabet_code::kNumberOfCodesInBlock;

				best_scores[extension_id] = best_score;
#if 0
// debug /////////////////////////////////////
				if (debug_db_p == sequence1_position
						|| blockIdx.x * blockDim.x + threadIdx.x
						== debug_thread_id) {
					printf("forward gapped extend score : %d\n", best_score);
				}
/////////////////////////////////////////////
#endif
			}
		}
		__syncthreads();
	}
	return;
}

__global__ void ConvertToGappedExtensionSeeds(const uint32_t size, bool reverse,
		const uint32_t *seed_positions,
		uint32_t* ungapped_extension_sequence_positions,
		uint32_t* gapped_extension_sequence_positions) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int skip = blockDim.x * gridDim.x;
	for (int i = idx; i < size; i += skip) {
		uint32_t seed_position = seed_positions[i];
		uint32_t sequence_position =
				ungapped_extension_sequence_positions[seed_position];
		gapped_extension_sequence_positions[i] =
				reverse ?
						sequence_position:
						sequence_position;
	}
}
}

GappedExtenderGpu::GappedExtenderGpu() :
		number_letters_(0), gap_open_(0), gap_extention_(0), cutoff_(0), sequence_delimiter_(
				0), d_database_sequence_(NULL), d_concatenated_query_sequence_(
				0), d_score_matrix_(0) {

}

GappedExtenderGpu::~GappedExtenderGpu() {

}

int GappedExtenderGpu::SetQueries(AlphabetCoder::Code sequence_delimiter,
		const packed_alphabet_code::PackedAlphabetCode *d_concatenated_query_sequence) {
	sequence_delimiter_ = sequence_delimiter;
	d_concatenated_query_sequence_ = d_concatenated_query_sequence;
	return 0;
}
int GappedExtenderGpu::SetDatabase(
		const packed_alphabet_code::PackedAlphabetCode *d_database_sequence) {
	d_database_sequence_ = d_database_sequence;
	return 0;
}
int GappedExtenderGpu::SetScoreParameters(const int *d_score_matrix,
		uint32_t number_letters, int gap_open, int gap_extention, int cutoff) {
	number_letters_ = number_letters;
	d_score_matrix_ = d_score_matrix;
	gap_open_ = gap_open;
	gap_extention_ = gap_extention;
	cutoff_ = cutoff;
	return 0;
}

int GappedExtenderGpu::ConvertToGappedExtensionSeedsAsync(size_t size,
		bool reverse, uint32_t* seed_positions,
		uint32_t* query_concatenated_positions, uint32_t* database_positions,
		uint32_t* d_seed_positions, uint32_t* d_temp_array,
		uint32_t* d_query_concatenated_positions,
		uint32_t* d_database_positions, cudaStream_t &stream) const {

	CUDA_CHECK_RETURN(
			cudaMemcpyAsync(d_seed_positions, seed_positions,
					sizeof(seed_positions[0]) * size, cudaMemcpyDefault,
					stream));

	CUDA_CHECK_RETURN(
			cudaMemcpyAsync(d_temp_array, query_concatenated_positions,
					sizeof(query_concatenated_positions[0]) * size,
					cudaMemcpyDefault, stream));

	gapped_extender_gpu_kernel::ConvertToGappedExtensionSeeds<<<128, 256, 0, stream>>>(
			size, reverse, d_seed_positions, d_temp_array, d_query_concatenated_positions
	);

	CUDA_CHECK_RETURN(
			cudaMemcpyAsync(d_temp_array, database_positions,
					sizeof(database_positions[0]) * size, cudaMemcpyDefault,
					stream));

	gapped_extender_gpu_kernel::ConvertToGappedExtensionSeeds<<<128, 256, 0, stream>>>(
			size, reverse, d_seed_positions, d_temp_array, d_database_positions
	);

	return 0;
}

int GappedExtenderGpu::ExtendOneSideScoreOnlyAsync(size_t size, bool reverse,
		uint32_t* query_concatenated_positions, uint32_t* database_positions,
		int* scores, uint32_t* d_query_concatenated_positions,
		uint32_t* d_database_positions, int* d_scores,
		cudaStream_t &stream) const {

	/*
	 CUDA_CHECK_RETURN(
	 cudaMemcpyAsync(d_query_concatenated_positions,
	 query_concatenated_positions,
	 sizeof(query_concatenated_positions[0]) * size,
	 cudaMemcpyDefault, stream));

	 CUDA_CHECK_RETURN(
	 cudaMemcpyAsync(d_database_positions, database_positions,
	 sizeof(database_positions[0]) * size, cudaMemcpyDefault,
	 stream));
	 */
	size_t loader_shared_memory_size =
			sizeof(uint32_t) * gapped_extender_gpu_kernel::kNumberThreads
					+ GroupLoader<
							const packed_alphabet_code::PackedAlphabetCode *,
							gapped_extender_gpu_kernel::kLoadLength>::GetTotalSharedMemorySize(
							gapped_extender_gpu_kernel::kNumberThreads);
	size_t dp_row_shared_memory_size =
			gapped_extender_gpu_kernel::kNumberThreads
					* (gapped_extender_gpu_kernel::kDpRowLength * sizeof(int)
							* 2);
	size_t shared_memory_size = max(loader_shared_memory_size,
			dp_row_shared_memory_size);
	gapped_extender_gpu_kernel::ExtendOneSideScoreOnlyKernel<<<
	gapped_extender_gpu_kernel::kNumberBlocks,
	gapped_extender_gpu_kernel::kNumberThreads, shared_memory_size,
	stream>>>(d_concatenated_query_sequence_, d_database_sequence_,
			size, reverse, sequence_delimiter_, d_score_matrix_,
			number_letters_ + 1, gap_open_, gap_extention_, cutoff_,
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
