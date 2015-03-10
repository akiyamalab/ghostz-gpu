/*
 * gapped_extender_gpu_ref.cu
 *
 *  Created on: 2014/08/23
 *      Author: shu
 */

#ifndef GAPPED_EXTENDER_GPU_REF_CU_
#define GAPPED_EXTENDER_GPU_REF_CU_

/*
#include "gapped_extender_gpu.h"

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/system/cuda/experimental/pinned_allocator.h>
#include <thrust/copy.h>
#include <vector>
#include <assert.h>

#include "score_matrix.h"
#include "cuda_common.h"

using namespace std;

__device__ int InitCachedSequence(const AlphabetCoder::Code* sequence,
		const AlphabetCoder::Code sequence_delimiter, int increment,
		AlphabetCoder::Code* cached_sequence) {
	AlphabetCoder::Code c = 0;
	int p = 0;
	for (c = 0, p = 0; c != sequence_delimiter; p += increment) {
		c = sequence[p];
		cached_sequence[p] = c;
	}
	return 0;
}

__device__ int InitScoreArray(const AlphabetCoder::Code* sequence,
		const AlphabetCoder::Code sequence_delimiter, int gap_init,
		int gap_extention, int default_cutoff, int increment,
		GappedExtenderGpu::DpCell* score_array) {
	int score = -gap_init;
	score_array[0].best = 0;
	score_array[0].best_gap = -gap_init;
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
		score_array[array_index].best = score;
		score_array[array_index].best_gap = score - gap_extention;
		score -= gap_extention;
	}

	return array_index;
}

__device__ int UpdateDpCells(const AlphabetCoder::Code* sequence0,
		const AlphabetCoder::Code s1_c, const int sequence1_position,
		const AlphabetCoder::Code sequence_delimiter, int increment, int cutoff,
		const int *score_matrix_row, const int gap_init,
		const int gap_extention, GappedExtenderGpu::DpCell* score_array,
		int *array_start_ptr, int *array_end_ptr,
		int* max_score_sequence0_position, int* max_score_sequence1_position,
		int* max_score_ptr) {

	int score = GappedExtenderGpu::kInitScore;
	int score_gap_row = GappedExtenderGpu::kInitScore;
	int max_score = *max_score_ptr;
	int array_start = *array_start_ptr;
	int array_end = *array_end_ptr;
	int sequence0_position = array_start * increment;
	int prev_score = score_array[array_start].best;
	score_array[array_start].best = score_array[array_start].best_gap;
	score_array[array_start].best_gap -= gap_extention;
	int array_last_index = array_start;

	for (int array_index = array_start + 1; array_index < array_end;
			++array_index, sequence0_position += increment) {
		score = prev_score + score_matrix_row[sequence0[sequence0_position]];
		prev_score = score_array[array_index].best;
		int score_gap_column = score_array[array_index].best_gap;

		score = score < score_gap_column ? score_gap_column : score;
		score = score < score_gap_row ? score_gap_row : score;

		if (max_score - score > cutoff) {
			array_start += array_start + 1 == array_index ? 1 : 0;
			score_array[array_index].best = GappedExtenderGpu::kInitScore;
			score_array[array_index].best_gap = GappedExtenderGpu::kInitScore;
			score_gap_row = GappedExtenderGpu::kInitScore;
		} else {
			array_last_index = array_index;
			if (score > max_score) {
				max_score = score;
				*max_score_sequence0_position = sequence0_position;
				*max_score_sequence1_position = sequence1_position;
			}

			score_array[array_index].best_gap = max(score - gap_init,
					score_gap_column - gap_extention);
			score_gap_row = max(score - gap_init,
					score_gap_row - gap_extention);
			score_array[array_index].best = score;
		}
	}

	if (array_start + 1 != array_end) {
		if (array_last_index < array_end - 1) {
			array_end = array_last_index + 1;
		} else {
			while (score_gap_row >= (max_score - cutoff)
					&& sequence0[sequence0_position] != sequence_delimiter) {
				score_array[array_end].best = score_gap_row;
				score_array[array_end].best_gap = score_gap_row - gap_init;
				score_gap_row -= gap_extention;
				++array_end;
				sequence0_position += increment;
			}
			if (sequence0[sequence0_position] != sequence_delimiter) {
				score_array[array_end].best = GappedExtenderGpu::kInitScore;
				score_array[array_end].best_gap = GappedExtenderGpu::kInitScore;
				++array_end;
			}
		}
	}

	*array_start_ptr = array_start;
	*array_end_ptr = array_end;
	*max_score_ptr = max_score;
	return 0;
}

__device__ int ExtendOneSideScoreOnlyDevice(
		const AlphabetCoder::Code* concatenated_sequence0,
		const uint32_t sequence_0_offset,
		const AlphabetCoder::Code* concatenated_sequence1,
		const uint32_t sequence_1_offset, const bool reverse,
		const AlphabetCoder::Code sequence_delimiter, const int* score_matrix,
		const uint32_t number_letters, const int gap_open,
		const int gap_extention, const int cutoff,
		GappedExtenderGpu::DpCell* dp_cells, uint32_t* best_sequence0_position,
		uint32_t* best_sequence1_position, int* best_score) {
	int increment = reverse ? -1 : 1;
	int max_score = 0;
	int max_score_sequence0_position = -increment;
	int max_score_sequence1_position = -increment;

	const AlphabetCoder::Code* sequence0 = concatenated_sequence0
			+ sequence_0_offset;
	const AlphabetCoder::Code* sequence1 = concatenated_sequence1
			+ sequence_1_offset;

	AlphabetCoder::Code cached_sequence0_mem[GappedExtenderGpu::kMaxSequence0Length];
	AlphabetCoder::Code *cached_sequence0 =
			reverse ?
					&cached_sequence0_mem[GappedExtenderGpu::kMaxSequence0Length
							- 1] :
					&cached_sequence0_mem[0];
	InitCachedSequence(sequence0, sequence_delimiter, increment,
			cached_sequence0);
	int sequence1_position = 0;
	GappedExtenderGpu::DpCell *score_array = dp_cells;

	int array_start = 0;
	int array_end = 0;
	int gap_init = gap_open + gap_extention;
	array_end = InitScoreArray(cached_sequence0, sequence_delimiter, gap_init,
			gap_extention, cutoff, increment, score_array);

#if 0
	printf("\n");
	printf("      ");
	for (int x = 0; sequence0[x] != sequence_delimiter; x += increment) {
		printf("%3d", sequence0[x]);
	}
	printf("\n");
#endif

	bool stop_flag = false;
	while (!stop_flag) {
#if 0
		printf("%3d", sequence1[sequence1_position - increment]);
		for (int x = 0; x < array_start; ++x) {
			printf("   ");
		}
		for (int x = array_start; x < array_end; ++x) {
			printf("%3d", score_array[x].best);
			//fprintf(stderr, "%3d", insertion_sequence1_row[x]);
		}
		printf("\n");
#endif
		AlphabetCoder::Code s1_c = sequence1[sequence1_position];
		if (s1_c == sequence_delimiter) {
			stop_flag = true;
		} else {
			UpdateDpCells(cached_sequence0, s1_c, sequence1_position,
					sequence_delimiter, increment, cutoff,
					score_matrix + s1_c * number_letters, gap_init,
					gap_extention, score_array, &array_start, &array_end,
					&max_score_sequence0_position,
					&max_score_sequence1_position, &max_score);
			sequence1_position += increment;
		}
#if 0
		// debug //////////////////
		if (sequence1_position >= 10) {
			stop_flag = true;
			break;
		}
		//////////////////////////
#endif
		if (array_start + 1 == array_end) {
			stop_flag = true;
			break;
		}

	}

	*best_score = max_score;
	*best_sequence0_position = max_score_sequence0_position;
	*best_sequence1_position = max_score_sequence1_position;
	return 0;
}

__global__ void
__launch_bounds__(128, 1) ExtendOneSideScoreOnlyKernel(
		const AlphabetCoder::Code* concatenated_sequence0,
		const AlphabetCoder::Code* concatenated_sequence1,
		const uint32_t number_extensions, const bool reverse,
		const AlphabetCoder::Code sequence_delimiter, const int* score_matrix,
		const uint32_t number_letters, const int gap_open,
		const int gap_extention, const int cutoff,
		uint32_t* sequence0_positions, uint32_t* sequence1_positions,
		int* best_scores) {

	const uint32_t thread_id = blockDim.x * blockIdx.x + threadIdx.x;
	const uint32_t thread_id_skip = gridDim.x * blockDim.x;
	GappedExtenderGpu::DpCell dp_cells[GappedExtenderGpu::kMaxSequence0Length];

	for (uint32_t i = thread_id; i < number_extensions; i += thread_id_skip) {
		uint32_t sequence_0_offset = sequence0_positions[i]
				+ cuda_common::kOneSideMarginSize;
		uint32_t sequence_1_offset = sequence1_positions[i]
				+ cuda_common::kOneSideMarginSize;
		uint32_t best_sequence_0_p = 0;
		uint32_t best_sequence_1_p = 0;
		int best_score = 0;
		ExtendOneSideScoreOnlyDevice(concatenated_sequence0, sequence_0_offset,
				concatenated_sequence1, sequence_1_offset, reverse,
				sequence_delimiter, score_matrix, number_letters, gap_open,
				gap_extention, cutoff, &dp_cells[0], &best_sequence_0_p,
				&best_sequence_1_p, &best_score);
		sequence0_positions[i] = sequence_0_offset + best_sequence_0_p
				- cuda_common::kOneSideMarginSize;
		sequence1_positions[i] = sequence_1_offset + best_sequence_1_p
				- cuda_common::kOneSideMarginSize;
		best_scores[i] = best_score;
	}
	return;
}

// call
/*
 ExtendOneSideScoreOnlyKernel<<<1024, 128, 0, stream>>>(d_concatenated_query_sequence_,
			d_database_sequence_,
			size, reverse, sequence_delimiter_, d_score_matrix_, number_letters_ + 1, gap_open_, gap_extention_, cutoff_,
			d_query_concatenated_positions,
			d_database_positions,
			d_scores);
 *
 */

#endif /* GAPPED_EXTENDER_GPU_REF_CU_ */

