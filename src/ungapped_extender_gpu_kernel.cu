/*
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/system/cuda/experimental/pinned_allocator.h>
#include <thrust/copy.h>
#include <vector>
#include <assert.h>

#include "ungapped_extender_gpu.h"
#include "group_loader.h"
#include "packed_alphabet_code.h"
#include "score_matrix.h"
#include "cuda_common.h"

const int debug_q_p = 119;
const int debug_db_p = 8963 + cuda_common::kOneSideMarginSize;

using namespace std;

__device__ int UngappedExtendOneSideCacheRef(
		packed_alphabet_code::PackedAlphabetCode sequence0_multi_code_cache,
		packed_alphabet_code::PackedAlphabetCode sequence1_multi_code_cache,
		const AlphabetCoder::Code sequence_delimiter, bool reversed,
		const int*  __restrict__  score_matrix,
		const uint32_t number_letters, const int cutoff, const int trigger,
		int *score_ptr, int *best_score_ptr) {
	int score = *score_ptr;
	int best_score = *best_score_ptr;
	int threshold = best_score - cutoff;
	bool stop_flag = false;
	const uint32_t mask = (1 << 8) - 1;
	uint32_t sift_offset =
			reversed ? 8 * (packed_alphabet_code::kNumberCodesInBlock - 1) : 0;

	for (uint32_t i = 0; i < packed_alphabet_code::kNumberCodesInBlock; ++i) {
		const uint32_t shift_size = (8 * i);
		const uint32_t shift = reversed ? sift_offset - shift_size : shift_size;
		const AlphabetCoder::Code s0_c = (sequence0_multi_code_cache >> shift)
				& mask;
		const AlphabetCoder::Code s1_c = (sequence1_multi_code_cache >> shift)
				& mask;

		stop_flag = stop_flag || s0_c == sequence_delimiter
				|| s1_c == sequence_delimiter;
		score += stop_flag ? 0 : score_matrix[s0_c * number_letters + s1_c];
		best_score = score > best_score ? score : best_score;
		threshold = best_score - cutoff;
		stop_flag = stop_flag || score <= threshold || best_score > trigger;
		if (stop_flag) {
			break;
		}
	}
	*score_ptr = score;
	*best_score_ptr = best_score;
	return stop_flag;
	//return true;
}

__device__ int UngappedExtendOneSideDeviceRef(
		const packed_alphabet_code::PackedAlphabetCode* sequence0_multi_code,
		const uint32_t s0_start_p,
		const packed_alphabet_code::PackedAlphabetCode* sequence1_multi_code,
		const uint32_t s1_start_p, const AlphabetCoder::Code sequence_delimiter,
		bool reversed, const int* __restrict__ score_matrix,
		const uint32_t number_letters, const int cutoff, const int trigger,
		int current_best_score) {
	int score = current_best_score;
	int best_score = current_best_score;
	int increment = reversed ? -1 : 1;

	const uint32_t s0_multi_code_start_p = s0_start_p
			/ packed_alphabet_code::kNumberCodesInBlock;
	const uint32_t s1_multi_code_start_p = s1_start_p
			/ packed_alphabet_code::kNumberCodesInBlock;
	const uint32_t s0_multi_code_p_offset = s0_start_p
			- s0_multi_code_start_p * packed_alphabet_code::kNumberCodesInBlock;
	const uint32_t s1_multi_code_p_offset = s1_start_p
			- s1_multi_code_start_p * packed_alphabet_code::kNumberCodesInBlock;

#if 0
	if ((debug_db_p - 1) == s1_start_p || debug_db_p == s1_start_p) {
		printf("block p  %d, %d\n", s0_multi_code_start_p,
				s1_multi_code_start_p);
	}
#endif

	uint32_t s0_multi_code_p = s0_multi_code_start_p;
	uint32_t s1_multi_code_p = s1_multi_code_start_p;
	packed_alphabet_code::PackedAlphabetCode s0_cache[2];
	packed_alphabet_code::PackedAlphabetCode s1_cache[2];

	const uint32_t s0_shift_size = (8 * s0_multi_code_p_offset);
	const uint32_t s1_shift_size = (8 * s1_multi_code_p_offset);
	const uint32_t sift_offset =
			reversed ? 8 * (packed_alphabet_code::kNumberCodesInBlock - 1) : 0;
	const uint32_t s0_sift0 =
			reversed ? sift_offset - s0_shift_size : s0_shift_size;
	const uint32_t s0_sift1 = 8 * packed_alphabet_code::kNumberCodesInBlock
			- s0_sift0;
	const uint32_t s1_sift0 =
			reversed ? sift_offset - s1_shift_size : s1_shift_size;
	const uint32_t s1_sift1 = 8 * packed_alphabet_code::kNumberCodesInBlock
			- s1_sift0;

	s0_cache[0] = sequence0_multi_code[s0_multi_code_p];
	s0_multi_code_p += increment;
	s1_cache[0] = sequence1_multi_code[s1_multi_code_p];
	s1_multi_code_p += increment;

#if 0
	if ((debug_db_p - 1) == s1_start_p || debug_db_p == s1_start_p) {
		printf("cached0 %ld, %ld\n", s0_cache[0], s1_cache[0]);
	}
#endif

	if (!reversed) {
		s0_cache[0] >>= s0_sift0;
		s1_cache[0] >>= s1_sift0;
	} else {
		s0_cache[0] <<= s0_sift0;
		s1_cache[0] <<= s1_sift0;
	}
#if 0
	if ((blockDim.x * blockIdx.x + threadIdx.x) == 3) {
		printf("edited cached0 %d, %d\n", s0_cache[0], s1_cache[0]);
	}
#endif
	while (1) {
#if 0
		if (2321 <= s0_multi_code_p || 11959 <= s1_multi_code_p) {
			printf("job id : %d \n", (blockDim.x * blockIdx.x + threadIdx.x));
		}
#endif
		s0_cache[1] = sequence0_multi_code[s0_multi_code_p];
		s0_multi_code_p += increment;
		s1_cache[1] = sequence1_multi_code[s1_multi_code_p];
		s1_multi_code_p += increment;
#if 0
		if ((blockDim.x * blockIdx.x + threadIdx.x) == 3) {
			printf("cached1 %d, %d\n", s0_cache[1], s1_cache[1]);
		}
#endif
		packed_alphabet_code::PackedAlphabetCode tmp = 0;
		if (!reversed) {
			tmp = s0_cache[1] << s0_sift1;
			s0_cache[0] |= tmp;
			s0_cache[1] >>= s0_sift0;
			tmp = s1_cache[1] << s1_sift1;
			s1_cache[0] |= tmp;
			s1_cache[1] >>= s1_sift0;
		} else {
			tmp = s0_cache[1] >> s0_sift1;
			s0_cache[0] |= tmp;
			s0_cache[1] <<= s0_sift0;
			tmp = s1_cache[1] >> s1_sift1;
			s1_cache[0] |= tmp;
			s1_cache[1] <<= s1_sift0;
		}
#if 0
		if ((debug_db_p - 1) == s1_start_p || debug_db_p == s1_start_p) {
			if (reversed) {
				printf("reverse ");
			} else {
				printf("foward ");
			}
			printf("multi_code %ld, %ld\n", s0_cache[0], s1_cache[0]);
		}
#endif
		if (UngappedExtendOneSideCacheRef(s0_cache[0], s1_cache[0],
				sequence_delimiter, reversed, score_matrix, number_letters,
				cutoff, trigger, &score, &best_score)) {
			break;
		}
		s0_cache[0] = s0_cache[1];
		s1_cache[0] = s1_cache[1];
	}

	return best_score;
}

__global__ void UngappedExtendKernelRef(
		const packed_alphabet_code::PackedAlphabetCode* sequence0_multi_code,
		const packed_alphabet_code::PackedAlphabetCode* sequence1_multi_code,
		const uint32_t number_extensions, const uint32_t* query_ids,
		const uint32_t * sequence0_positions,
		const uint32_t* sequence1_positions,
		const AlphabetCoder::Code sequence_delimiter,
		const int* __restrict__ score_matrix,
		const uint32_t number_letters, const int* cutoffs, const int* triggers,
		int* best_scores) {

	const uint32_t thread_id = blockDim.x * blockIdx.x + threadIdx.x;
	const uint32_t thread_id_skip = gridDim.x * blockDim.x;

	for (uint32_t i = thread_id; i < number_extensions; i += thread_id_skip) {
		const uint32_t s0_p = sequence0_positions[i]
				+ cuda_common::kOneSideMarginSize;
		const uint32_t s1_p = sequence1_positions[i]
				+ cuda_common::kOneSideMarginSize;
		const uint32_t query_id = query_ids[i];
		const int cutoff = cutoffs[query_id];
		const int trigger = triggers[query_id];
		int current_best_score = 0;
		current_best_score = UngappedExtendOneSideDeviceRef(
				sequence0_multi_code, s0_p - 1, sequence1_multi_code, s1_p - 1,
				sequence_delimiter, true, score_matrix, number_letters, cutoff,
				trigger, current_best_score);

		if (current_best_score <= trigger) {
			current_best_score = UngappedExtendOneSideDeviceRef(
					sequence0_multi_code, s0_p, sequence1_multi_code, s1_p,
					sequence_delimiter, false, score_matrix, number_letters,
					cutoff, trigger, current_best_score);
		}
		best_scores[i] = current_best_score;
	}

	return;
}
*/
