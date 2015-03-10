#include "distance_calculator_gpu.h"
#include "cuda_common.h"
#include "packed_alphabet_code.h"

/*
namespace distance_calculator_gpu_ref_kernel {

__device__ DistanceCalculatorGpu::Distance CalculateDisntaceDeviceCach(
		packed_alphabet_code::PackedAlphabetCode sequence0_multi_code_cache,
		packed_alphabet_code::PackedAlphabetCode sequence1_multi_code_cache,
		const AlphabetCoder::Code *redueced_code_map,
		//const int* __restrict__ redueced_code_map,
		const uint32_t length) {
	DistanceCalculatorGpu::Distance distance = 0;
	const uint32_t mask = (1 << 8) - 1;

	for (uint32_t i = 0; i < length; ++i) {
		const uint32_t shift = (8 * i);
		const AlphabetCoder::Code s0_c = (sequence0_multi_code_cache >> shift)
				& mask;
		const AlphabetCoder::Code s1_c = (sequence1_multi_code_cache >> shift)
				& mask;
#if 0
		if ((blockDim.x * blockIdx.x + threadIdx.x) == 2) {
			printf("%d, %d\n", s0_c, s1_c);
		}
#endif
		distance +=
				(redueced_code_map[s0_c] != redueced_code_map[s1_c]) ? 1 : 0;
	}

	return distance;
}

__device__ int CalculateDisntaceDevice(
		const packed_alphabet_code::PackedAlphabetCode* sequence0_multi_code,
		const uint32_t s0_start_p,
		const packed_alphabet_code::PackedAlphabetCode* sequence1_multi_code,
		const uint32_t s1_start_p, const AlphabetCoder::Code *redueced_code_map,
		//const int* __restrict__ redueced_code_map,
		const uint32_t subsequence_length) {

	const uint32_t s0_multi_code_start_p = s0_start_p
			/ packed_alphabet_code::kNumberCodesInBlock;
	const uint32_t s1_multi_code_start_p = s1_start_p
			/ packed_alphabet_code::kNumberCodesInBlock;
	const uint32_t s0_multi_code_p_offset = s0_start_p
			- s0_multi_code_start_p * packed_alphabet_code::kNumberCodesInBlock;
	const uint32_t s1_multi_code_p_offset = s1_start_p
			- s1_multi_code_start_p * packed_alphabet_code::kNumberCodesInBlock;

	uint32_t s0_multi_code_p = s0_multi_code_start_p;
	uint32_t s1_multi_code_p = s1_multi_code_start_p;
	packed_alphabet_code::PackedAlphabetCode s0_cache[2];
	packed_alphabet_code::PackedAlphabetCode s1_cache[2];

	const uint32_t s0_shift_size = (8 * s0_multi_code_p_offset);
	const uint32_t s1_shift_size = (8 * s1_multi_code_p_offset);
	const uint32_t s0_sift0 = s0_shift_size;
	const uint32_t s0_sift1 = 8 * packed_alphabet_code::kNumberCodesInBlock
			- s0_sift0;
	const uint32_t s1_sift0 = s1_shift_size;
	const uint32_t s1_sift1 = 8 * packed_alphabet_code::kNumberCodesInBlock
			- s1_sift0;
#if 0
	if ((blockDim.x * blockIdx.x + threadIdx.x) == 2) {
		printf("s0 shift %d, %d\n", s0_sift0, s0_sift1);
	}

	if ((blockDim.x * blockIdx.x + threadIdx.x) == 2) {
		printf("s1 shift %d, %d\n", s1_sift0, s1_sift1);
	}
#endif
	s0_cache[0] = sequence0_multi_code[s0_multi_code_p];
	++s0_multi_code_p;
	s1_cache[0] = sequence1_multi_code[s1_multi_code_p];
	++s1_multi_code_p;
#if 0
	if ((blockDim.x * blockIdx.x + threadIdx.x) == 2) {
		printf("cached0 %d, %d\n", s0_cache[0], s1_cache[0]);
	}
#endif
	s0_cache[0] >>= s0_sift0;
	s1_cache[0] >>= s1_sift0;
#if 0
	if ((blockDim.x * blockIdx.x + threadIdx.x) == 2) {
		printf("edited cached0 %d, %d\n", s0_cache[0], s1_cache[0]);
	}
#endif
	DistanceCalculatorGpu::Distance distance = 0;
	for (uint32_t offset = 0; offset < subsequence_length; offset +=
			packed_alphabet_code::kNumberCodesInBlock) {
#if 0
		if (2321 <= s0_multi_code_p || 11959 <= s1_multi_code_p) {
			printf("job id : %d \n", (blockDim.x * blockIdx.x + threadIdx.x));
		}
#endif
		s0_cache[1] = sequence0_multi_code[s0_multi_code_p];
		++s0_multi_code_p;
		s1_cache[1] = sequence1_multi_code[s1_multi_code_p];
		++s1_multi_code_p;
#if 0
		if ((blockDim.x * blockIdx.x + threadIdx.x) == 2) {
			printf("cached1 %d, %d\n", s0_cache[1], s1_cache[1]);
		}
#endif
		packed_alphabet_code::PackedAlphabetCode tmp = 0;
		tmp = s0_cache[1] << s0_sift1;
		s0_cache[0] |= tmp;
		s0_cache[1] >>= s0_sift0;
		tmp = s1_cache[1] << s1_sift1;
		s1_cache[0] |= tmp;
		s1_cache[1] >>= s1_sift0;
#if 0
		if ((blockDim.x * blockIdx.x + threadIdx.x) == 2) {
			printf("input cached0 %d, %d\n", s0_cache[0], s1_cache[0]);
			printf("input cached1 %d, %d\n", s0_cache[1], s1_cache[1]);
		}
#endif
		uint32_t length =
				(offset + packed_alphabet_code::kNumberCodesInBlock)
						<= subsequence_length ?
						packed_alphabet_code::kNumberCodesInBlock :
						subsequence_length - offset;
		distance += CalculateDisntaceDeviceCach(s0_cache[0], s1_cache[0],
				redueced_code_map, length);

		s0_cache[0] = s0_cache[1];
		s1_cache[0] = s1_cache[1];
	}

	return distance;
}

__global__ void CalculateDisntacesKernel(
		const packed_alphabet_code::PackedAlphabetCode* sequence0_multi_code,
		const packed_alphabet_code::PackedAlphabetCode* sequence1_multi_code,
		const uint32_t number_distances, const uint32_t * sequence0_positions,
		const uint32_t* sequence1_positions,
		const AlphabetCoder::Code *redueced_code_map,
		//const int* __restrict__ redueced_code_map,
		const uint32_t subsequence_length,
		DistanceCalculatorGpu::Distance *distances) {
	const uint32_t foward_direction = subsequence_length / 2;
	const uint32_t thread_id = blockDim.x * blockIdx.x + threadIdx.x;
	const uint32_t thread_id_skip = gridDim.x * blockDim.x;

	for (uint32_t i = thread_id; i < number_distances; i += thread_id_skip) {
		const uint32_t s0_p = sequence0_positions[i]
				+ cuda_common::kOneSideMarginSize - foward_direction;
		const uint32_t s1_p = sequence1_positions[i]
				+ cuda_common::kOneSideMarginSize - foward_direction;
		distances[i] = CalculateDisntaceDevice(sequence0_multi_code, s0_p,
				sequence1_multi_code, s1_p, redueced_code_map,
				subsequence_length);
	}
	return;
}

}
*/
