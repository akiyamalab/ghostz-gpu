/*
 * packed_alphabet_code_block.h
 *
 *  Created on: Aug 21, 2014
 *      Author: shu
 */

#pragma once

#include <stdint.h>
#include "alphabet_coder.h"
#include "cuda_common.h"

namespace packed_alphabet_code {
typedef uint64_t PackedAlphabetCode;
static const uint32_t kPackedAlphabetCodeSize = sizeof(PackedAlphabetCode);
static const uint32_t kCodeBitSize = 5;
static const uint32_t kNumberOfCodesInBlock = (8 * kPackedAlphabetCodeSize)
		/ kCodeBitSize;
static const uint32_t kMask = (1 << kCodeBitSize) - 1;

size_t GetPackedAlphabetCodeSequenceLength(size_t code_sequence_length);
int PackAlphabetCodeSequence(const AlphabetCoder::Code *code_sequence,
		const size_t code_sequence_length,
		PackedAlphabetCode* packed_code_sequence);

__forceinline__ __host__           __device__ uint32_t GetBlockPosition(
		const uint32_t position) {
	return position / kNumberOfCodesInBlock;
}

template<cuda_common::Direction TDirection>
__forceinline__ __host__ __device__ int InitPackedAlphabetCode(
		const PackedAlphabetCode* code_blocks, uint32_t original_start,
		int *next_block_position_ptr, uint32_t *block_sift_ptr,
		uint32_t *next_block_sift_ptr,
		PackedAlphabetCode *temp_code_block_ptr) {
	return 1;
}

template<>
__forceinline__ __host__ __device__ int InitPackedAlphabetCode<
		cuda_common::kFoward>(const PackedAlphabetCode* code_blocks,
		uint32_t original_start, int *next_block_position_ptr,
		uint32_t *block_sift_ptr, uint32_t *next_block_sift_ptr,
		PackedAlphabetCode *temp_code_block_ptr) {
	const int increment = 1;
	const uint32_t block_start = GetBlockPosition(original_start);
	const uint32_t code_offset = original_start
			- block_start * kNumberOfCodesInBlock;

	const uint32_t shift_size = (kCodeBitSize * code_offset);

	const int block_shift = shift_size;
	*block_sift_ptr = block_shift;
	*next_block_sift_ptr = kCodeBitSize * kNumberOfCodesInBlock - block_shift;

	*temp_code_block_ptr = code_blocks[block_start];
	*next_block_position_ptr = block_start + increment;
	return 0;
}

template<>
__forceinline__ __host__ __device__ int InitPackedAlphabetCode<
		cuda_common::kReverse>(const PackedAlphabetCode* code_blocks,
		uint32_t original_start, int *next_block_position_ptr,
		uint32_t *block_sift_ptr, uint32_t *next_block_sift_ptr,
		PackedAlphabetCode *temp_code_block_ptr) {
	const int increment = -1;
	const uint32_t block_start = GetBlockPosition(original_start);
	const uint32_t code_offset = original_start
			- block_start * kNumberOfCodesInBlock;

	const uint32_t remain_bit_size = kPackedAlphabetCodeSize * 8
			- kNumberOfCodesInBlock * kCodeBitSize;
	const uint32_t block_shift = remain_bit_size
			+ (kCodeBitSize * (kNumberOfCodesInBlock - code_offset - 1));
	*block_sift_ptr = block_shift;
	*next_block_sift_ptr = kCodeBitSize * kNumberOfCodesInBlock - block_shift;

	*temp_code_block_ptr = code_blocks[block_start];
	*next_block_position_ptr = (int) block_start + increment;
	return 0;
}

template<cuda_common::Direction TDirection>
__forceinline__ __host__                              __device__ PackedAlphabetCode GetPackedAlphabetCode(
		const PackedAlphabetCode* code_bloks, const uint32_t block_shift,
		const uint32_t next_block_sift, int *block_position_ptr,
		PackedAlphabetCode *temp_code_block_ptr) {
	return 1;
}

template<>
__forceinline__ __host__                              __device__ PackedAlphabetCode GetPackedAlphabetCode<
		cuda_common::kFoward>(const PackedAlphabetCode* code_bloks,
		const uint32_t block_shift, const uint32_t next_block_sift,
		int *block_position_ptr, PackedAlphabetCode *temp_code_block_ptr) {
	const int increment = 1;
	PackedAlphabetCode ret = *temp_code_block_ptr;
	ret >>= block_shift;
	*temp_code_block_ptr = code_bloks[*block_position_ptr];
	*block_position_ptr += increment;
	ret |= *temp_code_block_ptr << next_block_sift;
	return ret;
}

template<>
__forceinline__ __host__                              __device__ PackedAlphabetCode GetPackedAlphabetCode<
		cuda_common::kReverse>(const PackedAlphabetCode* code_bloks,
		const uint32_t block_shift, const uint32_t next_block_sift,
		int *block_position_ptr, PackedAlphabetCode *temp_code_block_ptr) {
	const int increment = -1;
	PackedAlphabetCode ret = *temp_code_block_ptr;
	ret <<= block_shift;
	*temp_code_block_ptr = code_bloks[*block_position_ptr];
	*block_position_ptr += increment;
	ret |= *temp_code_block_ptr >> next_block_sift;
	return ret;
}

template<cuda_common::Direction TDirection>
__forceinline__ __host__              __device__ uint32_t GetAlphabetCode(
		const PackedAlphabetCode code_bloks, int i) {
	return 255;
}

template<>
__forceinline__ __host__                              __device__ uint32_t GetAlphabetCode<
		cuda_common::kFoward>(const PackedAlphabetCode code_bloks, int i) {
	return (uint32_t) (code_bloks >> (kCodeBitSize * i)) & kMask;
}

template<>
__forceinline__ __host__                              __device__ uint32_t GetAlphabetCode<
		cuda_common::kReverse>(const PackedAlphabetCode code_bloks, int i) {
	const uint32_t remain_bit_size = kPackedAlphabetCodeSize * 8
			- kNumberOfCodesInBlock * kCodeBitSize;
	const uint32_t sift = kCodeBitSize * (kNumberOfCodesInBlock - 1 - i)
			+ remain_bit_size;
	return (uint32_t) (code_bloks >> sift) & kMask;
}

}

