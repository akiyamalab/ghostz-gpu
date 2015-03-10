/*
 * distance_calculator_gpu.h
 *
 *  Created on: Jul 28, 2014
 *      Author: shu
 */

#ifndef DISTANCE_CALCULATOR_GPU_H_
#define DISTANCE_CALCULATOR_GPU_H_

#include "distance_calculator.h"
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/system/cuda/experimental/pinned_allocator.h>
#include <thrust/copy.h>
#include "packed_alphabet_code.h"

class DistanceCalculatorGpu {
public:
	typedef DistanceCalculator::Distance Distance;

	DistanceCalculatorGpu();
	virtual ~DistanceCalculatorGpu();

	int SetQueries(const packed_alphabet_code::PackedAlphabetCode *query_sequence);
	int SetDatabase(const packed_alphabet_code::PackedAlphabetCode *database_sequence);
	int SetReducedCodeMap(const AlphabetCoder::Code *redueced_code_map);
	int SetSubsequenceLength(uint32_t length);
	int CalculateDistancesAsync(size_t size, uint32_t* query_concatenated_centers,
			uint32_t* database_centers,
			Distance* distances,
			uint32_t* d_query_concatenated_centers,
			uint32_t* d_database_centers,
			Distance* d_distances, cudaStream_t &stream);
private:
	uint32_t subsequence_length_;
	const packed_alphabet_code::PackedAlphabetCode *d_concatenated_query_sequence_;
	const packed_alphabet_code::PackedAlphabetCode *d_database_sequence_;
	const AlphabetCoder::Code *d_redueced_code_map_;
};

#endif /* DISTANCE_CALCULATOR_GPU_H_ */
