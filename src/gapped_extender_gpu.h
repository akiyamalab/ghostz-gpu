/*
 * gapped_extender_gpu.h
 *
 *  Created on: Aug 11, 2014
 *      Author: shu
 */

#ifndef GAPPED_EXTENDER_GPU_H_
#define GAPPED_EXTENDER_GPU_H_

#include <stdint.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/system/cuda/experimental/pinned_allocator.h>
#include <thrust/copy.h>
#include "alphabet_coder.h"
#include "packed_alphabet_code.h"

class GappedExtenderGpu {
public:
	struct DpColumnCell {
		int h;
		int e;
	};
	static const size_t kMaxSequence0Length = 1 << 6;
	static const int kInitScore = -(1 << 10);

	GappedExtenderGpu();
	virtual ~GappedExtenderGpu();

	int SetQueries(AlphabetCoder::Code sequence_delimiter,
			const packed_alphabet_code::PackedAlphabetCode *d_concatenated_query_sequence);
	int SetDatabase(
			const packed_alphabet_code::PackedAlphabetCode *d_database_sequence);
	int SetScoreParameters(const int *d_score_matrix, uint32_t number_letters,
			int gap_open, int gap_extention, int cutoff);

	int ConvertToGappedExtensionSeedsAsync(size_t size, bool reverse,
			uint32_t* seed_positions, uint32_t* query_concatenated_positions,
			uint32_t* database_positions, uint32_t* d_seed_positions,
			uint32_t* d_temp_array, uint32_t* d_query_concatenated_positions,
			uint32_t* d_database_positions, cudaStream_t &stream) const;
	int ExtendOneSideScoreOnlyAsync(size_t size, bool reverse,
			uint32_t* query_concatenated_positions,
			uint32_t* database_positions, int* scores,
			uint32_t* d_query_concatenated_positions,
			uint32_t* d_database_positions, int* d_scores,
			cudaStream_t &stream) const;

private:
	uint32_t number_letters_;
	int gap_open_;
	int gap_extention_;
	int cutoff_;
	AlphabetCoder::Code sequence_delimiter_;
	const packed_alphabet_code::PackedAlphabetCode *d_database_sequence_;
	const packed_alphabet_code::PackedAlphabetCode *d_concatenated_query_sequence_;
	const int *d_score_matrix_;
	const int *d_ungapped_extension_cutoffs_;
};

#endif /* GAPPED_EXTENDER_GPU_H_ */
