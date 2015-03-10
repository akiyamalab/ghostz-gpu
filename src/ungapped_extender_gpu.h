/*
 * ungapped_extender_gpu.h
 *
 *  Created on: Jul 21, 2014
 *      Author: shu
 */

#ifndef UNGAPPED_EXTENDER_GPU_H_
#define UNGAPPED_EXTENDER_GPU_H_

#include <stdint.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/system/cuda/experimental/pinned_allocator.h>
#include <thrust/copy.h>
#include "alphabet_coder.h"
#include "packed_alphabet_code.h"

class ScoreMatrix;

class UngappedExtenderGpu {
public:
	UngappedExtenderGpu();
	virtual ~UngappedExtenderGpu();
	int SetQuerySeedDataList(const uint32_t* d_query_ids,
			const uint32_t* d_query_positions);
	int SetQueries(AlphabetCoder::Code sequence_delimiter,
			const packed_alphabet_code::PackedAlphabetCode *d_concatenated_query_sequence,
			const int *d_ungapped_extension_cutoffs,
			const int *d_gapped_extension_triggers);
	int SetDatabase(
			const packed_alphabet_code::PackedAlphabetCode *d_database_sequence);
	int SetScoreMatrix(const int *d_score_matrix, uint32_t number_letters);

	int SetQuerySeeds(size_t query_seeds_size, size_t size,
			uint32_t* query_seed_ids, uint32_t* query_seed_starts,
			size_t temp_storage_bytes, void* d_temp_storage,
			uint32_t* d_query_seed_ids, uint32_t* d_query_seed_starts,
			uint32_t* d_query_ids, uint32_t* d_query_concatenated_positions,
			uint32_t* temp_array, cudaStream_t &stream) const;

	int ExtendWithTriggerAsync(size_t size, uint32_t* query_ids,
			uint32_t* query_concatenated_positions,
			uint32_t* database_positions, char* flags, uint32_t* d_query_ids,
			uint32_t* d_query_concatenated_positions,
			uint32_t* d_database_positions, char* d_flags, int* d_temp_array,
			cudaStream_t &stream) const;

	int ExtendWithTriggerAsync(size_t query_seeds_size,
			uint32_t* query_seed_ids, uint32_t* query_seed_starts, size_t size,
			uint32_t* query_ids, uint32_t* query_concatenated_positions,
			uint32_t* database_positions, int* scores,
			uint32_t* d_query_seed_ids, uint32_t* d_query_seed_starts,
			size_t temp_storage_bytes, void* d_temp_storage,
			uint32_t* d_query_ids, uint32_t* d_query_concatenated_positions,
			uint32_t* d_database_positions, int* d_scores,
			cudaStream_t &stream) const;

	int ExtendOneSideAsync(size_t size, bool reverse, uint32_t* query_ids,
			uint32_t* query_concatenated_positions,
			uint32_t* database_positions, int* scores, uint32_t* d_query_ids,
			uint32_t* d_query_concatenated_positions,
			uint32_t* d_database_positions, int* d_scores,
			cudaStream_t &stream) const;

private:
	uint32_t number_letters_;
	AlphabetCoder::Code sequence_delimiter_;
	const uint32_t* d_query_seed_query_ids_;
	const uint32_t* d_query_seed_query_positions_;
	const packed_alphabet_code::PackedAlphabetCode *d_database_sequence_;
	const packed_alphabet_code::PackedAlphabetCode *d_concatenated_query_sequence_;
	const int *d_score_matrix_;
	const int *d_ungapped_extension_cutoffs_;
	const int *d_gapped_extension_triggers_;
};

#endif /* UNGAPPED_EXTENDER_GPU_H_ */
