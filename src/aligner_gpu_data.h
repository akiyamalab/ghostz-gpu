/*
 * aligner_gpu_data.h
 *
 *  Created on: Jul 30, 2014
 *      Author: shu
 */

#ifndef ALIGNER_GPU_DATA_H_
#define ALIGNER_GPU_DATA_H_

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/system/cuda/experimental/pinned_allocator.h>
#include <thrust/copy.h>
#include <vector>
#include "cuda_common.h"
#include "alphabet_coder.h"
#include "packed_alphabet_code.h"

class AlignerGpuData {
public:
	AlignerGpuData();
	AlignerGpuData(int gpu_id);
	virtual ~AlignerGpuData();

	void SetGpuId(int gpu_id);
	int GetGpuId();

	void SetGpuQueryIds(uint32_t *query_ids, size_t length);
	uint32_t* GetGpuQuerySeedQueryIds();

	void SetGpuConcatenatedQuerySequencePositions(
			uint32_t *concatenated_query_sequence_positions, size_t length);
	uint32_t* GetGpuQuerySeedQueryPositions();

	void SetGpuDatabaseSequence(const AlphabetCoder::Code* sequence,
			size_t length);
	packed_alphabet_code::PackedAlphabetCode* GetGpuDatabaseSequence();
	void SetGpuQueriesSequence(const AlphabetCoder::Code* sequence,
			size_t length);
	packed_alphabet_code::PackedAlphabetCode* GetGpuQueriesSequence();

	void SetGpuReducedCodeMap(const AlphabetCoder::Code* map, size_t length);
	AlphabetCoder::Code* GetGpuReducedCodeMap();

	void SetGpuUngappedExtensionCutoffs(const int* cutoffs, size_t number);
	int* GetGpuUngappedExtensionCutoffs();

	void SetGpuGappedExtensionTriggers(const int* triggers, size_t number);
	int* GetGpuGappedExtensionTriggers();

	void SetGpuScoreMatrix(const int* score_matrix, size_t score_matrix_size);
	int* GetGpuScoreMatrix();

private:
	int gpu_id_;
	std::vector<packed_alphabet_code::PackedAlphabetCode> temp_packed_alphabet_code_sequence_;
	thrust::device_vector<uint32_t> d_query_ids_;
	thrust::device_vector<uint32_t> d_concatenated_query_sequence_positions_;
	thrust::device_vector<packed_alphabet_code::PackedAlphabetCode> d_concatenated_query_sequence_;
	thrust::device_vector<packed_alphabet_code::PackedAlphabetCode> d_database_sequence_;
	thrust::device_vector<AlphabetCoder::Code> d_reduced_code_map_;
	thrust::device_vector<int> d_ungapped_extension_cutoffs_;
	thrust::device_vector<int> d_gapped_extension_triggers_;
	thrust::device_vector<int> d_score_matrix_;

	void SetPackedAlphabetCodeSequenceToGpu(size_t length,
			const AlphabetCoder::Code* sequence,
			thrust::device_vector<packed_alphabet_code::PackedAlphabetCode> &d_packed_alphabet_code_sequence);
};

#endif /* ALIGNER_GPU_DATA_H_ */
