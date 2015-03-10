/*
 * aligner_gpu_data.cpp
 *
 *  Created on: Jul 30, 2014
 *      Author: shu
 */

#include "aligner_gpu_data.h"
#include <vector>

using namespace std;

AlignerGpuData::AlignerGpuData() :
		gpu_id_(0) {
}

AlignerGpuData::AlignerGpuData(int gpu_id) :
		gpu_id_(gpu_id) {
	SetGpuId(gpu_id);
}

AlignerGpuData::~AlignerGpuData() {
}

void AlignerGpuData::SetGpuId(int gpu_id) {
	gpu_id_ = gpu_id;
	CUDA_CHECK_RETURN(cudaSetDevice(gpu_id_));
}

int AlignerGpuData::GetGpuId() {
	return gpu_id_;
}

void AlignerGpuData::SetGpuQueryIds(uint32_t *query_ids, size_t length) {
	d_query_ids_.resize(length);
	thrust::copy(query_ids, query_ids + length, d_query_ids_.begin());
}

uint32_t* AlignerGpuData::GetGpuQuerySeedQueryIds() {
	return thrust::raw_pointer_cast(d_query_ids_.data());
}

void AlignerGpuData::SetGpuConcatenatedQuerySequencePositions(
		uint32_t *concatenated_query_sequence_positions, size_t length) {
	d_concatenated_query_sequence_positions_.resize(length);
	thrust::copy(concatenated_query_sequence_positions,
			concatenated_query_sequence_positions + length,
			d_concatenated_query_sequence_positions_.begin());
}

uint32_t* AlignerGpuData::GetGpuQuerySeedQueryPositions() {
	return thrust::raw_pointer_cast(
			d_concatenated_query_sequence_positions_.data());
}

void AlignerGpuData::SetGpuDatabaseSequence(const AlphabetCoder::Code* sequence,
		size_t length) {
	SetPackedAlphabetCodeSequenceToGpu(length, sequence, d_database_sequence_);
	return;
}
packed_alphabet_code::PackedAlphabetCode* AlignerGpuData::GetGpuDatabaseSequence() {
	return thrust::raw_pointer_cast(d_database_sequence_.data());
}

void AlignerGpuData::SetGpuQueriesSequence(const AlphabetCoder::Code* sequence,
		size_t length) {
	SetPackedAlphabetCodeSequenceToGpu(length, sequence,
			d_concatenated_query_sequence_);
	return;
}

packed_alphabet_code::PackedAlphabetCode* AlignerGpuData::GetGpuQueriesSequence() {
	return thrust::raw_pointer_cast(d_concatenated_query_sequence_.data());
}

void AlignerGpuData::SetGpuUngappedExtensionCutoffs(const int* cutoffs,
		size_t number) {
	d_ungapped_extension_cutoffs_.resize(number);
	thrust::copy(cutoffs, cutoffs + number,
			d_ungapped_extension_cutoffs_.begin());
	return;
}

int* AlignerGpuData::GetGpuUngappedExtensionCutoffs() {
	return thrust::raw_pointer_cast(d_ungapped_extension_cutoffs_.data());
}

void AlignerGpuData::SetGpuGappedExtensionTriggers(const int* triggers,
		size_t number) {
	d_gapped_extension_triggers_.resize(number);
	thrust::copy(triggers, triggers + number,
			d_gapped_extension_triggers_.begin());
}

int* AlignerGpuData::GetGpuGappedExtensionTriggers() {
	return thrust::raw_pointer_cast(d_gapped_extension_triggers_.data());
}

void AlignerGpuData::SetGpuReducedCodeMap(const AlphabetCoder::Code* map,
		size_t length) {
	d_reduced_code_map_.resize(length);
	thrust::copy(map, map + length, d_reduced_code_map_.begin());
}
AlphabetCoder::Code* AlignerGpuData::GetGpuReducedCodeMap() {
	return thrust::raw_pointer_cast(d_reduced_code_map_.data());
}

void AlignerGpuData::SetGpuScoreMatrix(const int* score_matrix,
		size_t number_letters) {
	size_t number_letters_for_gpu = number_letters + 1;
	vector<int> score_matrix_for_gpu(
			number_letters_for_gpu * number_letters_for_gpu, 0);
	for (size_t c0 = 0; c0 < number_letters; ++c0) {
		for (size_t c1 = 0; c1 < number_letters; ++c1) {
			score_matrix_for_gpu[c0 * number_letters_for_gpu + c1] =
					score_matrix[c0 * number_letters + c1];
		}
	}
	d_score_matrix_.resize(number_letters_for_gpu * number_letters_for_gpu);
	thrust::copy(score_matrix_for_gpu.begin(), score_matrix_for_gpu.end(),
			d_score_matrix_.begin());
}
int* AlignerGpuData::GetGpuScoreMatrix() {
	return thrust::raw_pointer_cast(d_score_matrix_.data());
}

void AlignerGpuData::SetPackedAlphabetCodeSequenceToGpu(size_t length,
		const AlphabetCoder::Code* sequence,
		thrust::device_vector<packed_alphabet_code::PackedAlphabetCode> &d_packed_alphabet_code_sequence) {
	size_t packed_sequence_length =
			packed_alphabet_code::GetPackedAlphabetCodeSequenceLength(length);
	if (temp_packed_alphabet_code_sequence_.size() < packed_sequence_length) {
		temp_packed_alphabet_code_sequence_.resize(packed_sequence_length);
	}
	packed_alphabet_code::PackAlphabetCodeSequence(sequence, length,
			&temp_packed_alphabet_code_sequence_[0]);
	size_t new_d_packed_alphabet_code_sequence_size = packed_sequence_length
			+ (2 * cuda_common::kMaxLoadLength);
	if (d_packed_alphabet_code_sequence.size() < new_d_packed_alphabet_code_sequence_size) {
		d_packed_alphabet_code_sequence.resize(new_d_packed_alphabet_code_sequence_size);
	}
	thrust::copy(temp_packed_alphabet_code_sequence_.begin(),
			temp_packed_alphabet_code_sequence_.begin()
					+ packed_sequence_length,
			d_packed_alphabet_code_sequence.begin()
					+ cuda_common::kMaxLoadLength);
}
