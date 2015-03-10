/*
 * gapped_extender_test.cpp
 *
 *  Created on: 2012/11/20
 *      Author: shu
 */

#include <gtest/gtest.h>
#include <string>
#include <stdint.h>
#include <fstream>
#include <limits.h>
#include "../src/score_matrix.h"
#include "../src/alphabet_coder.h"
#include "../src/sequence_type.h"
#include "../src/protein_type.h"
#include "../src/dna_type.h"
#include "../src/score_matrix.h"
#include "../src/edit_blocks.h"
#include "../src/ungapped_extender_gpu.h"
#include "../src/cuda_common.h"
#include "../src/aligner_gpu_data.h"

#include <thrust/host_vector.h>
#include <thrust/system/cuda/experimental/pinned_allocator.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>

using namespace std;

class UngappedExtenderGpuTest: public ::testing::Test {
protected:
	virtual void SetUp() {
	}

	virtual void TearDown() {

	}
};

TEST_F(UngappedExtenderGpuTest, ExtendShortSequencesForwardWithMinXDrop) {
	DnaType type;
	AlphabetCoder coder(type);
	AlphabetCoder::Code delimiter_code = coder.GetMaxCode() + 1;
	ScoreMatrix score_matrix("test_matrix", 4, 2, -1);

	string seq0("AGCAC");
	vector<AlphabetCoder::Code> encoded_seq0(seq0.size() + 2);
	encoded_seq0[0] = delimiter_code;
	coder.Encode(&seq0[0], seq0.size(), &encoded_seq0[1]);
	encoded_seq0[encoded_seq0.size() - 1] = delimiter_code;

	string seq1("AGAAG");
	vector<AlphabetCoder::Code> encoded_seq1(seq1.size() + 2);
	encoded_seq1[0] = delimiter_code;
	coder.Encode(&seq1[0], seq1.size(), &encoded_seq1[1]);
	encoded_seq1[encoded_seq1.size() - 1] = delimiter_code;
	vector<int> sequence0_cutoff(1, 1);
	vector<int> sequence0_trigger(1, 3);
	AlignerGpuData gpu_data;
	gpu_data.SetGpuQueriesSequence(&encoded_seq0[0], encoded_seq0.size());
	gpu_data.SetGpuDatabaseSequence(&encoded_seq1[0], encoded_seq1.size());
	gpu_data.SetGpuScoreMatrix(score_matrix.GetMatrix(),
			score_matrix.GetNumberLetters());
	gpu_data.SetGpuUngappedExtensionCutoffs(&sequence0_cutoff[0],
			sequence0_cutoff.size());
	gpu_data.SetGpuGappedExtensionTriggers(&sequence0_trigger[0],
			sequence0_trigger.size());

	UngappedExtenderGpu e;

	e.SetQueries(delimiter_code, gpu_data.GetGpuQueriesSequence(),
			gpu_data.GetGpuUngappedExtensionCutoffs(),
			gpu_data.GetGpuGappedExtensionTriggers());
	e.SetDatabase(gpu_data.GetGpuDatabaseSequence());
	e.SetScoreMatrix(gpu_data.GetGpuScoreMatrix(),
			score_matrix.GetNumberLetters());

	size_t number_extensions = 1;
	thrust::host_vector<uint32_t,
			thrust::cuda::experimental::pinned_allocator<uint32_t> > ids(1, 0);
	thrust::host_vector<uint32_t,
			thrust::cuda::experimental::pinned_allocator<uint32_t> > sequence0_positions(
			number_extensions);
	thrust::host_vector<uint32_t,
			thrust::cuda::experimental::pinned_allocator<uint32_t> > sequence1_positions(
			number_extensions);
	thrust::host_vector<int, thrust::cuda::experimental::pinned_allocator<int> > scores(
			number_extensions, 0);
	sequence0_positions[0] = 1;
	sequence1_positions[0] = 1;
	scores[0] = 0;

	cudaStream_t stream;

	thrust::device_vector<uint32_t> d_ids(ids.size());
	thrust::device_vector<uint32_t> d_sequence0_positions(
			sequence0_positions.size());
	thrust::device_vector<uint32_t> d_sequence1_positions(
			sequence1_positions.size());
	thrust::device_vector<bool> d_scores(sequence1_positions.size());

	cudaStreamCreate(&stream);
	e.ExtendWithTriggerAsync(number_extensions, &ids[0], &sequence0_positions[0],
			&sequence1_positions[0], &scores[0],
			thrust::raw_pointer_cast(d_ids.data()),
			thrust::raw_pointer_cast(d_sequence0_positions.data()),
			thrust::raw_pointer_cast(d_sequence1_positions.data()),
			thrust::raw_pointer_cast(d_scores.data()), stream);
	cudaStreamSynchronize(stream);
	cudaStreamDestroy(stream);

	EXPECT_EQ(4, scores[0]);
}

TEST_F(UngappedExtenderGpuTest, ExtendShortSequencesForward) {
	DnaType type;
	AlphabetCoder coder(type);
	AlphabetCoder::Code delimiter_code = coder.GetMaxCode() + 1;
	ScoreMatrix score_matrix("test_matrix", 4, 2, -1);

	string seq0("AGCAC");
	vector<AlphabetCoder::Code> encoded_seq0(seq0.size() + 2);
	encoded_seq0[0] = delimiter_code;
	coder.Encode(&seq0[0], seq0.size(), &encoded_seq0[1]);
	encoded_seq0[encoded_seq0.size() - 1] = delimiter_code;

	string seq1("AGCAG");
	vector<AlphabetCoder::Code> encoded_seq1(seq1.size() + 2);
	encoded_seq1[0] = delimiter_code;
	coder.Encode(&seq1[0], seq1.size(), &encoded_seq1[1]);
	encoded_seq1[encoded_seq1.size() - 1] = delimiter_code;

	vector<int> sequence0_cutoff(1, 1024);
	vector<int> sequence0_trigger(1, 7);
	AlignerGpuData gpu_data;
	gpu_data.SetGpuQueriesSequence(&encoded_seq0[0], encoded_seq0.size());
	gpu_data.SetGpuDatabaseSequence(&encoded_seq1[0], encoded_seq1.size());
	gpu_data.SetGpuScoreMatrix(score_matrix.GetMatrix(),
			score_matrix.GetNumberLetters());
	gpu_data.SetGpuUngappedExtensionCutoffs(&sequence0_cutoff[0],
			sequence0_cutoff.size());
	gpu_data.SetGpuGappedExtensionTriggers(&sequence0_trigger[0],
			sequence0_trigger.size());

	UngappedExtenderGpu e;

	e.SetQueries(delimiter_code, gpu_data.GetGpuQueriesSequence(),
			gpu_data.GetGpuUngappedExtensionCutoffs(),
			gpu_data.GetGpuGappedExtensionTriggers());
	e.SetDatabase(gpu_data.GetGpuDatabaseSequence());
	e.SetScoreMatrix(gpu_data.GetGpuScoreMatrix(),
			score_matrix.GetNumberLetters());

	size_t number_extensions = 1;
	thrust::host_vector<uint32_t,
			thrust::cuda::experimental::pinned_allocator<uint32_t> > ids(1, 0);
	thrust::host_vector<uint32_t,
			thrust::cuda::experimental::pinned_allocator<uint32_t> > sequence0_positions(
			number_extensions);
	thrust::host_vector<uint32_t,
			thrust::cuda::experimental::pinned_allocator<uint32_t> > sequence1_positions(
			number_extensions);
	thrust::host_vector<int, thrust::cuda::experimental::pinned_allocator<int> > scores(
			number_extensions, 0);
	sequence0_positions[0] = 1;
	sequence1_positions[0] = 1;
	scores[0] = 0;

	cudaStream_t stream;
	thrust::device_vector<uint32_t> d_ids(ids.size());
	thrust::device_vector<uint32_t> d_sequence0_positions(
			sequence0_positions.size());
	thrust::device_vector<uint32_t> d_sequence1_positions(
			sequence1_positions.size());
	thrust::device_vector<int> d_scores(sequence1_positions.size());
	cudaStreamCreate(&stream);
	e.ExtendWithTriggerAsync(number_extensions, &ids[0], &sequence0_positions[0],
			&sequence1_positions[0], &scores[0],
			thrust::raw_pointer_cast(d_ids.data()),
			thrust::raw_pointer_cast(d_sequence0_positions.data()),
			thrust::raw_pointer_cast(d_sequence1_positions.data()),
			thrust::raw_pointer_cast(d_scores.data()), stream);
	cudaStreamSynchronize(stream);
	cudaStreamDestroy(stream);

	EXPECT_EQ(8, scores[0]);
}

TEST_F(UngappedExtenderGpuTest, ExtendShortSequencesReverse) {
	DnaType type;
	AlphabetCoder coder(type);
	AlphabetCoder::Code delimiter_code = coder.GetMaxCode() + 1;
	ScoreMatrix score_matrix("test_matrix", 4, 2, -1);

	string seq0("CAGCA");
	vector<AlphabetCoder::Code> encoded_seq0(seq0.size() + 2);
	encoded_seq0[0] = delimiter_code;
	coder.Encode(&seq0[0], seq0.size(), &encoded_seq0[1]);
	encoded_seq0[encoded_seq0.size() - 1] = delimiter_code;

	string seq1("GAGCA");
	vector<AlphabetCoder::Code> encoded_seq1(seq1.size() + 2);
	encoded_seq1[0] = delimiter_code;
	coder.Encode(&seq1[0], seq1.size(), &encoded_seq1[1]);
	encoded_seq1[encoded_seq1.size() - 1] = delimiter_code;

	vector<int> sequence0_cutoff(1, 1024);
	vector<int> sequence0_trigger(1, 7);
	AlignerGpuData gpu_data;
	gpu_data.SetGpuQueriesSequence(&encoded_seq0[0], encoded_seq0.size());
	gpu_data.SetGpuDatabaseSequence(&encoded_seq1[0], encoded_seq1.size());
	gpu_data.SetGpuScoreMatrix(score_matrix.GetMatrix(),
			score_matrix.GetNumberLetters());
	gpu_data.SetGpuUngappedExtensionCutoffs(&sequence0_cutoff[0],
			sequence0_cutoff.size());
	gpu_data.SetGpuGappedExtensionTriggers(&sequence0_trigger[0],
			sequence0_trigger.size());

	UngappedExtenderGpu e;

	e.SetQueries(delimiter_code, gpu_data.GetGpuQueriesSequence(),
			gpu_data.GetGpuUngappedExtensionCutoffs(),
			gpu_data.GetGpuGappedExtensionTriggers());
	e.SetDatabase(gpu_data.GetGpuDatabaseSequence());
	e.SetScoreMatrix(gpu_data.GetGpuScoreMatrix(),
			score_matrix.GetNumberLetters());

	size_t number_extensions = 1;
	thrust::host_vector<uint32_t,
			thrust::cuda::experimental::pinned_allocator<uint32_t> > ids(1, 0);
	thrust::host_vector<uint32_t,
			thrust::cuda::experimental::pinned_allocator<uint32_t> > sequence0_positions(
			number_extensions);
	thrust::host_vector<uint32_t,
			thrust::cuda::experimental::pinned_allocator<uint32_t> > sequence1_positions(
			number_extensions);
	thrust::host_vector<int, thrust::cuda::experimental::pinned_allocator<int> > scores(
			number_extensions, 0);
	sequence0_positions[0] = encoded_seq0.size() - 2;
	sequence1_positions[0] = encoded_seq1.size() - 2;
	scores[0] = 0;

	cudaStream_t stream;
	thrust::device_vector<uint32_t> d_ids(ids.size());
	thrust::device_vector<uint32_t> d_sequence0_positions(
			sequence0_positions.size());
	thrust::device_vector<uint32_t> d_sequence1_positions(
			sequence1_positions.size());
	thrust::device_vector<int> d_scores(sequence1_positions.size());
	cudaStreamCreate(&stream);
	e.ExtendWithTriggerAsync(number_extensions, &ids[0], &sequence0_positions[0],
			&sequence1_positions[0], &scores[0],
			thrust::raw_pointer_cast(d_ids.data()),
			thrust::raw_pointer_cast(d_sequence0_positions.data()),
			thrust::raw_pointer_cast(d_sequence1_positions.data()),
			thrust::raw_pointer_cast(d_scores.data()), stream);
	cudaStreamSynchronize(stream);
	cudaStreamDestroy(stream);

	EXPECT_EQ(8, scores[0]);
}

