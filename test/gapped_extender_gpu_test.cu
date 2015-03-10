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
#include "../src/gapped_extender_gpu.h"
#include "../src/cuda_common.h"
#include "../src/aligner_gpu_data.h"

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>

using namespace std;

class GappedExtenderGpuTest: public ::testing::Test {
protected:
	virtual void SetUp() {
	}

	virtual void TearDown() {

	}
};

TEST_F(GappedExtenderGpuTest, ExtendShortSequencesForwardWithMinXDrop) {
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

	vector<int> sequence0_cutoff(1, 1);

	AlignerGpuData gpu_data;
	gpu_data.SetGpuQueriesSequence(&encoded_seq0[0], encoded_seq0.size());
	gpu_data.SetGpuDatabaseSequence(&encoded_seq1[0], encoded_seq1.size());
	gpu_data.SetGpuScoreMatrix(score_matrix.GetMatrix(),
			score_matrix.GetNumberLetters());
	gpu_data.SetGpuUngappedExtensionCutoffs(&sequence0_cutoff[0],
			sequence0_cutoff.size());

	GappedExtenderGpu e;
	e.SetQueries(delimiter_code, gpu_data.GetGpuQueriesSequence());
	e.SetDatabase(gpu_data.GetGpuDatabaseSequence());
	e.SetScoreParameters(gpu_data.GetGpuScoreMatrix(),
			score_matrix.GetNumberLetters(), 0, 1, 0);

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

	e.ExtendOneSideScoreOnlyAsync(number_extensions, false,
			&sequence0_positions[0], &sequence1_positions[0], &scores[0],
			thrust::raw_pointer_cast(d_sequence0_positions.data()),
			thrust::raw_pointer_cast(d_sequence1_positions.data()),
			thrust::raw_pointer_cast(d_scores.data()), stream);

	cudaStreamSynchronize(stream);
	cudaStreamDestroy(stream);

	EXPECT_EQ(8, scores[0]);
	EXPECT_EQ(4, sequence0_positions[0]);
	EXPECT_EQ(4, sequence1_positions[0]);
}

TEST_F(GappedExtenderGpuTest, ExtendShortSequencesForward) {
	DnaType type;
	AlphabetCoder coder(type);
	AlphabetCoder::Code delimiter_code = coder.GetMaxCode() + 1;
	ScoreMatrix score_matrix("test_matrix", 4, 2, -1);

	string seq0("AGTCAC");
	vector<AlphabetCoder::Code> encoded_seq0(seq0.size() + 2);
	encoded_seq0[0] = delimiter_code;
	coder.Encode(&seq0[0], seq0.size(), &encoded_seq0[1]);
	encoded_seq0[encoded_seq0.size() - 1] = delimiter_code;

	string seq1("AGCAG");
	vector<AlphabetCoder::Code> encoded_seq1(seq1.size() + 2);
	encoded_seq1[0] = delimiter_code;
	coder.Encode(&seq1[0], seq1.size(), &encoded_seq1[1]);
	encoded_seq1[encoded_seq1.size() - 1] = delimiter_code;

	vector<int> sequence0_cutoff(1, 1);

	AlignerGpuData gpu_data;
	gpu_data.SetGpuQueriesSequence(&encoded_seq0[0], encoded_seq0.size());
	gpu_data.SetGpuDatabaseSequence(&encoded_seq1[0], encoded_seq1.size());
	gpu_data.SetGpuScoreMatrix(score_matrix.GetMatrix(),
			score_matrix.GetNumberLetters());
	gpu_data.SetGpuUngappedExtensionCutoffs(&sequence0_cutoff[0],
			sequence0_cutoff.size());

	GappedExtenderGpu e;
	e.SetQueries(delimiter_code, gpu_data.GetGpuQueriesSequence());
	e.SetDatabase(gpu_data.GetGpuDatabaseSequence());
	e.SetScoreParameters(gpu_data.GetGpuScoreMatrix(),
			score_matrix.GetNumberLetters(), 0, 1, INT_MAX);

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

	e.ExtendOneSideScoreOnlyAsync(number_extensions, false,
			&sequence0_positions[0], &sequence1_positions[0], &scores[0],
			thrust::raw_pointer_cast(d_sequence0_positions.data()),
			thrust::raw_pointer_cast(d_sequence1_positions.data()),
			thrust::raw_pointer_cast(d_scores.data()), stream);

	cudaStreamSynchronize(stream);
	cudaStreamDestroy(stream);

	EXPECT_EQ(7, scores[0]);
	EXPECT_EQ(5, sequence0_positions[0]);
	EXPECT_EQ(4, sequence1_positions[0]);
}

TEST_F(GappedExtenderGpuTest, ExtendShortSequencesReverse) {
	DnaType type;
	AlphabetCoder coder(type);
	AlphabetCoder::Code delimiter_code = coder.GetMaxCode() + 1;
	ScoreMatrix score_matrix("test_matrix", 4, 2, -1);

	string seq0("AGTCAC");
	vector<AlphabetCoder::Code> encoded_seq0(seq0.size() + 2);
	encoded_seq0[0] = delimiter_code;
	coder.Encode(&seq0[0], seq0.size(), &encoded_seq0[1]);
	encoded_seq0[encoded_seq0.size() - 1] = delimiter_code;

	string seq1("AGCAG");
	vector<AlphabetCoder::Code> encoded_seq1(seq1.size() + 2);
	encoded_seq1[0] = delimiter_code;
	coder.Encode(&seq1[0], seq1.size(), &encoded_seq1[1]);
	encoded_seq1[encoded_seq1.size() - 1] = delimiter_code;

	vector<int> sequence0_cutoff(1, 1);
	AlignerGpuData gpu_data;
	gpu_data.SetGpuQueriesSequence(&encoded_seq0[0], encoded_seq0.size());
	gpu_data.SetGpuDatabaseSequence(&encoded_seq1[0], encoded_seq1.size());
	gpu_data.SetGpuScoreMatrix(score_matrix.GetMatrix(),
			score_matrix.GetNumberLetters());
	gpu_data.SetGpuUngappedExtensionCutoffs(&sequence0_cutoff[0],
			sequence0_cutoff.size());

	GappedExtenderGpu e;
	e.SetQueries(delimiter_code, gpu_data.GetGpuQueriesSequence());
	e.SetDatabase(gpu_data.GetGpuDatabaseSequence());
	e.SetScoreParameters(gpu_data.GetGpuScoreMatrix(),
			score_matrix.GetNumberLetters(), 0, 1, INT_MAX);

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
	sequence0_positions[0] = encoded_seq0.size() - 3;
	sequence1_positions[0] = encoded_seq1.size() - 3;
	scores[0] = 0;

	cudaStream_t stream;
	thrust::device_vector<uint32_t> d_ids(ids.size());
	thrust::device_vector<uint32_t> d_sequence0_positions(
			sequence0_positions.size());
	thrust::device_vector<uint32_t> d_sequence1_positions(
			sequence1_positions.size());
	thrust::device_vector<int> d_scores(sequence1_positions.size());
	cudaStreamCreate(&stream);

	e.ExtendOneSideScoreOnlyAsync(number_extensions, true,
			&sequence0_positions[0], &sequence1_positions[0], &scores[0],
			thrust::raw_pointer_cast(d_sequence0_positions.data()),
			thrust::raw_pointer_cast(d_sequence1_positions.data()),
			thrust::raw_pointer_cast(d_scores.data()), stream);

	cudaStreamSynchronize(stream);
	cudaStreamDestroy(stream);

	EXPECT_EQ(7, scores[0]);
	EXPECT_EQ(1, sequence0_positions[0]);
	EXPECT_EQ(1, sequence1_positions[0]);
}

