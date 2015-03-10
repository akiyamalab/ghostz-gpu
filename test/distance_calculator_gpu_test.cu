#include <gtest/gtest.h>
#include <algorithm>
#include <string>
#include <sstream>
#include <stdint.h>
#include <fstream>
#include <vector>

#include "../src/protein_type.h"
#include "../src/alphabet_coder.h"
#include "../src/score_matrix_reader.h"
#include "../src/score_matrix.h"
#include "../src/reduced_alphabet_file_reader.h"
#include "../src/reduced_alphabet_coder.h"
#include "../src/reduced_alphabet_variable_hash_function.h"
#include "../src/distance_calculator_gpu.h"
#include "../src/cuda_common.h"
#include "../src/aligner_gpu_data.h"

#include <thrust/device_vector.h>
#include <thrust/copy.h>

using namespace std;

class DistanceCalculatorGpuTest: public ::testing::Test {
protected:
	virtual void SetUp() {
		ProteinType protein_type;
		const std::string reduced_alphabet = "A KR EDNQ C G H ILVM FYW P ST";

		std::istringstream in(reduced_alphabet);
		std::vector<std::string> alphabet_sets;
		ReducedAlphabetFileReader reduced_alphabet_reader;
		reduced_alphabet_reader.Read(in, alphabet_sets);
		coder_ = AlphabetCoder(protein_type);
		ReducedAlphabetCoder reduced_alphabet_coder(protein_type,
				alphabet_sets);
		AlphabetCoder::Code max_code = coder_.GetMaxRegularLetterCode();
		reduced_code_map_.resize(max_code + 1);
		for (AlphabetCoder::Code code = coder_.GetMinCode(); code <= max_code;
				++code) {
			char c = coder_.Decode(code);
			AlphabetCoder::Code reduced_code = reduced_alphabet_coder.Encode(c);
			reduced_code_map_[code] = reduced_code;
		}
		max_code_ = reduced_alphabet_coder.GetMaxRegularLetterCode();
		max_length_ = 2;

	}

	virtual void TearDown() {

	}
	AlphabetCoder::Code max_code_;
	uint32_t max_length_;
	int score_threshold_;
	std::vector<AlphabetCoder::Code> reduced_code_map_;
	std::vector<int> code_score_;
	AlphabetCoder coder_;
};

TEST_F(DistanceCalculatorGpuTest, CalculateDistances) {
	size_t subsequence_length = 4;
	string sequence0 = "IMEHQSIMEHQS";
	vector<AlphabetCoder::Code> coded_sequence0(sequence0.length());
	coder_.Encode(&sequence0[0], sequence0.length(), &coded_sequence0[0]);
	string sequence1 = "LMEHQALMEHQA";
	vector<AlphabetCoder::Code> coded_sequence1(sequence1.length());
	coder_.Encode(&sequence1[0], sequence1.length(), &coded_sequence1[0]);

	AlignerGpuData gpu_data;
	gpu_data.SetGpuQueriesSequence(&coded_sequence0[0], coded_sequence0.size());
	gpu_data.SetGpuDatabaseSequence(&coded_sequence1[0], coded_sequence1.size());
	gpu_data.SetGpuReducedCodeMap(&reduced_code_map_[0], reduced_code_map_.size());

	thrust::host_vector<uint32_t,
			thrust::cuda::experimental::pinned_allocator<uint32_t> > sequence0_positions;
	thrust::host_vector<uint32_t,
			thrust::cuda::experimental::pinned_allocator<uint32_t> > sequence1_positions;
	thrust::host_vector<DistanceCalculator::Distance,
			thrust::cuda::experimental::pinned_allocator<
					DistanceCalculator::Distance> > distances(3, 0);
	sequence0_positions.push_back(2);
	sequence1_positions.push_back(2);
	sequence0_positions.push_back(6);
	sequence1_positions.push_back(6);
	sequence0_positions.push_back(10);
	sequence1_positions.push_back(10);

	DistanceCalculatorGpu distance_calculator;
	cudaStream_t stream;
	thrust::device_vector<uint32_t> d_sequence0_positions(
			sequence0_positions.size());
	thrust::device_vector<uint32_t> d_sequence1_positions(
			sequence1_positions.size());
	thrust::device_vector<DistanceCalculatorGpu::Distance> d_distances(sequence0_positions.size());

	cudaStreamCreate(&stream);
	distance_calculator.SetQueries(gpu_data.GetGpuQueriesSequence());
	distance_calculator.SetDatabase(
			gpu_data.GetGpuDatabaseSequence());
	distance_calculator.SetReducedCodeMap(gpu_data.GetGpuReducedCodeMap());
	distance_calculator.SetSubsequenceLength(subsequence_length);
	distance_calculator.CalculateDistancesAsync(3, &sequence0_positions[0],
			&sequence1_positions[0], &distances[0],
			thrust::raw_pointer_cast(d_sequence0_positions.data()),
			thrust::raw_pointer_cast(d_sequence1_positions.data()),
			thrust::raw_pointer_cast(d_distances.data()), stream);
	cudaStreamSynchronize(stream);
	cudaStreamDestroy(stream);
	EXPECT_EQ(0, distances[0]);
	EXPECT_EQ(1, distances[1]);
	EXPECT_EQ(1, distances[2]);
}

