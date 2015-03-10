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
#include "../src/packed_alphabet_code.h"

using namespace std;

class PackedAlphabetCodeTest: public ::testing::Test {
protected:
	virtual void SetUp() {
	}

	virtual void TearDown() {

	}
};

TEST_F(PackedAlphabetCodeTest, GetPackedAlphabetCodeSequenceLength) {
	size_t length = 10;
	EXPECT_EQ(
			(length + packed_alphabet_code::kNumberOfCodesInBlock - 1)
					/ packed_alphabet_code::kNumberOfCodesInBlock,
			packed_alphabet_code::GetPackedAlphabetCodeSequenceLength(length));
}

TEST_F(PackedAlphabetCodeTest, GetAlphabetCodeFoward) {
	DnaType type;
	AlphabetCoder coder(type);
	AlphabetCoder::Code delimiter_code = coder.GetMaxCode() + 1;
	string seq0("ACGTACGTACGT");
	vector<AlphabetCoder::Code> encoded_seq0(seq0.size() + 2, 0);
	encoded_seq0[0] = delimiter_code;
	coder.Encode(&seq0[0], seq0.size(), &encoded_seq0[1]);
	encoded_seq0[encoded_seq0.size() - 1] = delimiter_code;

	vector<packed_alphabet_code::PackedAlphabetCode> packed_code_seq0(
			packed_alphabet_code::GetPackedAlphabetCodeSequenceLength(
					encoded_seq0.size()) + 2);
	packed_alphabet_code::PackAlphabetCodeSequence(&encoded_seq0[0],
			encoded_seq0.size(), &packed_code_seq0[1]);
	const packed_alphabet_code::PackedAlphabetCode *code_blocks =
			&packed_code_seq0[1];

	int next_bolock_position = 0;
	uint32_t block_sift;
	uint32_t next_block_sift;
	packed_alphabet_code::PackedAlphabetCode temp_code_block;
	uint32_t original_start = 1;

	packed_alphabet_code::InitPackedAlphabetCode<cuda_common::kFoward>(
			code_blocks, original_start, &next_bolock_position, &block_sift,
			&next_block_sift, &temp_code_block);
	packed_alphabet_code::PackedAlphabetCode block =
			packed_alphabet_code::GetPackedAlphabetCode<cuda_common::kFoward>(
					code_blocks, block_sift, next_block_sift,
					&next_bolock_position, &temp_code_block);
	EXPECT_EQ(0,
			packed_alphabet_code::GetAlphabetCode<cuda_common::kFoward>(block,
					0));
	EXPECT_EQ(3,
			packed_alphabet_code::GetAlphabetCode<cuda_common::kFoward>(block,
					3));
}

TEST_F(PackedAlphabetCodeTest, GetAlphabetCodeReverse) {
	DnaType type;
	AlphabetCoder coder(type);
	AlphabetCoder::Code delimiter_code = coder.GetMaxCode() + 1;
	string seq0("ACGTACGTACGT");
	vector<AlphabetCoder::Code> encoded_seq0(seq0.size() + 2, 0);
	encoded_seq0[0] = delimiter_code;
	coder.Encode(&seq0[0], seq0.size(), &encoded_seq0[1]);
	encoded_seq0[encoded_seq0.size() - 1] = delimiter_code;

	vector<packed_alphabet_code::PackedAlphabetCode> packed_code_seq0(
			packed_alphabet_code::GetPackedAlphabetCodeSequenceLength(
					encoded_seq0.size()) + 2);
	packed_alphabet_code::PackAlphabetCodeSequence(&encoded_seq0[0],
			encoded_seq0.size(), &packed_code_seq0[1]);
	const packed_alphabet_code::PackedAlphabetCode *code_blocks =
			&packed_code_seq0[1];

	int next_bolock_position = 0;
	uint32_t block_sift;
	uint32_t next_block_sift;
	packed_alphabet_code::PackedAlphabetCode temp_code_block;
	uint32_t original_start = 12;

	packed_alphabet_code::InitPackedAlphabetCode<cuda_common::kReverse>(
			code_blocks, original_start, &next_bolock_position, &block_sift,
			&next_block_sift, &temp_code_block);

	packed_alphabet_code::PackedAlphabetCode block =
			packed_alphabet_code::GetPackedAlphabetCode<cuda_common::kReverse>(
					code_blocks, block_sift, next_block_sift,
					&next_bolock_position, &temp_code_block);
	EXPECT_EQ(-1, next_bolock_position);
	EXPECT_EQ(3,
			packed_alphabet_code::GetAlphabetCode<cuda_common::kReverse>(block,
					0));
	EXPECT_EQ(0,
			packed_alphabet_code::GetAlphabetCode<cuda_common::kReverse>(block,
					3));
}

