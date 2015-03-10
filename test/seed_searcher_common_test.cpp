/*
 * seed_search_common_test.cpp
 *
 *  Created on: 2014/05/13
 *      Author: shu
 */

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
#include "../src/seed_searcher_common.h"

using namespace std;

class SeedSearcherCommonCalculateSimilarityTest: public ::testing::Test {
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

TEST_F(SeedSearcherCommonCalculateSimilarityTest, Equal) {
	string sequence0 = "GMEHQS";
	vector<AlphabetCoder::Code> coded_sequence0(sequence0.length());
	coder_.Encode(&sequence0[0], sequence0.length(), &coded_sequence0[0]);
	string sequence1 = "GMEHQS";
	vector<AlphabetCoder::Code> coded_sequence1(sequence1.length());
	coder_.Encode(&sequence1[0], sequence1.length(), &coded_sequence1[0]);
	size_t center_position = sequence0.length() / 2;
	SeedSearcherCommon::Distance s = SeedSearcherCommon::CalculateDistance(
			&coded_sequence0[center_position],
			&coded_sequence1[center_position], sequence0.length(),
			reduced_code_map_);
	EXPECT_EQ(0, s);
}

TEST_F(SeedSearcherCommonCalculateSimilarityTest, 1Mismatch) {
	string sequence0 = "IMEHQS";
	vector<AlphabetCoder::Code> coded_sequence0(sequence0.length());
	coder_.Encode(&sequence0[0], sequence0.length(), &coded_sequence0[0]);
	string sequence1 = "LMEHQA";
	vector<AlphabetCoder::Code> coded_sequence1(sequence1.length());
	coder_.Encode(&sequence1[0], sequence1.length(), &coded_sequence1[0]);
	size_t center_position = sequence0.length() / 2;
	SeedSearcherCommon::Distance s = SeedSearcherCommon::CalculateDistance(
			&coded_sequence0[center_position],
			&coded_sequence1[center_position], sequence0.length(),
			reduced_code_map_);
	EXPECT_EQ(1, s);
}

class SeedSearcherCommonCalculateSeedPositionTest: public ::testing::Test {
protected:
	virtual void SetUp() {

	}

	virtual void TearDown() {

	}
};

TEST_F(SeedSearcherCommonCalculateSeedPositionTest, SeedLength2) {
	uint32_t hashed_sequence_length = 2;
	uint32_t hashed_sequence_start = 0;
	uint32_t seed_position = SeedSearcherCommon::CalculateSeedPosition(
			hashed_sequence_start, hashed_sequence_length);
	EXPECT_EQ(hashed_sequence_start, seed_position);

}
