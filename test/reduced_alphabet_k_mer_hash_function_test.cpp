/*
 * reduced_alphabet_k_mer_hash_function_test.cpp
 *
 *  Created on: 2013/04/22
 *      Author: shu
 */

#include <gtest/gtest.h>
#include <string>
#include <sstream>
#include <stdint.h>
#include <fstream>
#include <vector>
#include "../src/protein_type.h"
#include "../src/alphabet_coder.h"
#include "../src/reduced_alphabet_file_reader.h"
#include "../src/reduced_alphabet_coder.h"
#include "../src/reduced_alphabet_k_mer_hash_function.h"

using namespace std;

class ReducedAlphabetKMerHashFunctionTest: public ::testing::Test {
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
		k_mer_length_ = 2;
		hash_function_ = ReducedAlphabetKMerHashFunction(max_code_,
				k_mer_length_, reduced_code_map_);
	}

	virtual void TearDown() {

	}
	AlphabetCoder::Code max_code_;
	uint32_t k_mer_length_;
	std::vector<AlphabetCoder::Code> reduced_code_map_;
	AlphabetCoder coder_;
	ReducedAlphabetKMerHashFunction hash_function_;
};

TEST_F(ReducedAlphabetKMerHashFunctionTest, VoidConstructor) {
	ReducedAlphabetKMerHashFunction hash_function;
	EXPECT_FALSE(hash_function.CanCalculateHash());
}

TEST_F(ReducedAlphabetKMerHashFunctionTest, GetKMerLength) {
	EXPECT_EQ(k_mer_length_, hash_function_.GetKMerLength());
}

TEST_F(ReducedAlphabetKMerHashFunctionTest, GetMaxHash) {
	uint32_t shift_size = 0;
	for (AlphabetCoder::Code c = 1; c < max_code_; c <<= 1, ++shift_size) {
	}

	ReducedAlphabetKMerHashFunction::Hash max_hash = 0;
	for (uint32_t i = 0; i < k_mer_length_; ++i) {
		max_hash <<= shift_size;
		max_hash += max_code_;
	}
	EXPECT_EQ(max_hash, hash_function_.GetMaxHash());
}

TEST_F(ReducedAlphabetKMerHashFunctionTest, CalculateHash) {
  uint32_t shift_size = 0;
  for (AlphabetCoder::Code c = 1; c < max_code_; c <<= 1, ++shift_size) {
  }
  string k_mer_string = "GMEHQS";
  vector<AlphabetCoder::Code> k_mer_coded_string(k_mer_string.length());
  coder_.Encode(&k_mer_string[0], k_mer_string.length(), &k_mer_coded_string[0]);

  ReducedAlphabetKMerHashFunction::Hash expected_hash = 0;
  for (uint32_t i = 0; i < k_mer_length_; ++i) {
    expected_hash <<= shift_size;
    expected_hash += reduced_code_map_[k_mer_coded_string[i]];
  }

  ReducedAlphabetKMerHashFunction::Hash hash = 0;
  EXPECT_EQ(0, hash_function_.CalculateHash(&k_mer_coded_string[0], &hash));
  EXPECT_EQ(expected_hash, hash);
  cout << hash << endl;
}

TEST_F(ReducedAlphabetKMerHashFunctionTest, SaveAndLoad) {
	string filename = "hash_function_test";
	ofstream ofs(filename.c_str());
	hash_function_.Save(ofs);
	ofs.close();
	ifstream ifs(filename.c_str());
	ReducedAlphabetKMerHashFunction loaded_hash_function(ifs);
	string k_mer_string = "AK";
	vector<AlphabetCoder::Code> k_mer_coded_string(k_mer_string.length());
	coder_.Encode(&k_mer_string[0], k_mer_string.length(),
			&k_mer_coded_string[0]);

	ReducedAlphabetKMerHashFunction::Hash expected_hash = 0;
	int ret = hash_function_.CalculateHash(&k_mer_coded_string[0],
			&expected_hash);

	ReducedAlphabetKMerHashFunction::Hash hash = 0;
	EXPECT_EQ(ret,
			loaded_hash_function.CalculateHash(&k_mer_coded_string[0], &hash));
	EXPECT_EQ(expected_hash, hash);
}

