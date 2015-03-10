/*
 * one_mismatch_hash_generator_test.cpp
 *
 *  Created on: 2013/05/03
 *      Author: shu
 */

#include <gtest/gtest.h>
#include <string>
#include <stdint.h>
#include <fstream>
#include <stdint.h>
#include "../src/protein_type.h"
#include "../src/reduced_alphabet_file_reader.h"
#include "../src/reduced_alphabet_coder.h"
#include "../src/reduced_alphabet_k_mer_hash_function.h"
#include "../src/alphabet_coder.h"
#include "../src/one_mismatch_hash_generator.h"

using namespace std;

class OneMismatchHashGeneratorTest: public ::testing::Test {
protected:
  virtual void SetUp() {
    ProteinType protein_type;
    const std::string reduced_alphabet = "A KR EDNQ C G H ILVM FYW P ST";
    std::istringstream in(reduced_alphabet);
    std::vector<std::string> alphabet_sets;
    ReducedAlphabetFileReader reduced_alphabet_reader;
    reduced_alphabet_reader.Read(in, alphabet_sets);
    coder_ = AlphabetCoder(protein_type);
    ReducedAlphabetCoder reduced_alphabet_coder(protein_type, alphabet_sets);
    AlphabetCoder::Code max_code = coder_.GetMaxRegularLetterCode();
    reduced_code_map_.resize(max_code + 1);
    for (AlphabetCoder::Code code = coder_.GetMinCode(); code <= max_code; ++code) {
      char c = coder_.Decode(code);
      AlphabetCoder::Code reduced_code = reduced_alphabet_coder.Encode(c);
      reduced_code_map_[code] = reduced_code;
    }
    max_code_ = reduced_alphabet_coder.GetMaxRegularLetterCode();
    k_mer_length_ = 2;
    hash_function_ = ReducedAlphabetKMerHashFunction(max_code_, k_mer_length_, reduced_code_map_);
  }

  virtual void TearDown() {

  }
  AlphabetCoder::Code max_code_;
  uint32_t k_mer_length_;
  std::vector<AlphabetCoder::Code> reduced_code_map_;
  AlphabetCoder coder_;
  ReducedAlphabetKMerHashFunction hash_function_;
};

TEST_F(OneMismatchHashGeneratorTest, GenerateOneMismatchHash) {
  ReducedAlphabetKMerHashFunction hash_function;
  OneMismatchHashGenerator one_mismatch_hash_generator;
  ReducedAlphabetKMerHashFunction::Hash basic_hash;
  string k_mer_string = "AK";
  vector<AlphabetCoder::Code> k_mer_coded_string(k_mer_string.length());
  coder_.Encode(&k_mer_string[0], k_mer_string.length(), &k_mer_coded_string[0]);
  hash_function.CalculateHash(&k_mer_coded_string[0], &basic_hash);
  vector<ReducedAlphabetKMerHashFunction::Hash> onemismatch_hashs;
  one_mismatch_hash_generator.GenerateOneMismatchHash(basic_hash, k_mer_length_, max_code_,
      hash_function_.GetBitShiftSize(), &onemismatch_hashs);

  EXPECT_EQ(k_mer_length_*max_code_, onemismatch_hashs.size());
  uint32_t bit_mask = (1 << hash_function_.GetBitShiftSize()) - 1;
  for (size_t i = 0; i < onemismatch_hashs.size(); ++i) {
    uint32_t mismatch_count = 0;
    ReducedAlphabetKMerHashFunction::Hash generated_hash = onemismatch_hashs[i];
    ReducedAlphabetKMerHashFunction::Hash tmp_basic_hash = basic_hash;
    for (uint32_t hash_i = 0; hash_i < k_mer_length_;
        ++hash_i, tmp_basic_hash >>= hash_function_.GetBitShiftSize(), generated_hash >>=
            hash_function_.GetBitShiftSize()) {
      if ((tmp_basic_hash & bit_mask) != (generated_hash & bit_mask)) {
        ++mismatch_count;
      }

    }
    EXPECT_EQ(1, mismatch_count);
  }
}

