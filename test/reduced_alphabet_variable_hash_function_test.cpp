/*
 * reduced_alphabet_variable_hash_function_test.cpp
 *
 *  Created on: May 15, 2013
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

using namespace std;

class ReducedAlphabetVariableHashFunctionTest: public ::testing::Test {
protected:
  virtual void SetUp() {
    ProteinType protein_type;
    const std::string reduced_alphabet = "A KR EDNQ C G H ILVM FYW P ST";
    const string default_protein_matrix_name = "BLOSUM62";
    const string default_protein_matrix =
        "   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *\nA  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4\n R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4\n N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4\n D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4\n C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4\n Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4\n E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4\n G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4\n H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4\n I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4\n L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4\n K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4\n M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4\n F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4\n P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4\n S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4\n T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4\n W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4\n Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4\n V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4\n B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4\n Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4\n X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4\n * -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1\n";

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
    max_length_ = 2;
    ScoreMatrixReader score_matrix_reader;
    vector<int> matrix;
    unsigned int number_letters;
    istringstream default_protein_score_matrix_is(default_protein_matrix);
    score_matrix_reader.Read(default_protein_score_matrix_is, protein_type, matrix, number_letters);
    ScoreMatrix score_matrix = ScoreMatrix(default_protein_matrix_name, &matrix[0], number_letters);
    score_threshold_ = 14;
    code_score_.resize(max_code_ + 1);
    for (uint32_t i = 0; i < number_letters; ++i) {
      if (i < reduced_code_map_.size()) {
        AlphabetCoder::Code c = reduced_code_map_[i];
        if (c <= max_code_) {
          code_score_[c] = std::max(code_score_[c],
              score_matrix.GetMatrix()[i * number_letters + i]);
        }
      }
    }
    hash_function_ = ReducedAlphabetVariableHashFunction(max_code_, score_threshold_, max_length_,
        reduced_code_map_, code_score_);
  }

  virtual void TearDown() {

  }
  AlphabetCoder::Code max_code_;
  uint32_t max_length_;
  int score_threshold_;
  std::vector<AlphabetCoder::Code> reduced_code_map_;
  std::vector<int> code_score_;
  AlphabetCoder coder_;
  ReducedAlphabetVariableHashFunction hash_function_;
};

TEST_F(ReducedAlphabetVariableHashFunctionTest, VoidConstructor) {
  ReducedAlphabetVariableHashFunction hash_function;
  EXPECT_FALSE(hash_function.CanCalculateHash());
}

TEST_F(ReducedAlphabetVariableHashFunctionTest, GetMaxLength) {
  EXPECT_EQ(max_length_, hash_function_.GetMaxLength());
}

TEST_F(ReducedAlphabetVariableHashFunctionTest, GetMaxHash) {
  uint32_t shift_size = 0;
  for (AlphabetCoder::Code c = 1; c < max_code_; c <<= 1, ++shift_size) {
  }

  ReducedAlphabetVariableHashFunction::Hash max_hash = 0;
  for (uint32_t i = 0; i < max_length_; ++i) {
    max_hash <<= shift_size;
    max_hash += max_code_ + 1;
  }
  EXPECT_EQ(max_hash, hash_function_.GetMaxHash());
}

TEST_F(ReducedAlphabetVariableHashFunctionTest, CalculateHash) {
  uint32_t shift_size = 0;
  for (AlphabetCoder::Code c = 1; c < max_code_; c <<= 1, ++shift_size) {
  }
  string k_mer_string = "GMEHQS";
  vector<AlphabetCoder::Code> k_mer_coded_string(k_mer_string.length() + 1);
  coder_.Encode(&k_mer_string[0], k_mer_string.length(), &k_mer_coded_string[0]);
  k_mer_coded_string[k_mer_string.length()] = reduced_code_map_.size();

  ReducedAlphabetVariableHashFunction::Hash expected_hash = 0;
  int score = 0;
  for (uint32_t i = 0;
      (i < max_length_) && (score < score_threshold_)
          && k_mer_coded_string[i] < reduced_code_map_.size(); ++i) {
    expected_hash <<= shift_size;
    expected_hash += reduced_code_map_[k_mer_coded_string[i]] + 1;
    score += code_score_[reduced_code_map_[k_mer_coded_string[i]]];
  }

  ReducedAlphabetVariableHashFunction::Hash hash = 0;
  EXPECT_EQ(0, hash_function_.CalculateHash(&k_mer_coded_string[0], &hash));
  EXPECT_EQ(expected_hash, hash);
  cout << hash << endl;
}

TEST_F(ReducedAlphabetVariableHashFunctionTest, SaveAndLoad) {
  string filename = "hash_function_test";
  ofstream ofs(filename.c_str());
  hash_function_.Save(ofs);
  ofs.close();
  ifstream ifs(filename.c_str());
  ReducedAlphabetVariableHashFunction loaded_hash_function(ifs);
  string k_mer_string = "GMEHQS";
  vector<AlphabetCoder::Code> k_mer_coded_string(k_mer_string.length());
  coder_.Encode(&k_mer_string[0], k_mer_string.length(), &k_mer_coded_string[0]);

  ReducedAlphabetVariableHashFunction::Hash expected_hash = 0;
  int ret = hash_function_.CalculateHash(&k_mer_coded_string[0], &expected_hash);

  ReducedAlphabetVariableHashFunction::Hash hash = 0;
  EXPECT_EQ(ret, loaded_hash_function.CalculateHash(&k_mer_coded_string[0], &hash));
  EXPECT_EQ(expected_hash, hash);
}

