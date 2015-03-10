/*
 * ungapped_extender_test.cpp
 *
 *  Created on: 2012/11/19
 *      Author: shu
 *
 */

#include "../src/dna_type.h"
#include "../src/alphabet_coder.h"
#include "../src/score_matrix.h"
#include "../src/edit_blocks.h"
#include "../src/ungapped_extender.h"
#include <gtest/gtest.h>

using namespace std;

TEST(UngappedExtendTest, ExtendOneSide) {
  DnaType type;
  AlphabetCoder coder(type);
  AlphabetCoder::Code delimiter_code = coder.GetMaxCode() + 1;
  ScoreMatrix score_matrix("test_matrix", 4, 2, -1);

  string seq0("AAAAAAC");
  vector<AlphabetCoder::Code> encoded_seq0(seq0.size() + 2);
  encoded_seq0[0] = delimiter_code;
  coder.Encode(&seq0[0], seq0.size(), &encoded_seq0[1]);
  encoded_seq0[encoded_seq0.size() - 1] = delimiter_code;

  string seq1("AACAAAG");
  vector<AlphabetCoder::Code> encoded_seq1(seq1.size() + 2);
  encoded_seq1[0] = delimiter_code;
  coder.Encode(&seq1[0], seq1.size(), &encoded_seq1[1]);
  encoded_seq1[encoded_seq1.size() - 1] = delimiter_code;
  EditBlocks edit_blocks;
  UngappedExtender e;
  int best_score;
  int best_sequence0_position;
  int best_sequence1_position;
  bool ret = e.ExtendOneSide(&encoded_seq0[1], &encoded_seq1[1], delimiter_code,
      false, score_matrix, INT_MAX, &best_score, &best_sequence0_position,
      &best_sequence1_position, &edit_blocks);
  EXPECT_EQ(true, ret);
  EXPECT_EQ(9, best_score);
  EXPECT_EQ(5, best_sequence0_position);
  EXPECT_EQ(5, best_sequence1_position);

}

TEST(UngappedExtendTest, ExtendOneSideReverse) {
  DnaType type;
  AlphabetCoder coder(type);
  AlphabetCoder::Code delimiter_code = coder.GetMaxCode() + 1;
  ScoreMatrix score_matrix("test_matrix", 4, 2, -1);

  string seq0("AAAAAAC");
  vector<AlphabetCoder::Code> encoded_seq0(seq0.size() + 2);
  encoded_seq0[0] = delimiter_code;
  coder.Encode(&seq0[0], seq0.size(), &encoded_seq0[1]);
  encoded_seq0[encoded_seq0.size() - 1] = delimiter_code;

  string seq1("AACAAAG");
  vector<AlphabetCoder::Code> encoded_seq1(seq1.size() + 2);
  encoded_seq1[0] = delimiter_code;
  coder.Encode(&seq1[0], seq1.size(), &encoded_seq1[1]);
  encoded_seq1[encoded_seq1.size() - 1] = delimiter_code;
  EditBlocks edit_blocks;
  UngappedExtender e;
  int best_score;
  int best_sequence0_position;
  int best_sequence1_position;
  bool ret = e.ExtendOneSide(&encoded_seq0[1], &encoded_seq1[1], delimiter_code,
      false, score_matrix, INT_MAX, &best_score, &best_sequence0_position,
      &best_sequence1_position, &edit_blocks);
  int reverse_best_score;
  ret = e.ExtendOneSide(&encoded_seq0[encoded_seq0.size() - 3], &encoded_seq1[encoded_seq1.size() - 3], delimiter_code,
      true, score_matrix, INT_MAX, &reverse_best_score, &best_sequence0_position,
      &best_sequence1_position, &edit_blocks);
  EXPECT_EQ(true, ret);
  EXPECT_EQ(best_score, reverse_best_score);
  EXPECT_EQ(-5, best_sequence0_position);
  EXPECT_EQ(-5, best_sequence1_position);

}

TEST(UngappedExtendTest, ExtendOneSideCutoffTest) {
  DnaType type;
  AlphabetCoder coder(type);
  AlphabetCoder::Code delimiter_code = coder.GetMaxCode() + 1;
  ScoreMatrix score_matrix("test_matrix", 4, 2, -1);

  string seq0("AAAAAAC");
  vector<AlphabetCoder::Code> encoded_seq0(seq0.size() + 2);
  encoded_seq0[0] = delimiter_code;
  coder.Encode(&seq0[0], seq0.size(), &encoded_seq0[1]);
  encoded_seq0[encoded_seq0.size() - 1] = delimiter_code;

  string seq1("AACAAAG");
  vector<AlphabetCoder::Code> encoded_seq1(seq1.size() + 2);
  encoded_seq1[0] = delimiter_code;
  coder.Encode(&seq1[0], seq1.size(), &encoded_seq1[1]);
  encoded_seq1[encoded_seq1.size() - 1] = delimiter_code;
  EditBlocks edit_blocks;
  UngappedExtender e;
  int best_score;
  int best_sequence0_position;
  int best_sequence1_position;
  bool ret = e.ExtendOneSide(&encoded_seq0[1], &encoded_seq1[1], delimiter_code,
      false, score_matrix, 1, &best_score, &best_sequence0_position,
      &best_sequence1_position, &edit_blocks);
  EXPECT_EQ(true, ret);
  EXPECT_EQ(4, best_score);
  EXPECT_EQ(1, best_sequence0_position);
  EXPECT_EQ(1, best_sequence1_position);

}

