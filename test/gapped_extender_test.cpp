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
#include "../src/gapped_extender.h"

using namespace std;

TEST(GappedExtenderTest, Extender)
{
  DnaType type;
  AlphabetCoder coder( type );
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

  int result_score = 0;
  int seq0_pos;
  int seq1_pos;
  EditBlocks edit_blocks;

  GappedExtender e;
  e.ExtendOneSide(&encoded_seq0[1], encoded_seq1.size() - 1, &encoded_seq1[1], delimiter_code, false, score_matrix, 0, 1,
      INT_MAX, &result_score, &seq0_pos, &seq1_pos, &edit_blocks);

  EXPECT_EQ(7, result_score);
  EXPECT_EQ(4, seq0_pos);
  EXPECT_EQ(3, seq1_pos);

  vector<EditBlocks::EditOpType> edit_array =  edit_blocks.ToVector();

#if 0
  int seq0_position = seq0_pos;
  int seq1_position = seq1_pos;
  for (int i = (int)edit_array.size() - 1; i >= 0; --i) {
    EditBlocks::EditOpType op = edit_array[i];
    switch (op) {
    case EditBlocks::kSubstitution:
      if (seq0[seq0_position] == seq1[seq1_position]) {
        cout << seq0[seq0_position] << " - " << seq1[seq1_position];
      } else {
        cout << seq0[seq0_position] << "   " << seq1[seq1_position];
      }
      --seq0_position;
      --seq1_position;
      break;
    case EditBlocks::kGapInSeq0:
      cout << "|   " << seq1[seq1_position];
      --seq1_position;
      break;
    case EditBlocks::kGapInSeq1:
      cout << seq0[seq0_position] << "   |";
      --seq0_position;
      break;
    default:
      assert(!"invalid edit op");
      break;
    }
    cout << endl;
  }
#endif
  EXPECT_EQ(EditBlocks::kSubstitution, edit_array[0]);
  EXPECT_EQ(EditBlocks::kSubstitution, edit_array[1]);
  EXPECT_EQ(EditBlocks::kGapInSeq1, edit_array[2]);
  EXPECT_EQ(EditBlocks::kSubstitution, edit_array[3]);
  EXPECT_EQ(EditBlocks::kSubstitution, edit_array[4]);
}
