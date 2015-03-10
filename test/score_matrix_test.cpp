/*
 * score_matrix_test.cpp
 *
 *  Created on: 2010/09/12
 *      Author: shu
 */

#include <gtest/gtest.h>
#include <string>
#include <stdint.h>
#include <fstream>
#include "../src/score_matrix.h"
#include "../src/alphabet_coder.h"
#include "../src/sequence_type.h"
#include "../src/protein_type.h"
#include "../src/dna_type.h"
#include "../src/score_matrix_reader.h"

using namespace std;

TEST(ScoreMatrixTest, ReadingFile)
{
  DnaType type;
  ScoreMatrixReader reader;
  ifstream in("./test/DNA");
  vector<int> matrix;
  uint32_t number_letters;
  reader.Read(in, type, matrix, number_letters);
  EXPECT_EQ(5, number_letters);
  for (uint32_t i = 0; i < number_letters -1; ++i) {
    for (uint32_t j = 0; j < number_letters-1; ++j) {
      if (i == j) {
        EXPECT_EQ(1, matrix[number_letters*i + j]);
      } else {
        EXPECT_EQ(-3, matrix[number_letters*i + j]);
      }
    }
  }
}
