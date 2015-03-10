/*
 * translator_test.cpp
 *
 *  Created on: 2010/10/03
 *      Author: shu
 */

#include <gtest/gtest.h>
#include <string>
#include <stdint.h>
#include "../src/translator.h"
#include "../src/sequence_type.h"
#include "../src/protein_type.h"
#include "../src/dna_sequence.h"
#include "../src/protein_sequence.h"

using namespace std;

TEST(TranslatorTest, translate)
{
  DnaSequence seq0(string("test0"), string("GCCCGCCACCT"));
  ProteinType protein;
  Translator t;
  vector<ProteinSequence> proteins;

  t.Translate(seq0, proteins);
  EXPECT_EQ(string("ARH"),proteins[0].GetSequenceData());
  EXPECT_EQ(string("PAT"),proteins[1].GetSequenceData());
  EXPECT_EQ(string("PPP"),proteins[2].GetSequenceData());
  EXPECT_EQ(string("RWR"),proteins[3].GetSequenceData());
  EXPECT_EQ(string("GGG"),proteins[4].GetSequenceData());
  EXPECT_EQ(string("VAG"),proteins[5].GetSequenceData());

  DnaSequence seq1(string("test1"), string("GCCCGCCACC"));

  t.Translate(seq1, proteins);
  EXPECT_EQ(string("ARH"),proteins[0].GetSequenceData());
  EXPECT_EQ(string("PAT"),proteins[1].GetSequenceData());
  EXPECT_EQ(string("PP"),proteins[2].GetSequenceData());
  EXPECT_EQ(string("GGG"),proteins[3].GetSequenceData());
  EXPECT_EQ(string("VAG"),proteins[4].GetSequenceData());
  EXPECT_EQ(string("WR"),proteins[5].GetSequenceData());
}
