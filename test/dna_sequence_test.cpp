/*
 * dna_seuqence_test.cpp
 *
 *  Created on: 2010/09/12
 *      Author: shu
 */

#include <gtest/gtest.h>
#include <string>
#include <stdint.h>
#include "../src/sequence.h"
#include "../src/dna_sequence.h"

using namespace std;

TEST(DnaSequneceTest, GetComplementaryStrand)
{
  string seq("GCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTATGAAAGGTGCTTATTGTCCTCTGAATGAT");
  DnaSequence dna(string("test"),seq);

  string comp_seq("ATCATTCAGAGGACAATAAGCACCTTTCATAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGC");
  EXPECT_EQ(comp_seq, dna.GetComplementaryStrand());
}
