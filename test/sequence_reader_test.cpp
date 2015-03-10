/*
 * sequence_reader_test.cpp
 *
 *  Created on: 2010/09/17
 *      Author: shu
 */

#include <gtest/gtest.h>
#include <stdint.h>
#include <iostream>
#include "../src/fasta_sequence_reader.h"
#include "../src/sequence.h"

using namespace std;

TEST(SequenceReaderTest, Read) {
  ifstream in("./test/test_protein.fa");
  if (!in) {
    cout << "open error" << endl;
    exit(1);
  }
  FastaSequenceReader reader(in);
  bool ret;
  string name;
  string sequence;

  ret = reader.Read(name, sequence);
  EXPECT_EQ("test0", name);
  EXPECT_EQ("SESEWWAKKPILQNFMKGAYCPLND",sequence);
  ret = reader.Read(name, sequence);
  EXPECT_EQ("test1", name);
  EXPECT_EQ("TLFVTTSNMFMIQSMATNQLGLRPL",sequence);
  ret = reader.Read(name, sequence);
  EXPECT_EQ("test2", name);
  EXPECT_EQ("MKTPQALSQLDRWEQLLESQYSFMF",sequence);
  ret = reader.Read(name, sequence);
  EXPECT_EQ("test3", name);
  EXPECT_EQ("CVDQPQIQRCNFEWYHRMQVWGRCD",sequence);
  ret = reader.Read(name, sequence);
  EXPECT_EQ("test4", name);
  EXPECT_EQ("ENLVSSSEYRADNWWYVLWIIQLFD",sequence);
  ret = reader.Read(name, sequence);
  EXPECT_EQ("test5", name);
  EXPECT_EQ("AAAAAAAAAAAAAAAAAAAAAAAAA",sequence);
  ret = reader.Read(name, sequence);
  EXPECT_EQ("test6", name);
  EXPECT_EQ("AAAAAAAAAAAAAAAMKGAYCPLND",sequence);
  ret = reader.Read(name, sequence);
  EXPECT_EQ("test7", name);
  EXPECT_EQ("ENLVSSSEYRADNYVLWIIQLFD",sequence);
  ret = reader.Read(name, sequence);
  EXPECT_FALSE(ret);

  in.close();
}
