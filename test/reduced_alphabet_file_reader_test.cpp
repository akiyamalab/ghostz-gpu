#include <gtest/gtest.h>
#include <string>
#include <sstream>
#include <stdint.h>
#include <fstream>
#include <vector>
#include "../src/reduced_alphabet_file_reader.h"


using namespace std;


TEST(RedicedAlphabetFileReaderTest, Read)
{
  string reduced_alphabet =
       "AR NDCQEG HIL KMFPST WYV B Z X";
  istringstream in(reduced_alphabet);
  vector<string> alphabet_sets;
  ReducedAlphabetFileReader reduced_alphabet_reader;
  reduced_alphabet_reader.Read(in, alphabet_sets);
  vector<string> ideal_alphabet_sets(0);
  ideal_alphabet_sets.push_back("AR");
  ideal_alphabet_sets.push_back("NDCQEG");
  ideal_alphabet_sets.push_back("HIL");
  ideal_alphabet_sets.push_back("KMFPST");
  ideal_alphabet_sets.push_back("WYV");
  ideal_alphabet_sets.push_back("B");
  ideal_alphabet_sets.push_back("Z");
  ideal_alphabet_sets.push_back("X");
  for (size_t i = 0; i < ideal_alphabet_sets.size(); ++i) {
    string alphabet_set = ideal_alphabet_sets[i];
    size_t j = 0;
    for (; j < alphabet_sets.size(); ++j) {
      if (alphabet_set == alphabet_sets[j]) {
        break;
      }
    }
    EXPECT_NE(alphabet_sets.size(), j);
  }
}
