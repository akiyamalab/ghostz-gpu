#include <gtest/gtest.h>
#include <string>
#include <sstream>
#include <stdint.h>
#include <fstream>
#include <vector>
#include "../src/protein_type.h"
#include "../src/reduced_alphabet_file_reader.h"
#include "../src/reduced_alphabet_coder.h"
#include "../src/alphabet_coder.h"

using namespace std;

TEST(RedicedAlphabetCoderTest, Read) {
  const string alphabet_sets_string = "A KR EDNQ C G H JLVM FYW P ST";
  istringstream in(alphabet_sets_string);
  vector<string> alphabet_sets;
  ReducedAlphabetFileReader reduced_alphabet_reader;
  reduced_alphabet_reader.Read(in, alphabet_sets);
  ProteinType type;
  AlphabetCoder alphabet_coder(type);
  ReducedAlphabetCoder reduced_alphabet_coder(type, alphabet_sets);
  for (AlphabetCoder::Code code1 = alphabet_coder.GetMinCode(); code1 <= alphabet_coder.GetMaxCode();
      ++code1) {
    char c1 = toupper(alphabet_coder.Decode(code1));
    size_t set_id = 0;
    for (; set_id < alphabet_sets.size(); ++set_id) {
      size_t count = 0;
      for (; count < alphabet_sets[set_id].size(); ++count) {
        if (toupper(c1) == toupper(alphabet_sets[set_id][count])) {
          break;
        }
      }
      if (count < alphabet_sets[set_id].size()) {
        break;
      }
    }
    for (AlphabetCoder::Code code2 = alphabet_coder.GetMinCode(); code2 <= alphabet_coder.GetMaxCode();
        ++code2) {
      char c2 = toupper(alphabet_coder.Decode(code2));
      if (set_id < alphabet_sets.size()) {
        size_t count = 0;
        for (; count < alphabet_sets[set_id].size(); ++count) {
          if (c2 == alphabet_sets[set_id][count]) {
            break;
          }
        }
        if (count < alphabet_sets[set_id].size()) {
          EXPECT_EQ(reduced_alphabet_coder.Encode(c1), reduced_alphabet_coder.Encode(c2));
        } else {
          EXPECT_NE(reduced_alphabet_coder.Encode(c1), reduced_alphabet_coder.Encode(c2));
        }
      } else {
        if (c1 != c2) {
          EXPECT_NE(reduced_alphabet_coder.Encode(c1), reduced_alphabet_coder.Encode(c2));
        }
      }
    }
  }
}
