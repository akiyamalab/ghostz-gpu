#include <gtest/gtest.h>
#include <string>
#include <stdint.h>
#include <fstream>
#include "../src/alphabet_coder.h"
#include "../src/sequence_type.h"
#include "../src/dna_type.h"
#include "../src/perfect_hash_inverted_index.h"

using namespace std;

class PerfectHashInvertedIndexTest: public ::testing::Test {
protected:
  typedef uint32_t Value;
  typedef uint32_t Key;
  static const Key kMaxKey = 2;
  typedef PerfectHashInvertedIndex<Key, Value> TestIndex;
  virtual void SetUp() {
    key_value_pairs_.push_back(make_pair((Key)1, (Value)0));
    key_value_pairs_.push_back(make_pair((Key)2, (Value)1));
    key_value_pairs_.push_back(make_pair((Key)2, (Value)2));
  }

  virtual void TearDown() {

  }

  std::vector< std::pair<Key, Value> > key_value_pairs_;
  TestIndex inverted_index_;
};

TEST_F(PerfectHashInvertedIndexTest, Build) {
  EXPECT_EQ(0, inverted_index_.Build(key_value_pairs_.begin(), key_value_pairs_.end()));
}

TEST_F(PerfectHashInvertedIndexTest, GetValues) {
  inverted_index_.Build(key_value_pairs_.begin(), key_value_pairs_.end());
  const TestIndex::Value *values;
  size_t length = 0;
  Key hash = 0;
  EXPECT_EQ(1, inverted_index_.GetValues(hash, &values, &length));
  hash = 1;
  EXPECT_EQ(0, inverted_index_.GetValues(hash, &values, &length));
  EXPECT_EQ(1, length);
  EXPECT_EQ(0, values[0]);
  hash = 2;
  EXPECT_EQ(0, inverted_index_.GetValues(hash, &values, &length));
  EXPECT_EQ(2, length);
  EXPECT_EQ(1, values[0]);
  EXPECT_EQ(2, values[1]);
}

TEST_F(PerfectHashInvertedIndexTest, LoadAndSave) {
  string test_filename = "inverted_index_test";
  ofstream ofs(test_filename.c_str(), ios::binary);
  inverted_index_.Build(key_value_pairs_.begin(), key_value_pairs_.end());
  EXPECT_EQ(0, inverted_index_.Save(ofs));
  ofs.close();
  ifstream ifs(test_filename.c_str(), ios::binary);
  TestIndex loaded_inverted_index(ifs);
  for (TestIndex::Key hash = 0; hash <= kMaxKey; ++hash) {
    const TestIndex::Value *values;
    size_t length = 0;
    const TestIndex::Value *loaded_index_values;
    size_t loaded_index_length = 0;
    int ret = inverted_index_.GetValues(hash, &values, &length);
    EXPECT_EQ(ret,
        loaded_inverted_index.GetValues(hash, &loaded_index_values, &loaded_index_length));
    if (ret == 0) {
      EXPECT_EQ(length, loaded_index_length);
      for (size_t values_i = 0; values_i < length; ++values_i) {
        EXPECT_EQ(values[values_i], loaded_index_values[values_i]);
      }
    }
  }
}

