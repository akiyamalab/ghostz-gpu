/*
 * k_mer_sequences_index.cpp
 *
 *  Created on: 2013/04/23
 *      Author: shu
 */

#include "k_mer_sequences_index.h"
#include <assert.h>
#include <algorithm>

using namespace std;

KMerSequencesIndex::KMerSequencesIndex() {
}


KMerSequencesIndex::KMerSequencesIndex(std::istream &is) {
  Load(is);
}

KMerSequencesIndex::~KMerSequencesIndex() {

}

int KMerSequencesIndex::Build(std::vector< std::pair<HashKey, Value> >::iterator key_value_pairs_first,
    std::vector<std::pair<HashKey, Value> >::iterator key_value_pairs_last, const HashFunction &hash_function) {
  hash_function_ = hash_function;
  index_.Build(key_value_pairs_first, key_value_pairs_last);
  return 0;
}

int KMerSequencesIndex::GetHashKey(const AlphabetCoder::Code *sequence,
    KMerSequencesIndex::HashKey *hash_key) const {
  return hash_function_.CalculateHash(sequence, hash_key);
}

int KMerSequencesIndex::Load(std::istream &is) {
  hash_function_.Load(is);
  index_.Load(is);
  return 0;
}
int KMerSequencesIndex::Save(std::ostream &os) {
  hash_function_.Save(os);
  index_.Save(os);
  return 0;
}

