/*
 * variable_hash_sequences_index.cpp
 *
 *  Created on: 2013/05/21
 *      Author: shu
 */

#include <assert.h>
#include <algorithm>
#include "alphabet_coder.h"
#include "variable_hash_sequences_index.h"

using namespace std;

VariableHashSequencesIndex::VariableHashSequencesIndex() {
}

VariableHashSequencesIndex::VariableHashSequencesIndex(std::istream &is) {
	Load(is);
}

VariableHashSequencesIndex::~VariableHashSequencesIndex() {

}

int VariableHashSequencesIndex::Build(
		std::vector<std::pair<HashKey, Value> >::iterator key_value_pairs_first,
		std::vector<std::pair<HashKey, Value> >::iterator key_value_pairs_last,
		const HashFunction &hash_function) {
	hash_function_ = hash_function;
	index_.Build(key_value_pairs_first, key_value_pairs_last);
	return 0;
}

int VariableHashSequencesIndex::GetHashKey(const AlphabetCoder::Code *sequence,
		VariableHashSequencesIndex::HashKey *hash_key,
		uint32_t *hash_length) const {
	return hash_function_.CalculateHash(sequence, hash_key, hash_length);
}

int VariableHashSequencesIndex::GetHashKey(const AlphabetCoder::Code *sequence,
		VariableHashSequencesIndex::HashKey *hash_key) const {
	return hash_function_.CalculateHash(sequence, hash_key);
}

int VariableHashSequencesIndex::Load(std::istream &is) {
	hash_function_.Load(is);
	index_.Load(is);
	return 0;
}
int VariableHashSequencesIndex::Save(std::ostream &os) {
	hash_function_.Save(os);
	index_.Save(os);
	return 0;
}
