/*
 * variable_hash_sequences_index.h
 *
 *  Created on: 2013/05/21
 *      Author: shu
 */

#ifndef VARIABLE_HASH_SEQUENCES_INDEX_H_
#define VARIABLE_HASH_SEQUENCES_INDEX_H_

#include "seed_searcher_common.h"
#include "dense_hash_inverted_index.h"

class VariableHashSequencesIndex {
public:
	typedef SeedSearcherCommon::HashFunction HashFunction;
	typedef HashFunction::Hash HashKey;
	typedef uint32_t Value;

	VariableHashSequencesIndex();
	VariableHashSequencesIndex(std::istream &is);
	virtual ~VariableHashSequencesIndex();

	int Build(
			std::vector<std::pair<HashKey, Value> >::iterator key_value_pairs_first,
			std::vector<std::pair<HashKey, Value> >::iterator key_value_pairs_last,
			const HashFunction &hash_function);

	int GetHashKey(const AlphabetCoder::Code *sequence, HashKey *hash_key,
			uint32_t *hash_length) const;

	int GetHashKey(const AlphabetCoder::Code *sequence,
			HashKey *hash_key) const;
	int GetValues(const HashKey &hash_key, Value const **values,
			size_t *length) const;
	int GetValues(const AlphabetCoder::Code *sequence, Value const **values,
			size_t *length) const;
	int Load(std::istream &is);
	int Save(std::ostream &os);
private:
	typedef DenseHashInvertedIndex<HashFunction::Hash, Value> Index;

	HashFunction hash_function_;
	Index index_;
};

inline int VariableHashSequencesIndex::GetValues(
		const VariableHashSequencesIndex::HashKey &hash_key,
		Value const **values, size_t *length) const {
	return index_.GetValues(hash_key, values, length);
}

inline int VariableHashSequencesIndex::GetValues(
		const AlphabetCoder::Code *sequence, Value const **values,
		size_t *length) const {
	HashFunction::Hash h = 0;
	int ret = GetHashKey(sequence, &h);
	if (ret) {
		return 1;
	}
	return GetValues(h, values, length);
}

#endif /* VARIABLE_HASH_SEQUENCES_INDEX_H_ */

