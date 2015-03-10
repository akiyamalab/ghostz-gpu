/*
 * k_mer_sequences_index.h
 *
 *  Created on: 2013/04/23
 *      Author: shu
 */

#ifndef K_MER_SEQUENCES_INDEX_H_
#define K_MER_SEQUENCES_INDEX_H_

#include <vector>
#include <fstream>
#include "reduced_alphabet_k_mer_hash_function.h"
#include "dense_hash_inverted_index.h"

class KMerSequencesIndex {
public:
	typedef ReducedAlphabetKMerHashFunction HashFunction;
	typedef HashFunction::Hash HashKey;
	typedef uint32_t Value;

	KMerSequencesIndex();
	KMerSequencesIndex(std::istream &is);
	virtual ~KMerSequencesIndex();
	int Build(
			std::vector<std::pair<HashKey, Value> >::iterator key_value_pairs_first,
			std::vector<std::pair<HashKey, Value> >::iterator key_value_pairs_last,
			const HashFunction &hash_function);
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

inline int KMerSequencesIndex::GetValues(
		const KMerSequencesIndex::HashKey &hash_key, Value const **values,
		size_t *length) const {
	return index_.GetValues(hash_key, values, length);
}

inline int KMerSequencesIndex::GetValues(const AlphabetCoder::Code *sequence,
		Value const **values, size_t *length) const {
	HashFunction::Hash h = 0;
	int ret = GetHashKey(sequence, &h);
	if (ret) {
		return 1;
	}
	return GetValues(h, values, length);
}

#endif /* K_MER_SEQUENCES_INDEX_H_ */
