/*
 * perfect_hash_inverted_index.h
 *
 *  Created on: 2012/11/07
 *      Author: shu
 */

#ifndef PERFECT_HASH_INVERTED_INDEX_H_
#define PERFECT_HASH_INVERTED_INDEX_H_

#include <vector>
#include <fstream>
#include <stdint.h>
#include <algorithm>
#include <assert.h>
#include "alphabet_coder.h"

template<typename TKey, typename TValue>
class PerfectHashInvertedIndex {
public:
	typedef TKey Key;
	typedef TValue Value;
	PerfectHashInvertedIndex();
	template<typename TInputIterator>
	PerfectHashInvertedIndex(TInputIterator first, TInputIterator last);
	PerfectHashInvertedIndex(std::istream &is);
	~PerfectHashInvertedIndex();
	template<typename TInputIterator>
	int Build(TInputIterator first, TInputIterator last);
	int GetValues(Key key, Value const **values, size_t *length) const;
	int Load(std::istream &is);
	int Save(std::ostream &os);
private:
	Key max_key_;
	std::vector<uint32_t> offsets_;
	std::vector<Value> values_;

};

template<typename TKey, typename TValue>
PerfectHashInvertedIndex<TKey, TValue>::PerfectHashInvertedIndex() :
		max_key_(0), offsets_(0), values_(0) {
}

template<typename TKey, typename TValue>
template<typename TInputIterator>
PerfectHashInvertedIndex<TKey, TValue>::PerfectHashInvertedIndex(TInputIterator first,
		TInputIterator last) :
		max_key_(0), offsets_(0), values_(0) {
	Build(first, last);
}

template<typename TKey, typename TValue>
PerfectHashInvertedIndex<TKey, TValue>::PerfectHashInvertedIndex(std::istream &is) :
		max_key_(0), offsets_(0), values_(0) {
	Load(is);
}

template<typename TKey, typename TValue>
PerfectHashInvertedIndex<TKey, TValue>::~PerfectHashInvertedIndex() {
}

template<typename TKey, typename TValue>
template<typename TInputIterator>
int PerfectHashInvertedIndex<TKey, TValue>::Build(TInputIterator first,
		TInputIterator last) {
	max_key_ = 0;
	for (TInputIterator it = first; it != last; ++it) {
		max_key_ = std::max(max_key_, it->first);
	}
	offsets_.clear();
	offsets_.resize(max_key_ + 2, 0);
	std::vector<uint32_t> counters(max_key_ + 1);
	for (TInputIterator it = first; it != last; ++it) {
		++counters[it->first];
	}
	uint32_t sum = 0;
	for (size_t hash = 0; hash <= max_key_; ++hash) {
		offsets_[hash] = sum;
		sum += counters[hash];
		counters[hash] = 0;
	}
	offsets_[max_key_ + 1] = sum;
	size_t key_values_length = std::distance(first, last);
	values_.resize(key_values_length);
	for (TInputIterator it = first; it != last; ++it) {
		const Key &key = it->first;
		values_[offsets_[key] + counters[key]] = it->second;
		++counters[key];
	}
	return 0;
}

template<typename TKey, typename TValue>
int PerfectHashInvertedIndex<TKey, TValue>::GetValues(Key key,
		Value const **values, size_t *length) const {
	if (key <= max_key_) {
		*values = &values_[offsets_[key]];
		*length = offsets_[key + 1] - offsets_[key];
#if 0
		assert((*values + *length) <= &values_[0] + values_.size());
#endif
		return (*length == 0) ? 1 : 0;
	} else {
		*values = NULL;
		*length = 0;
		return 1;
	}
}

template<typename TKey, typename TValue>
int PerfectHashInvertedIndex<TKey, TValue>::Load(std::istream &is) {
	is.read((char *) &max_key_, sizeof(max_key_));
	size_t offsets_size = 0;
	is.read((char *) &offsets_size, sizeof(offsets_size));
	offsets_.resize(offsets_size);
	is.read((char *) &offsets_[0], sizeof(offsets_[0]) * offsets_size);
	size_t positions_size = 0;
	is.read((char *) &positions_size, sizeof(positions_size));
	values_.resize(positions_size);
	is.read((char *) &values_[0], sizeof(values_[0]) * positions_size);
	return 0;
}

template<typename TKey, typename TValue>
int PerfectHashInvertedIndex<TKey, TValue>::Save(std::ostream &os) {
	os.write((char *) &max_key_, sizeof(max_key_));
	size_t offsets_size = offsets_.size();
	os.write((char *) &offsets_size, sizeof(offsets_size));
	os.write((char *) &offsets_[0], sizeof(offsets_[0]) * offsets_size);
	size_t positions_size = values_.size();
	os.write((char *) &positions_size, sizeof(positions_size));
	os.write((char *) &values_[0], sizeof(values_[0]) * positions_size);
	return 0;
}

#endif /* PERFECT_HASH_INVERTED_INDEX_H_ */
