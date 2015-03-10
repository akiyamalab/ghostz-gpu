/*
 * dense_hash_inverted_index.h
 *
 *  Created on: 2012/11/07
 *      Author: shu
 */

#ifndef DENSE_HASH_INVERTED_INDEX_H_
#define DENSE_HASH_INVERTED_INDEX_H_

#include <vector>
#include <fstream>
#include <stdint.h>
#include <algorithm>
#include <assert.h>
#include <vector>
#include <iostream>
#include <tr1/unordered_map>
#include "alphabet_coder.h"

template<typename TKey, typename TValue>
class DenseHashInvertedIndex {
public:
  typedef TKey Key;
  typedef TValue Value;
  DenseHashInvertedIndex();
  template<typename TInputIterator>
  DenseHashInvertedIndex(TInputIterator first, TInputIterator last);
  DenseHashInvertedIndex(std::istream &is);
  ~DenseHashInvertedIndex();
  template<typename TInputIterator>
  int Build(TInputIterator first, TInputIterator last);
  int GetValues(Key key, Value const ** values, size_t *length) const;
  int Load(std::istream &is);
  int Save(std::ostream &os);
private:
  std::tr1::unordered_map<Key, uint32_t> key_map_;
  std::vector<uint32_t> offsets_;
  std::vector<Value> values_;

};

template<typename TKey, typename TValue>
DenseHashInvertedIndex<TKey, TValue>::DenseHashInvertedIndex() :
    key_map_(0), offsets_(0), values_(0) {
}

template<typename TKey, typename TValue>
template<typename TInputIterator>
DenseHashInvertedIndex<TKey, TValue>::DenseHashInvertedIndex(TInputIterator first,
    TInputIterator last) :
    key_map_(0), offsets_(0), values_(0) {
  Build(first, last);
}

template<typename TKey, typename TValue>
DenseHashInvertedIndex<TKey, TValue>::DenseHashInvertedIndex(std::istream &is) :
    key_map_(0), offsets_(0), values_(0) {
  Load(is);
}

template<typename TKey, typename TValue>
DenseHashInvertedIndex<TKey, TValue>::~DenseHashInvertedIndex() {
}

template<typename TKey, typename TValue>
template<typename TInputIterator>
int DenseHashInvertedIndex<TKey, TValue>::Build(TInputIterator first, TInputIterator last) {
  key_map_.clear();
  std::vector<uint32_t> counters(0);
  uint32_t next_key_index = 0;
  typename std::tr1::unordered_map<Key, uint32_t>::iterator find_it;
  for (TInputIterator it = first; it != last; ++it) {
    find_it = key_map_.find(it->first);
    if (find_it == key_map_.end()) {
      key_map_.insert(std::make_pair(it->first, next_key_index));
      ++next_key_index;
      counters.push_back(1);
    } else {
      const uint32_t &key_index = find_it->second;
      ++counters[key_index];
    }
  }

  offsets_.resize(counters.size() + 1, 0);
  uint32_t sum = 0;
  for (size_t i = 0; i < counters.size(); ++i) {
    offsets_[i] = sum;
    sum += counters[i];
    counters[i] = 0;
  }
  offsets_[offsets_.size() - 1] = sum;

  size_t key_values_length = std::distance(first, last);
  values_.resize(key_values_length);

  for (TInputIterator it = first; it != last; ++it) {
    const uint32_t &key_index = key_map_.find(it->first)->second;
    values_[offsets_[key_index] + counters[key_index]] = it->second;
    ++counters[key_index];
  }
  return 0;
}

template<typename TKey, typename TValue>
inline int DenseHashInvertedIndex<TKey, TValue>::GetValues(Key key, Value const ** values, size_t *length) const {
  typename std::tr1::unordered_map<Key, uint32_t>::const_iterator find_it = key_map_.find(key);
  if (find_it != key_map_.end()) {
    *values = &values_[offsets_[find_it->second]];
    *length = offsets_[find_it->second + 1] - offsets_[find_it->second];
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
int DenseHashInvertedIndex<TKey, TValue>::Load(std::istream &is) {
  size_t key_map_size = 0;
  is.read((char *) &key_map_size, sizeof(key_map_size));
  std::vector<std::pair<Key, uint32_t> > key_map_array(key_map_size);
  for (typename std::vector<std::pair<Key, uint32_t> >::iterator it = key_map_array.begin();
      it != key_map_array.end(); ++it) {
    is.read((char *) &(it->first), sizeof(it->first));
    is.read((char *) &(it->second), sizeof(it->second));
  }
  key_map_.clear();
  key_map_.insert(key_map_array.begin(), key_map_array.end());
  size_t offsets_size = 0;
  is.read((char *) &offsets_size, sizeof(offsets_size));
  offsets_.resize(offsets_size);
  is.read((char *) &offsets_[0], sizeof(offsets_[0]) * offsets_.size());
  size_t values_size = 0;
  is.read((char *) &values_size, sizeof(values_size));
  values_.resize(values_size);
  is.read((char *) &values_[0], sizeof(values_[0]) * values_.size());
  return 0;
}

template<typename TKey, typename TValue>
int DenseHashInvertedIndex<TKey, TValue>::Save(std::ostream &os) {
  size_t key_map_size = key_map_.size();
  os.write((char *) &key_map_size, sizeof(key_map_size));
  for (typename std::tr1::unordered_map<Key, uint32_t>::iterator it = key_map_.begin();
      it != key_map_.end(); ++it) {
    os.write((char *) &(it->first), sizeof(it->first));
    os.write((char *) &(it->second), sizeof(it->second));
  }
  size_t offsets_size = offsets_.size();
  os.write((char *) &offsets_size, sizeof(offsets_size));
  os.write((char *) &offsets_[0], sizeof(offsets_[0]) * offsets_.size());
  size_t values_size = values_.size();
  os.write((char *) &values_size, sizeof(values_size));
  os.write((char *) &values_[0], sizeof(values_[0]) * values_.size());
  return 0;
}

#endif /* DENSE_HASH_INVERTED_INDEX_H_ */
