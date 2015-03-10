/*
 * reduced_alphabet_k_mer_hash_function.cpp
 *
 *  Created on: 2013/04/22
 *      Author: shu
 */

#include <fstream>
#include <assert.h>
#include "reduced_alphabet_k_mer_hash_function.h"

ReducedAlphabetKMerHashFunction::ReducedAlphabetKMerHashFunction() :
    can_calculate_hash_flag_(false), max_code_(0), k_mer_length_(0), shift_size_(0), max_hash_(0), reduced_code_map_(
        0) {
}

ReducedAlphabetKMerHashFunction::ReducedAlphabetKMerHashFunction(AlphabetCoder::Code max_code,
    uint32_t k_mer_length, std::vector<AlphabetCoder::Code> &reduced_code_map) :
    can_calculate_hash_flag_(false), max_code_(max_code), k_mer_length_(k_mer_length), shift_size_(
        0), max_hash_(0), reduced_code_map_(reduced_code_map) {
  bool normal_flag = true;
  if (normal_flag && SetShiftSize(max_code_)) {
    normal_flag = false;
  }
  if (normal_flag && SetMaxHash(max_code_, k_mer_length_, shift_size_)) {
    normal_flag = false;
  }
  can_calculate_hash_flag_ = normal_flag;
}

ReducedAlphabetKMerHashFunction::ReducedAlphabetKMerHashFunction(std::istream &is) :
    can_calculate_hash_flag_(false), max_code_(0), k_mer_length_(0), shift_size_(0), max_hash_(0), reduced_code_map_(
        0) {
  Load(is);
}

ReducedAlphabetKMerHashFunction::~ReducedAlphabetKMerHashFunction() {

}

int ReducedAlphabetKMerHashFunction::CalculateHash(const AlphabetCoder::Code *sequence,
    Hash *hash) const {
  Hash ret_hash = 0;
  for (uint32_t i = 0; i < k_mer_length_; ++i) {
    if (sequence[i] < reduced_code_map_.size() && reduced_code_map_[sequence[i]] <= max_code_) {
      ret_hash <<= shift_size_;
      ret_hash += reduced_code_map_[sequence[i]];
    } else {
      *hash = max_hash_ + 1;
      return 1;
    }
  }
  *hash = ret_hash;
  return 0;
}

uint32_t ReducedAlphabetKMerHashFunction::GetKMerLength() const {
  return k_mer_length_;
}
bool ReducedAlphabetKMerHashFunction::CanCalculateHash() const {
  return can_calculate_hash_flag_;
}

ReducedAlphabetKMerHashFunction::Hash ReducedAlphabetKMerHashFunction::GetMaxHash() const {
  return max_hash_;
}

uint32_t ReducedAlphabetKMerHashFunction::GetBitShiftSize() const {
  return shift_size_;
}

int ReducedAlphabetKMerHashFunction::SetShiftSize(AlphabetCoder::Code max_code) {
  shift_size_ = 0;
  for (AlphabetCoder::Code c = 1; c < max_code; c <<= 1, ++shift_size_) {
  }
  return 0;
}

int ReducedAlphabetKMerHashFunction::SetMaxHash(AlphabetCoder::Code max_code, uint32_t k_mer_length,
    uint32_t shift_size) {
  max_hash_ = 0;
  if ((shift_size * k_mer_length > sizeof(Hash) * 8)
      || ((shift_size * k_mer_length == sizeof(Hash) * 8) && (max_code == (1 << shift_size)))) {
    return 1;
  }
  for (uint32_t i = 0; i < k_mer_length; ++i) {
    max_hash_ <<= shift_size;
    max_hash_ += max_code;
  }
  return 0;
}

int ReducedAlphabetKMerHashFunction::Save(std::ostream &os) const {
  os.write((char *) &can_calculate_hash_flag_, sizeof(can_calculate_hash_flag_));
  os.write((char *) &max_code_, sizeof(max_code_));
  os.write((char *) &k_mer_length_, sizeof(k_mer_length_));
  os.write((char *) &shift_size_, sizeof(shift_size_));
  os.write((char *) &max_hash_, sizeof(max_hash_));
  size_t reduced_code_map_size = reduced_code_map_.size();
  os.write((char *) &reduced_code_map_size, sizeof(reduced_code_map_size));
  os.write((char *) &reduced_code_map_[0], sizeof(reduced_code_map_[0]) * reduced_code_map_.size());
  return 0;
}
int ReducedAlphabetKMerHashFunction::Load(std::istream &is) {
  is.read((char *) &can_calculate_hash_flag_, sizeof(can_calculate_hash_flag_));
  is.read((char *) &max_code_, sizeof(max_code_));
  is.read((char *) &k_mer_length_, sizeof(k_mer_length_));
  is.read((char *) &shift_size_, sizeof(shift_size_));
  is.read((char *) &max_hash_, sizeof(max_hash_));
  size_t reduced_code_map_size = 0;
  is.read((char *) &reduced_code_map_size, sizeof(reduced_code_map_size));
  reduced_code_map_.resize(reduced_code_map_size);
  is.read((char *) &reduced_code_map_[0], sizeof(reduced_code_map_[0]) * reduced_code_map_.size());

  return 0;
}
