/*
 * one_mismatch_hash_generator.cpp
 *
 *  Created on: 2013/05/03
 *      Author: shu
 */

#include "one_mismatch_hash_generator.h"

OneMismatchHashGenerator::OneMismatchHashGenerator() {

}

OneMismatchHashGenerator::~OneMismatchHashGenerator() {

}

int OneMismatchHashGenerator::GenerateOneMismatchHash(const Hash basic_hash, const uint32_t k_mer_length, const AlphabetCoder::Code max_code,
      const uint32_t bit_shift_size, std::vector<Hash> *one_mismatch_hashs) {
  uint32_t hash_mask = (1 << bit_shift_size) - 1;

  for (uint32_t i = 0; i < k_mer_length; ++i) {
    uint32_t i_bit_shift_size = (i*bit_shift_size);
    AlphabetCoder::Code basic_c = (basic_hash >> i_bit_shift_size) & hash_mask;
    Hash hash_base = basic_hash & ~(hash_mask << i_bit_shift_size);
    for (AlphabetCoder::Code c = 0; c <= max_code; ++c) {
      if (c != basic_c) {
        one_mismatch_hashs->push_back(hash_base + (c << i_bit_shift_size));
      }
    }
  }
  return 0;
}

