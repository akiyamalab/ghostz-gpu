/*
 * one_mismatch_hash_generator.h
 *
 *  Created on: 2013/05/03
 *      Author: shu
 */

#ifndef ONE_MISMATCH_HASH_GENERATOR_H_
#define ONE_MISMATCH_HASH_GENERATOR_H_

#include <stdint.h>
#include "reduced_alphabet_k_mer_hash_function.h"
#include "alphabet_coder.h"

class OneMismatchHashGenerator {
public:
  typedef ReducedAlphabetKMerHashFunction::Hash Hash;
  OneMismatchHashGenerator();
  virtual ~OneMismatchHashGenerator();

  int GenerateOneMismatchHash(const Hash basic_hash, const uint32_t k_mer_length, const AlphabetCoder::Code max_code,
      const uint32_t bit_shift_size, std::vector<Hash> *one_mismatch_hashs);
};

#endif /* ONE_MISMATCH_HASH_GENERATOR_H_ */
