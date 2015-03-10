/*
 * reduced_alphabet_k_mer_hash_function.h
 *
 *  Created on: 2013/04/22
 *      Author: shu
 */

#ifndef REDUCED_ALPHABET_K_MER_HASH_FUNCTION_H_
#define REDUCED_ALPHABET_K_MER_HASH_FUNCTION_H_

#include <stdint.h>
#include "alphabet_coder.h"

class ReducedAlphabetKMerHashFunction {
public:
  typedef uint32_t Hash;
  ReducedAlphabetKMerHashFunction();
  ReducedAlphabetKMerHashFunction(AlphabetCoder::Code max_code, uint32_t k_mer_length,
      std::vector<AlphabetCoder::Code> &reduced_code_map);
  ReducedAlphabetKMerHashFunction(std::istream &is);
  virtual ~ReducedAlphabetKMerHashFunction();

  int CalculateHash(const AlphabetCoder::Code *sequence, Hash *hash) const;

  uint32_t GetKMerLength() const;
  bool CanCalculateHash() const ;
  Hash GetMaxHash() const ;
  uint32_t GetBitShiftSize() const ;

  int Save(std::ostream &os) const;
  int Load(std::istream &is);

private:
  int SetShiftSize(AlphabetCoder::Code max_code);
  int SetMaxHash(AlphabetCoder::Code max_code, uint32_t k_mer_length, uint32_t shift_size);
  bool can_calculate_hash_flag_;
  AlphabetCoder::Code max_code_;
  uint32_t k_mer_length_;
  uint32_t shift_size_;
  Hash max_hash_;
  std::vector<AlphabetCoder::Code> reduced_code_map_;
};

#endif /* REDUCED_ALPHABET_K_MER_HASH_FUNCTION_H_ */
