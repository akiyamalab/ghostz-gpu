/*
 * reduced_alphabet_variable_hash_function.h
 *
 *  Created on: May 15, 2013
 *      Author: shu
 */

#ifndef REDUCED_ALPHABET_VARIABLE_HASH_FUNCTION_H_
#define REDUCED_ALPHABET_VARIABLE_HASH_FUNCTION_H_

#include <stdint.h>
#include "alphabet_coder.h"
#include "score_matrix.h"

class ReducedAlphabetVariableHashFunction {
public:
  typedef uint32_t Hash;
  static const uint32_t kMinHashLength = 6;
  ReducedAlphabetVariableHashFunction();
  ReducedAlphabetVariableHashFunction(AlphabetCoder::Code max_code, int score_threshold,
      uint32_t max_length, std::vector<AlphabetCoder::Code> &reduced_code_map,
      std::vector<int> code_scores);
  ReducedAlphabetVariableHashFunction(std::istream &is);
  virtual ~ReducedAlphabetVariableHashFunction();

  int CalculateHash(const AlphabetCoder::Code *sequence, Hash *hash, uint32_t *hash_length) const;
  int CalculateHash(const AlphabetCoder::Code *sequence, Hash *hash) const;

  uint32_t GetMaxLength() const;
  bool CanCalculateHash() const;
  Hash GetMaxHash() const;
  uint32_t GetBitShiftSize() const;

  int Save(std::ostream &os) const;
  int Load(std::istream &is);

private:
  int SetShiftSize(AlphabetCoder::Code max_code);
  int SetMaxHash(AlphabetCoder::Code max_code, uint32_t k_mer_length, uint32_t shift_size);
  bool can_calculate_hash_flag_;
  AlphabetCoder::Code max_code_;
  int score_threshold_;
  uint32_t max_length_;
  uint32_t shift_size_;
  Hash max_hash_;
  std::vector<AlphabetCoder::Code> reduced_code_map_;
  std::vector<int> code_scores_;
};

#endif /* REDUCED_ALPHABET_VARIABLE_HASH_FUNCTION_H_ */
