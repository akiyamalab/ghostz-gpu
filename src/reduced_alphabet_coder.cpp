/*
 * reduced_alphabet_coder.cpp
 *
 *  Created on: 2012/11/29
 *      Author: shu
 */

#include <iostream>
#include <stdint.h>
#include <vector>
#include <string.h>
#include <string>
#include <algorithm>
#include <limits.h>
#include "sequence_type.h"
#include "alphabet_coder.h"
#include "reduced_alphabet_coder.h"

using namespace std;

ReducedAlphabetCoder::ReducedAlphabetCoder()
{}

ReducedAlphabetCoder::ReducedAlphabetCoder(const SequenceType &type) :
    AlphabetCoder(type) {
}

ReducedAlphabetCoder::ReducedAlphabetCoder(const SequenceType &type, const std::vector<std::string> &alphabet_sets) {
  Set(type, alphabet_sets);
}

ReducedAlphabetCoder::~ReducedAlphabetCoder()
{}

bool ReducedAlphabetCoder::Set(const SequenceType &type, const std::vector<std::string> &alphabet_sets) {
  vector<string> upper_alphabet_sets(alphabet_sets.size());
  for (size_t i = 0; i < upper_alphabet_sets.size(); ++i) {
    upper_alphabet_sets[i].resize(alphabet_sets[i].size());
    transform(alphabet_sets[i].begin(), alphabet_sets[i].end(), upper_alphabet_sets[i].begin(), ToUpper());
  }

  min_regular_letter_code_ = 0;
  uint32_t number_codes = 0;

  string letters("");
  letters += type.GetRegularLetters();
  size_t ambiguous_letters_offset = letters.size();
  letters += type.GetAmbiguousLetters();
  size_t unkown_letter_offset = letters.size();
  letters += type.GetUnknownLetter();
  unknown_letter_ = toupper(type.GetUnknownLetter());

  transform(letters.begin(), letters.end(), letters.begin(), ToUpper());
  for (size_t i = 0; i < ambiguous_letters_offset; ++i) {
    char representative_alphabet = GetRepresentativeAlphabet(upper_alphabet_sets, letters[i]);
    if(representative_alphabet == letters[i]) {
      ++number_codes;
    }
  }
  max_regular_letter_code_ = number_codes - 1;
  for (size_t i = ambiguous_letters_offset; i < unkown_letter_offset; ++i) {
    char representative_alphabet = GetRepresentativeAlphabet(upper_alphabet_sets, letters[i]);
    if(representative_alphabet == letters[i]) {
      ++number_codes;
    }
  }
  char representative_alphabet = GetRepresentativeAlphabet(upper_alphabet_sets, letters[unkown_letter_offset]);
  unknown_code_ = representative_alphabet;
  if (representative_alphabet == letters[unkown_letter_offset]) {
    ++number_codes;
  }

  max_code_ = number_codes - 1;
  for (int i = 0; i < UCHAR_MAX; ++i) {
    code_map_[i] = unknown_code_;
  }

  Code code = 0;
  decode_map_.resize(max_code_ + 1,unknown_letter_);
  for (uint32_t i = 0; i < letters.length(); ++i) {
    representative_alphabet = GetRepresentativeAlphabet(upper_alphabet_sets, letters[i]);
    if (representative_alphabet == letters[i]) {
      code_map_[static_cast<int>(letters[i])] = code;
      code_map_[tolower(letters[i])] = code;
      decode_map_[code] = toupper(letters[i]);
      ++code;
    }
  }
  for (uint32_t i = 0; i < letters.length(); ++i) {
    representative_alphabet = GetRepresentativeAlphabet(upper_alphabet_sets, letters[i]);
    if (representative_alphabet != letters[i]) {
      code = code_map_[static_cast<int>(representative_alphabet)];
      code_map_[static_cast<int>(letters[i])] = code;
      code_map_[tolower(letters[i])] = code;
    }
  }
  return true;
}

char ReducedAlphabetCoder::GetRepresentativeAlphabet(const std::vector<std::string> &alphabet_sets, const char c) {
  for (size_t i = 0; i < alphabet_sets.size(); ++i) {
    size_t count = 0;
    for (; count < alphabet_sets[i].size(); ++count) {
      if (alphabet_sets[i][count] == c) {
        return alphabet_sets[i][0];
      }
    }
  }
  return c;
}

