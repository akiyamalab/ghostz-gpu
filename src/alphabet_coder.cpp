/*
 * coder.cpp
 *
 *  Created on: 2010/09/14
 *      Author: shu
 */

#include "alphabet_coder.h"
#include "sequence_type.h"
#include <stdint.h>
#include <vector>
#include <string.h>
#include <string>
#include <algorithm>
#include <limits.h>

using namespace std;

AlphabetCoder::AlphabetCoder(const SequenceType &type) :
		min_regular_letter_code_(0), max_regular_letter_code_(0), min_code_(0), max_code_(
				0), decode_map_(0), unknown_code_(0), unknown_letter_(0) {
	Set(type);
}

bool AlphabetCoder::Set(const SequenceType &type) {
	string letters("");

	letters += type.GetRegularLetters();
	letters += type.GetAmbiguousLetters();
	transform(letters.begin(), letters.end(), letters.begin(), ToUpper());

	min_regular_letter_code_ = 0;
	max_regular_letter_code_ =
			static_cast<Code>(type.GetRegularLetters().length() - 1);
	unknown_code_ = static_cast<Code>(letters.length());
	unknown_letter_ = tolower(type.GetUnknownLetter());
	max_code_ = unknown_code_;
	for (int i = 0; i < UCHAR_MAX; ++i) {
		code_map_[i] = unknown_code_;
	}
	Code code = 0;

	decode_map_.resize(max_code_ + 1, unknown_letter_);
	for (uint32_t i = 0; i < letters.length(); ++i, ++code) {
		code_map_[static_cast<int>(letters[i])] = code;
		code_map_[tolower(letters[i])] = code;
		decode_map_[code] = toupper(letters[i]);
	}
	return true;
}

AlphabetCoder::Code AlphabetCoder::Encode(char letter) const {
	return code_map_[static_cast<int>(letter)];
}

void AlphabetCoder::Encode(const char *sequence, unsigned int length,
		Code *encoded_sequence) const {
	for (uint32_t i = 0; i < length; ++i) {
		encoded_sequence[i] = code_map_[static_cast<int>(sequence[i])];
	}
}

char AlphabetCoder::Decode(Code code) const {
	return decode_map_[code];
}

void AlphabetCoder::Decode(const Code *encoded_sequence, unsigned int length,
		char *sequence) const {
	for (uint32_t i = 0; i < length; ++i) {
		sequence[i] = decode_map_[encoded_sequence[i]];
	}
}

bool AlphabetCoder::IsUnknown(char letter) const {
	char l = tolower(letter);
	return (code_map_[static_cast<int>(l)] == unknown_code_)
			&& (l != unknown_letter_);
}
