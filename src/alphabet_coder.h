/*
 * alphabet_coder.h
 *
 *  Created on: 2010/09/12
 *      Author: shu
 */

#ifndef ALPHABET_CODER_H_
#define ALPHABET_CODER_H_

#include <string>
#include <vector>
#include <limits.h>

class SequenceType;
class Sequence;

class AlphabetCoder {
public:
	typedef unsigned char Code;

	AlphabetCoder() :
			min_regular_letter_code_(0), max_regular_letter_code_(0), min_code_(
					0), max_code_(0), decode_map_(0), unknown_code_(0), unknown_letter_(
					0) {
	}
	AlphabetCoder(const SequenceType &type);

	virtual ~AlphabetCoder() {
	}

	bool Set(const SequenceType &type);

	Code Encode(char letter) const;

	void Encode(const char *sequence_data, unsigned int length,
			Code *encoded_sequence) const;

	char Decode(Code code) const;

	void Decode(const Code *encoded_sequence, unsigned int length,
			char *sequence_data) const;

	bool IsUnknown(char letter) const;

	Code GetMinRegularLetterCode() const {
		return min_regular_letter_code_;
	}

	Code GetMaxRegularLetterCode() const {
		return max_regular_letter_code_;
	}

	Code GetMinCode() const {
		return min_code_;
	}

	Code GetMaxCode() const {
		return max_code_;
	}

	Code GetUnknownCode() const {
		return unknown_code_;
	}

protected:
	struct ToUpper {
		char operator()(char c) {
			return toupper(c);
		}
	};
	Code min_regular_letter_code_;
	Code max_regular_letter_code_;
	Code min_code_;
	Code max_code_;
	Code code_map_[UCHAR_MAX];
	std::vector<char> decode_map_;
	Code unknown_code_;
	char unknown_letter_;

};

#endif /* ALPHABET_CODER_H_ */
