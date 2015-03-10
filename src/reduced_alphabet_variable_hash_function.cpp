/*
 * reduced_alphabet_variable_hash_function.cpp
 *
 *  Created on: May 15, 2013
 *      Author: shu
 */

#include <iostream>
#include <fstream>
#include <assert.h>
#include "reduced_alphabet_variable_hash_function.h"

ReducedAlphabetVariableHashFunction::ReducedAlphabetVariableHashFunction() :
		can_calculate_hash_flag_(false), max_code_(0), score_threshold_(0), max_length_(
				0), shift_size_(0), max_hash_(0), reduced_code_map_(0) {
}

ReducedAlphabetVariableHashFunction::ReducedAlphabetVariableHashFunction(
		AlphabetCoder::Code max_code, int score_threshold, uint32_t max_length,
		std::vector<AlphabetCoder::Code> &reduced_code_map,
		std::vector<int> code_scores) :
		can_calculate_hash_flag_(false), max_code_(max_code + 1), score_threshold_(
				score_threshold), max_length_(max_length), shift_size_(0), max_hash_(
				0), reduced_code_map_(reduced_code_map), code_scores_(
				code_scores) {
	bool normal_flag = true;
	if (normal_flag && SetShiftSize(max_code_)) {
		normal_flag = false;
	}
	if (normal_flag && SetMaxHash(max_code_, max_length_, shift_size_)) {
		normal_flag = false;
	}
	can_calculate_hash_flag_ = normal_flag;

	code_scores_.resize(code_scores_.size() + 1);
	for (int i = code_scores_.size() - 1; i > 0; --i) {
		code_scores_[i] = code_scores_[i - 1];
	}
	code_scores_[0] = 0;
	for (size_t i = 0; i < reduced_code_map_.size(); ++i) {
		++reduced_code_map_[i];
	}
}

ReducedAlphabetVariableHashFunction::ReducedAlphabetVariableHashFunction(
		std::istream &is) :
		can_calculate_hash_flag_(false), max_code_(0), max_length_(0), shift_size_(
				0), max_hash_(0), reduced_code_map_(0) {
	Load(is);
}

ReducedAlphabetVariableHashFunction::~ReducedAlphabetVariableHashFunction() {

}

int ReducedAlphabetVariableHashFunction::CalculateHash(
		const AlphabetCoder::Code *sequence, Hash *hash,
		uint32_t *hash_length) const {
	Hash ret_hash = 0;
	int score = 0;
	uint32_t ret_hash_length = 0;
	for (ret_hash_length = 0; ret_hash_length < max_length_;
			++ret_hash_length) {
		if (sequence[ret_hash_length] < reduced_code_map_.size()
				&& reduced_code_map_[sequence[ret_hash_length]] <= max_code_) {
			if (ret_hash_length >= kMinHashLength
					&& score >= score_threshold_) {
				break;
			}
			ret_hash <<= shift_size_;
			AlphabetCoder::Code c = reduced_code_map_[sequence[ret_hash_length]];
			ret_hash += c;
			score += code_scores_[c];
#if 0
			std::cout << "c : " << (int) c << " score : " << score << std::endl;
#endif
		} else {
			*hash = max_hash_ + 1;
			*hash_length = 0;
			return 1;
		}
	}
	*hash = ret_hash;
	*hash_length = ret_hash_length;
	return 0;
}

int ReducedAlphabetVariableHashFunction::CalculateHash(
		const AlphabetCoder::Code *sequence, Hash *hash) const {
	uint32_t hash_length = 0;
	return CalculateHash(sequence, hash, &hash_length);
}

uint32_t ReducedAlphabetVariableHashFunction::GetMaxLength() const {
	return max_length_;
}
bool ReducedAlphabetVariableHashFunction::CanCalculateHash() const {
	return can_calculate_hash_flag_;
}

ReducedAlphabetVariableHashFunction::Hash ReducedAlphabetVariableHashFunction::GetMaxHash() const {
	return max_hash_;
}

uint32_t ReducedAlphabetVariableHashFunction::GetBitShiftSize() const {
	return shift_size_;
}

int ReducedAlphabetVariableHashFunction::SetShiftSize(
		AlphabetCoder::Code max_code) {
	shift_size_ = 0;
	for (AlphabetCoder::Code c = 1; c < max_code; c <<= 1, ++shift_size_) {
	}
	return 0;
}

int ReducedAlphabetVariableHashFunction::SetMaxHash(
		AlphabetCoder::Code max_code, uint32_t max_length,
		uint32_t shift_size) {
	max_hash_ = 0;
	if ((shift_size * max_length > sizeof(Hash) * 8)
			|| ((shift_size * max_length == sizeof(Hash) * 8)
					&& ((max_code + 1) == (1 << shift_size)))) {
		return 1;
	}
	for (uint32_t i = 0; i < max_length; ++i) {
		max_hash_ <<= shift_size;
		max_hash_ += max_code;
	}
	return 0;
}

int ReducedAlphabetVariableHashFunction::Save(std::ostream &os) const {
	os.write((char *) &can_calculate_hash_flag_,
			sizeof(can_calculate_hash_flag_));
	os.write((char *) &max_code_, sizeof(max_code_));
	os.write((char *) &score_threshold_, sizeof(score_threshold_));
	os.write((char *) &max_length_, sizeof(max_length_));
	os.write((char *) &shift_size_, sizeof(shift_size_));
	os.write((char *) &max_hash_, sizeof(max_hash_));
	size_t reduced_code_map_size = reduced_code_map_.size();
	os.write((char *) &reduced_code_map_size, sizeof(reduced_code_map_size));
	os.write((char *) &reduced_code_map_[0],
			sizeof(reduced_code_map_[0]) * reduced_code_map_.size());
	size_t code_scores_size = code_scores_.size();
	os.write((char *) &code_scores_size, sizeof(code_scores_size));
	os.write((char *) &code_scores_[0],
			sizeof(code_scores_[0]) * code_scores_.size());
	return 0;
}
int ReducedAlphabetVariableHashFunction::Load(std::istream &is) {
	is.read((char *) &can_calculate_hash_flag_,
			sizeof(can_calculate_hash_flag_));
	is.read((char *) &max_code_, sizeof(max_code_));
	is.read((char *) &score_threshold_, sizeof(score_threshold_));
	is.read((char *) &max_length_, sizeof(max_length_));
	is.read((char *) &shift_size_, sizeof(shift_size_));
	is.read((char *) &max_hash_, sizeof(max_hash_));
	size_t reduced_code_map_size = 0;
	is.read((char *) &reduced_code_map_size, sizeof(reduced_code_map_size));
	reduced_code_map_.resize(reduced_code_map_size);
	is.read((char *) &reduced_code_map_[0],
			sizeof(reduced_code_map_[0]) * reduced_code_map_.size());
	size_t code_scores_size = 0;
	is.read((char *) &code_scores_size, sizeof(code_scores_size));
	code_scores_.resize(code_scores_size);
	is.read((char *) &code_scores_[0],
			sizeof(code_scores_[0]) * code_scores_.size());

	return 0;
}
