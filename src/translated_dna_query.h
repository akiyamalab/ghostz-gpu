/*
 * translated_dna_query.h
 *
 *  Created on: 2011/04/12
 *      Author: shu
 */

#ifndef TRANSLATED_DNA_QUERY_H_
#define TRANSLATED_DNA_QUERY_H_

#include "query.h"
#include "translator.h"
#include "sequence.h"
#include "protein_sequence.h"
#include "translator.h"
#include "sequence_filter_interface.h"
#include <stdint.h>
#include <string>
#include <vector>
#include <iostream>

class TranslatedDnaQuery: public Query {
public:
	static uint32_t GetSequenceLength(uint32_t length) {
		if (length < 3) {
			return 0;
		} else {
			uint32_t t = length / 3;
			return 7 + t * 6;
		}
	}

	TranslatedDnaQuery(const Sequence &sequence,
			AlphabetCoder::Code sequence_delimiter, Translator &translator,
			SequenceFilterInterface *filter, ScoreMatrix score_matrix,
			Statistics::KarlinParameters &ungapped_ideal_karlin_parameters);

	// return the max length of translated sequence
	uint32_t GetRealSequenceLength() {
		return dna_length_ / 3;
	}

	uint32_t GetRealStart(uint32_t pos) {
		uint32_t base_length = dna_length_ / 3 + 1;
		uint32_t f = (pos - 1) / base_length;
		if (f < 3) {
			return ((pos % base_length) - 1) * 3 + f + 1;
		} else {
			return dna_length_ - (((pos % base_length) - 1) * 3 + f - 2) + 1;
		}
	}

	uint32_t GetRealEnd(uint32_t pos) {
		uint32_t base_length = dna_length_ / 3 + 1;
		uint32_t f = (pos - 1) / base_length;
		if (f < 3) {
			return GetRealStart(pos) + 2;
		} else {
			return GetRealStart(pos) - 2;
		}
	}

	uint32_t GetDistanceFromStartDelimiter(uint32_t pos) {
		uint32_t base_length = dna_length_ / 3 + 1;
		uint32_t f = (pos - 1) / base_length;
		return pos - f * base_length;
	}

	uint32_t GetDistanceFromEndDelimiter(uint32_t pos) {
		uint32_t base_length = dna_length_ / 3 + 1;
		uint32_t f = (pos - 1) / base_length;
		return (f + 1) * base_length - pos;
	}

protected:
	void SetName(const Sequence &sequence);
	void SetEncodedSequence(const Sequence &sequence,
			const AlphabetCoder &coder, Translator &translator,
			SequenceFilterInterface *filter);
	void SetStaticalParameters(ScoreMatrix &score_matrix,
			Statistics::KarlinParameters &ungapped_ideal_karlin_parameters);

private:
	uint32_t dna_length_;
};

#endif /* TRANSLATED_DNA_QUERY_H_ */
