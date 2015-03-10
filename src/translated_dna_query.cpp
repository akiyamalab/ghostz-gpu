/*
 * translated_dna_query.cpp
 *
 *  Created on: 2011/04/12
 *      Author: shu
 */

#include "query.h"
#include "translator.h"
#include "translated_dna_query.h"
#include "sequence.h"
#include "dna_sequence.h"
#include "protein_type.h"
#include "protein_sequence.h"
#include "statistics.h"
#include "logger.h"
#include "translator.h"
#include "../ext/seg/src/seg.h"
#include <string>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <sstream>

using namespace std;

TranslatedDnaQuery::TranslatedDnaQuery(const Sequence &sequence,
		AlphabetCoder::Code sequence_delimiter, Translator &translator,
		SequenceFilterInterface *filter, ScoreMatrix score_matrix,
		Statistics::KarlinParameters &ungapped_ideal_karlin_parameters) :
		Query(), dna_length_(sequence.GetSequenceData().size()) {
	sequence_delimiter_ = sequence_delimiter;
	ProteinType t;
	AlphabetCoder coder(t);
	SetName(sequence);
	SetEncodedSequence(sequence, coder, translator, filter);
	SetStaticalParameters(score_matrix, ungapped_ideal_karlin_parameters);
}

void TranslatedDnaQuery::SetName(const Sequence &sequence) {
	name_ = sequence.GetName();
}

void TranslatedDnaQuery::SetEncodedSequence(const Sequence &sequence,
		const AlphabetCoder &coder, Translator &translator,
		SequenceFilterInterface *filter) {
	DnaSequence dna(sequence.GetName(), sequence.GetSequenceData());
	vector<ProteinSequence> proteins;
	translator.Translate(dna, proteins);
	sequence_.resize(
			TranslatedDnaQuery::GetSequenceLength(
					dna.GetSequenceData().size()));
	sequence_[0] = sequence_delimiter_;
	uint32_t offset = 1;
	uint32_t base_length = sequence.GetSequenceData().size() / 3;
	for (int i = 0; i < 6; ++i) {
		string s = proteins[i].GetSequenceData();
		string masked_s;
		if (!filter->Filter(s, &masked_s)) {
			coder.Encode(&masked_s[0], s.size(), &sequence_[offset]);
		} else {
			for (uint32_t sequence_i = 0; sequence_i < s.size(); ++sequence_i) {
				sequence_[offset + sequence_i] = coder.GetUnknownCode();
			}
		}
		offset += proteins[i].GetSequenceData().size();
		if (proteins[i].GetSequenceData().size() < base_length) {
			sequence_[offset] = sequence_delimiter_;
			++offset;
		}
		sequence_[offset] = sequence_delimiter_;
		++offset;
	}
}

void TranslatedDnaQuery::SetStaticalParameters(ScoreMatrix &score_matrix,
		Statistics::KarlinParameters &ungapped_ideal_karlin_parameters) {
	ProteinType t;
	Statistics statistics(t);
	bool ret = statistics.CalculateUngappedKarlinParameters(&sequence_[0],
			sequence_.size(), score_matrix, &ungapped_karlin_parameters_);
	if (!ret) {
		ungapped_karlin_parameters_ = ungapped_ideal_karlin_parameters;
	}
}
