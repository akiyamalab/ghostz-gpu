/*
 * normal_query.cpp
 *
 *  Created on: 2011/04/13
 *      Author: shu
 */

#include "alphabet_coder.h"
#include "sequence.h"
#include "statistics.h"
#include "query.h"
#include "protein_type.h"
#include "protein_query.h"
#include "sequence_filter_interface.h"
#include <string>
#include <stdexcept>
#include <vector>
#include <iostream>

using namespace std;

ProteinQuery::ProteinQuery(const Sequence &sequence,
		AlphabetCoder::Code sequence_delimiter, SequenceFilterInterface *filter,
		ScoreMatrix score_matrix,
		Statistics::KarlinParameters &ungapped_ideal_karlin_parameters) :
		Query() {
	sequence_delimiter_ = sequence_delimiter;
	ProteinType t;
	AlphabetCoder coder(t);
	SetName(sequence);
	SetEncodedSequence(sequence, coder, filter);
	SetStaticalParameters(score_matrix, ungapped_ideal_karlin_parameters);
}

void ProteinQuery::SetName(const Sequence &sequence) {
	name_ = sequence.GetName();
}

void ProteinQuery::SetEncodedSequence(const Sequence &sequence,
		const AlphabetCoder &coder,  SequenceFilterInterface *filter) {
	sequence_.resize(sequence.GetSequenceData().length() + 2);
	sequence_[0] = sequence_delimiter_;
	string s = sequence.GetSequenceData();
	string masked_s;
	if (!filter->Filter(s, &masked_s)) {
		coder.Encode(&masked_s[0], masked_s.length(), &sequence_[1]);
		sequence_[sequence_.size() - 1] = sequence_delimiter_;
	}
}

void ProteinQuery::SetStaticalParameters(ScoreMatrix &score_matrix,
		Statistics::KarlinParameters &ungapped_ideal_karlin_parameters) {
	ProteinType t;
	Statistics statistics(t);
	bool ret = statistics.CalculateUngappedKarlinParameters(&sequence_[0],
			sequence_.size(), score_matrix, &ungapped_karlin_parameters_);
	if (!ret) {
		ungapped_karlin_parameters_ = ungapped_ideal_karlin_parameters;
	}
}

