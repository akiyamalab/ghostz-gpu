/*
 * chain_filter.cpp
 *
 *  Created on: 2011/01/27
 *      Author: shu
 */

#include <vector>
#include <math.h>
#include <algorithm>
#include "seed_searcher.h"
#include "statistics.h"
#include "chain_filter.h"
#include "query.h"

using namespace std;

ChainFilter::ChainFilter(AlphabetCoder::Code sequence_delimiter,
		ScoreMatrix &score_matrix) :
		sequence_delimiter_(sequence_delimiter), score_matrix_(score_matrix), chains_(
				0) {
}

void ChainFilter::Filter(const AlphabetCoder::Code *query_sequence,
		uint32_t query_sequence_length, int cutoff,
		const AlphabetCoder::Code *database_sequence, vector<Hit> &hits) {
	Chain init_chain;
	init_chain.diagonal = 0;
	init_chain.query_end = 0;
	chains_.resize(query_sequence_length);
	for (vector<Chain>::iterator it = chains_.begin(); it < chains_.end();
			++it) {
		*it = init_chain;
	}
	sort(hits.begin(), hits.end(), DbStartComp());
	uint32_t n = 0;
	for (vector<Hit>::iterator hits_it = hits.begin(); hits_it != hits.end();
			++hits_it) {
		uint32_t diagonal = query_sequence_length - 1
				+ hits_it->database_sequence_position
				- hits_it->query_sequence_position;
		uint32_t offset = diagonal % query_sequence_length;
		if (chains_[offset].diagonal != diagonal
				|| !Connected(query_sequence, database_sequence,
						chains_[offset].query_end, *hits_it, cutoff)) {
			chains_[offset].query_end = hits_it->query_sequence_position;
			chains_[offset].diagonal = query_sequence_length - 1
					+ hits_it->database_sequence_position
					- hits_it->query_sequence_position;
			hits[n] = *hits_it;
			++n;
		} else {
			chains_[offset].query_end = hits_it->query_sequence_position;
		}
	}
	hits.resize(n);

}

bool ChainFilter::Connected(const AlphabetCoder::Code *query_sequence,
		const AlphabetCoder::Code *database_sequence, uint32_t chain_query_end,
		const Hit &hit, int cutoff) {
	if (hit.query_sequence_position <= chain_query_end + 1) {
		return true;
	}
	uint32_t q_offset = hit.query_sequence_position;
	uint32_t d_offset = hit.database_sequence_position;
	//assert(q_offset > 0);
	//assert(d_offset > 0);
	//assert(query_sequence[0] == sequence_delimiter_);
	//assert(database_sequence[0] == sequence_delimiter_);
	int score = 0;
	int best_score = 0;
	int c = -cutoff;
	uint32_t number_letters = score_matrix_.GetNumberLetters();
	const int *score_matrix = score_matrix_.GetMatrix();

	do {
		--d_offset;
		--q_offset;
		if (database_sequence[d_offset] == sequence_delimiter_
				|| query_sequence[q_offset] == sequence_delimiter_
				|| q_offset <= chain_query_end) {
			break;
		}
		AlphabetCoder::Code d = database_sequence[d_offset];
		AlphabetCoder::Code q = query_sequence[q_offset];
		score += score_matrix[d * number_letters + q];
		if (best_score < score) {
			best_score = score;
			c = best_score - cutoff;
		}
	} while (score > c);
	if (q_offset <= chain_query_end) {
		return true;
	} else {
		return false;
	}
}
