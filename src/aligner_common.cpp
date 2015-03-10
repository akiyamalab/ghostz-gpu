/*
 * aligner_common.cpp
 *
 *  Created on: Aug 8, 2014
 *      Author: shu
 */

#include <vector>
#include "aligner_common.h"
using namespace std;

int AlignerCommon::AlignmentComp(const int a_score, const Coordinate &a_start,
		const Coordinate &a_end, const int b_score, const Coordinate &b_start,
		const Coordinate &b_end) {
	if (a_score > b_score) {
		return 1;
	} else if (a_score < b_score) {
		return -1;
	}
	uint32_t a_q_length = a_end.query_position - a_start.query_position;
	uint32_t b_q_length = b_end.query_position - b_start.query_position;

	if (a_q_length > b_q_length) {
		return 1;
	} else if (a_q_length < b_q_length) {
		return -1;
	}

	uint32_t a_s_length = a_end.database_position - a_start.database_position;
	uint32_t b_s_length = b_end.database_position - b_start.database_position;

	if (a_s_length > b_s_length) {
		return 1;
	} else if (a_s_length < b_s_length) {
		return -1;
	}

#if 1
	if (a_start.database_position < b_start.database_position) {
		return 1;
	} else {
		return -1;
	}
	if (a_start.query_position < b_start.query_position) {
		return 1;
	} else {
		return -1;
	}
#endif

	return 0;
}

bool AlignerCommon::IsContained(Coordinate &a_start, Coordinate &a_end,
		Coordinate &b_start, Coordinate &b_end) {
	if ((a_start.query_position == b_start.query_position
			&& a_start.database_position == b_start.database_position) // overlapped start
			|| (a_end.query_position == b_end.query_position
					&& a_end.database_position == b_end.database_position) // overlapped end
			|| (a_start.query_position <= b_start.query_position // a contain b
			&& b_end.query_position <= a_end.query_position
					&& a_start.database_position <= b_start.database_position
					&& b_end.database_position <= a_end.database_position)) {
		return true;
	} else {
		return false;
	}
}

AlphabetCoder::Code AlignerCommon::GetSequenceDelimiter(SequenceType &type) {
	AlphabetCoder coder(type);
	return coder.GetMaxCode() + 1;
}

void AlignerCommon::BuildQueriesParameters(AligningCommonParameters &parameters,
		Queries::Parameters &queries_parameters) {
	Statistics statistics(*parameters.aligning_sequence_type_ptr);
	queries_parameters.filter = parameters.filter;
	queries_parameters.file_sequence_type_ptr =
			parameters.queries_file_sequence_type_ptr;
	queries_parameters.aligning_sequence_type_ptr =
			parameters.aligning_sequence_type_ptr;
	statistics.CalculateUngappedIdealKarlinParameters(parameters.score_matrix,
			&queries_parameters.ungapped_karlin_parameters);
	queries_parameters.chunk_size = parameters.queries_chunk_size;
	queries_parameters.score_matrix = parameters.score_matrix;
	queries_parameters.sequence_delimiter = AlignerCommon::GetSequenceDelimiter(
			*parameters.aligning_sequence_type_ptr);
	queries_parameters.number_threads = parameters.number_threads;
}


