/*
 * distance_calculator.cpp
 *
 *  Created on: Jul 25, 2014
 *      Author: shu
 */

#include "distance_calculator.h"

DistanceCalculator::DistanceCalculator() :
		concatenated_query_sequence_(NULL), database_sequence_(NULL), redueced_code_map_(
				NULL) {

}

DistanceCalculator::~DistanceCalculator() {

}

int DistanceCalculator::SetQueries(const AlphabetCoder::Code *query_sequence) {
	concatenated_query_sequence_ = query_sequence;
	return 0;
}
int DistanceCalculator::SetDatabase(
		const AlphabetCoder::Code *database_sequence) {
	database_sequence_ = database_sequence;
	return 0;
}
int DistanceCalculator::SetReducedCodeMap(
		const AlphabetCoder::Code *redueced_code_map) {
	redueced_code_map_ = redueced_code_map;
	return 0;
}

DistanceCalculator::Distance DistanceCalculator::CalculateDistance(
		uint32_t query_concatenated_center, uint32_t database_center,
		uint32_t subsequence_length) const {
	const uint32_t foward_direction = subsequence_length / 2;
	Distance mismatch_count = 0;
	const AlphabetCoder::Code *subsequence0 = concatenated_query_sequence_
			+ query_concatenated_center - foward_direction;
	const AlphabetCoder::Code *subsequence1 = database_sequence_
			+ database_center - foward_direction;
	for (uint32_t i = 0; i < subsequence_length; ++i) {
#if 0
		std::cout << int(redueced_code_map[subsequence0[i]]) << ","
		<< int(redueced_code_map[subsequence1[i]]) << std::endl;
#endif
		mismatch_count +=
				(redueced_code_map_[subsequence0[i]]
						!= redueced_code_map_[subsequence1[i]]) ? 1 : 0;
	}
	return mismatch_count;
}

int DistanceCalculator::CalculateDistances(size_t size,
		uint32_t *query_concatenated_centers, uint32_t *database_centers,
		uint32_t subsequence_length, Distance *distances) const {
	for (size_t i = 0; i < size; ++i) {
		distances[i] = CalculateDistance(query_concatenated_centers[i],
				database_centers[i], subsequence_length);
	}
	return 0;
}
