/*
 * distance_calculator.h
 *
 *  Created on: Jul 25, 2014
 *      Author: shu
 */

#ifndef DISTANCE_CALCULATOR_H_
#define DISTANCE_CALCULATOR_H_

#include <stdint.h>
#include "alphabet_coder.h"

class DistanceCalculator {
public:
	typedef char Distance;

	DistanceCalculator();
	virtual ~DistanceCalculator();

	int SetQueries(const AlphabetCoder::Code *query_sequence);
	int SetDatabase(const AlphabetCoder::Code *database_sequence);
	int SetReducedCodeMap(const AlphabetCoder::Code *redueced_code_map);
	Distance CalculateDistance(uint32_t query_concatenated_position,
			uint32_t database_position, uint32_t subsequence_length) const;
	int CalculateDistances(size_t size, uint32_t *query_concatenated_centers,
			uint32_t *database_centers, uint32_t subsequence_length,
			Distance *distances) const;
private:
	const AlphabetCoder::Code *concatenated_query_sequence_;
	const AlphabetCoder::Code *database_sequence_;
	const AlphabetCoder::Code *redueced_code_map_;
};

#endif /* DISTANCE_CALCULATOR_H_ */
