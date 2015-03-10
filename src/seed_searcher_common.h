/*
 * seed_searcher_common.h
 *
 *  Created on: 2014/05/13
 *      Author: shu
 */

#ifndef SEED_SEARCHER_COMMON_H_
#define SEED_SEARCHER_COMMON_H_

#include "distance_calculator.h"
#include "reduced_alphabet_variable_hash_function.h"

class SeedSearcherCommon {
public:
	typedef DistanceCalculator::Distance Distance;
	typedef ReducedAlphabetVariableHashFunction HashFunction;
	typedef HashFunction::Hash Hash;

	typedef struct {
		uint32_t query_sequence_position;
		uint32_t database_sequence_position;
		int k_mer_count;
	} Hit;

	static SeedSearcherCommon::Distance CalculateDistance(
			const AlphabetCoder::Code *subsequence0_center,
			const AlphabetCoder::Code *subsequence1_center,
			const uint32_t sequence_length,
			const std::vector<AlphabetCoder::Code> &redueced_code_map);

	static uint32_t CalculateSeedPosition(uint32_t hashed_sequence_start,
			uint32_t hashed_sequence_length);

private:
};

inline SeedSearcherCommon::Distance SeedSearcherCommon::CalculateDistance(
		const AlphabetCoder::Code *subsequence0_center,
		const AlphabetCoder::Code *subsequence1_center,
		const uint32_t sequence_length,
		const std::vector<AlphabetCoder::Code> &redueced_code_map) {
	DistanceCalculator distance_calculator;
	distance_calculator.SetQueries(subsequence0_center);
	distance_calculator.SetDatabase(subsequence1_center);
	distance_calculator.SetReducedCodeMap(&redueced_code_map[0]);
	return distance_calculator.CalculateDistance(0, 0, sequence_length);
}

#endif /* SEED_SEARCHER_COMMON_H_ */
