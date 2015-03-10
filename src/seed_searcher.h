/*
 * seed_searcher.h
 *
 *  Created on: 2013/02/07
 *      Author: shu
 */

#ifndef SEED_SEARCHER_H_
#define SEED_SEARCHER_H_

#include "alphabet_coder.h"
#include "reduced_alphabet_k_mer_hash_function.h"
#include "seed_searcher_database_parameters.h"
#include "seed_searcher_query_parameters.h"
#include "ungapped_extender.h"
#include <stdint.h>
#include <list>
#include <iostream>
#include <iomanip>
#include <map>
#include <tr1/memory>

class SeedSearcher {
public:
	typedef SeedSearcherQueryParameters QueryParameters;
	typedef SeedSearcherDatabaseParameters DatabaseParameters;

	typedef SeedSearcherCommon::Hit Hit;
	typedef SeedSearcherCommon::Distance Distance;

	struct TempChainFilterSeedBucket {
		std::vector<uint32_t> query_ids;
		std::vector<uint32_t> query_positions;
		std::vector<uint32_t> database_positions;
	};

	struct TempChainFilterSeedList {
		std::vector<TempChainFilterSeedBucket> buckets;
	};

	SeedSearcher();
	virtual ~SeedSearcher();

	void SetNumberThreads(uint32_t number_threads) {
		number_of_threads_ = number_threads;
	}

	uint32_t Reset(SeedSearcherQueryParameters &query_parameters,
			SeedSearcherDatabaseParameters &database_parameters);

	bool Search(const uint32_t hash_position_data_list_i,
			std::vector<uint32_t> &hitting_query_position_data_i_list,
			TempChainFilterSeedList *temp_ungapped_extension_seed_data);
	int DumpSearchLog();

private:
	int SearchFromClusteringSubsequences(
			const std::vector<AlphabetCoder::Code *> &query_sequences,
			const AlphabetCoder::Code *database_sequence,
			const AlphabetCoder::Code sequence_delimiter,
			const SeedSearcherQueryParameters::HashPositionDataList &hash_position_data_list,
			const SeedSearcherDatabaseParameters::ClusteringSequencesIndex &clustering_sequences_index,
			const uint32_t subsequence_length,
			const std::vector<AlphabetCoder::Code> &redueced_code_map,
			const Distance max_representation_member_distance,
			const Distance distance_threshold, const ScoreMatrix &score_matrix,
			const std::vector<int> &ungapped_extension_cutoffs,
			const std::vector<int> &gapped_extension_triggers,
			std::vector<uint32_t> &hitting_query_position_data_i_list,
			TempChainFilterSeedList *temp_ungapped_extension_seed_data);
	int SearchFromNonClusteringSubsequences(
			const std::vector<AlphabetCoder::Code *> &query_sequences,
			const AlphabetCoder::Code *database_sequence,
			const AlphabetCoder::Code sequence_delimiter,
			const SeedSearcherQueryParameters::HashPositionDataList &hash_position_data_list,
			const SeedSearcherDatabaseParameters::SequencesIndex &sequences_index,
			const ScoreMatrix &score_matrix,
			const std::vector<int> &ungapped_extension_cutoffs,
			const std::vector<int> &gapped_extension_triggers,
			TempChainFilterSeedList *temp_ungapped_extension_seed_data);
	int UngappedExtend(const AlphabetCoder::Code *query_sequence,
			uint32_t query_position,
			const AlphabetCoder::Code *database_sequence,
			uint32_t database_position, AlphabetCoder::Code sequence_delimiter,
			const ScoreMatrix &score_matrix, int cutoff,
			int gapped_extension_trigger); // __attribute__((noinline));
	bool CheckSimilarity(AlphabetCoder::Code *query_sequence,
			uint32_t query_position, AlphabetCoder::Code *database_sequence,
			uint32_t database_position, uint32_t subsequence_length,
			std::vector<AlphabetCoder::Code> &redueced_code_map,
			Distance representation_member_similarity,
			Distance hit_similarity_threshold) const;
	void AddSeed(uint32_t query_id, uint32_t query_position,
			uint32_t database_position,
			TempChainFilterSeedList *temp_ungapped_extension_seed_data) const;
	bool IsSufficientSimilarityCluster(const int query_representation_distance,
			const int representation_member_distance,
			const int distance_threshold) const;

	bool used_clustering_;
	AlphabetCoder::Code sequence_delimiter_;
	float hit_similarity_threshold_rate_;
	uint32_t number_of_threads_;
	uint32_t subsequence_length_;
	int max_representation_member_distance_;
	Distance distance_threshold_;
	DatabaseParameters *database_parameters_ptr_;
	QueryParameters *query_parameters_ptr_;
	DatabaseParameters::SequencesIndex *sequences_index_;
	DatabaseParameters::ClusteringSequencesIndex *clustering_sequences_index_;
	uint32_t *query_seed_list_offsets_;
	std::vector<AlphabetCoder::Code> *redueced_code_map_;
	std::vector<SeedSearcherQueryParameters::HashPositionDataList> *hash_position_data_lists_;
	std::vector<AlphabetCoder::Code *> *query_sequences_;
	AlphabetCoder::Code *database_sequence_;
	ScoreMatrix *score_matrix_;
	std::vector<int> *ungapped_extension_cutoffs_;
	std::vector<int> *gapped_extension_triggers_;

#if 1
	uint64_t hash_hits_count;
	uint64_t calculate_distance_count;
	uint64_t filter_out_from_distance_count;
	uint64_t representation_hit_count;
	uint64_t representation_no_similarity_check_hit_count;
	uint64_t member_hit_count;
	uint64_t member_no_similarity_check_hit_count;
	uint64_t ungapped_extension_count;
#endif

};

inline int SeedSearcher::UngappedExtend(
		const AlphabetCoder::Code *query_sequence, uint32_t query_position,
		const AlphabetCoder::Code *database_sequence,
		uint32_t database_position, AlphabetCoder::Code sequence_delimiter,
		const ScoreMatrix &score_matrix, int cutoff,
		int gapped_extension_trigger) {
	//++ungapped_extension_count;

	int sum_score = 0;
	int score = 0;
	UngappedExtender::ExtendOneSideScoreOnly(
			query_sequence + query_position - 1,
			database_sequence + database_position - 1, sequence_delimiter, true,
			score_matrix, cutoff, &score);
	sum_score += score;

	if (sum_score <= gapped_extension_trigger) {
		UngappedExtender::ExtendOneSideScoreOnly(
				query_sequence + query_position,
				database_sequence + database_position, sequence_delimiter,
				false, score_matrix, cutoff, &score);
		sum_score += score;
	}

	return sum_score;
}
#endif /* SEED_SEARCHER_H_ */
