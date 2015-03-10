/*
 * seed_searcher_gpu.h
 *
 *  Created on: Jul 10, 2014
 *      Author: shu
 */

#ifndef SEED_SEARCHER_GPU_H_
#define SEED_SEARCHER_GPU_H_

#include "alphabet_coder.h"
#include "reduced_alphabet_k_mer_hash_function.h"
#include "seed_searcher_database_parameters.h"
#include "seed_searcher_gpu_query_parameters.h"
#include "distance_calculation_seed_list.h"
#include "ungapped_extension_with_trigger_seed_list.h"
#include <stdint.h>
#include <list>
#include <iostream>
#include <iomanip>
#include <map>
#include <tr1/memory>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/system/cuda/experimental/pinned_allocator.h>
#include <thrust/copy.h>

class SeedSearcherGpu {
public:
	typedef SeedSearcherGpuQueryParameters QueryParameters;
	typedef SeedSearcherDatabaseParameters DatabaseParameters;

	typedef SeedSearcherCommon::Distance Distance;
	typedef SeedSearcherCommon::Hit Hit;

	SeedSearcherGpu();
	virtual ~SeedSearcherGpu();

	uint32_t Reset(QueryParameters &query_parameters,
			DatabaseParameters &database_parameters);

	int SearchSeedsFromRepresentativesForSimilarityFiltering(
			const uint32_t hash_position_data_list_id,
			DistanceCalculationSeedList *distance_calculation_seed_list) const;
	int SearchSeedsFromClusteringSubsequence(
			const uint32_t hash_position_data_list_ids_i,
			DistanceCalculationSeedList &distance_calculation_seed_list,
			std::vector<uint32_t> &passing_distance_calculation_i_list,
			UngappedExtensionWithTriggerSeedList *ungapped_extension_seed_list) const;

	int SearchSeedsFromNonClusteringSubsequence(
			const uint32_t hash_position_data_list_ids_i,
			DistanceCalculationSeedList &distance_calculation_seed_list,
			UngappedExtensionWithTriggerSeedList *ungapped_extension_seed_list) const;
	int DumpSearchLog();

private:
	int SearchSeedsFromRepresentativesForSimilarityFiltering(
			const SeedSearcherGpuQueryParameters::HashPositionDataList &hash_position_data_list,
			const DatabaseParameters::ClusteringSequencesIndex &clustering_sequences_index,
			DistanceCalculationSeedList *distance_calculation_seed_list) const;

	int SearchSeedsFromClusteringSubsequences(
			const uint32_t hash_position_data_list_ids_i,
			const uint32_t query_seed_offset,
			const SeedSearcherGpuQueryParameters::HashPositionDataList &hash_position_data_list,
			const bool completed_hash_position_data_list_flag,
			const DatabaseParameters::ClusteringSequencesIndex &clustering_sequences_index,
			DistanceCalculationSeedList &distance_calculation_seed_list,
			const Distance max_representation_member_distance,
			const Distance distance_threshold,
			std::vector<uint32_t> &passing_distance_calculation_i_list,
			UngappedExtensionWithTriggerSeedList *ungapped_extension_seed_list) const;

	int SearchSeedsFromNonClusteringSubsequences(
			const uint32_t query_seed_offset,
			const SeedSearcherGpuQueryParameters::HashPositionDataList &hash_position_data_list,
			const DatabaseParameters::SequencesIndex &sequences_index,
			UngappedExtensionWithTriggerSeedList *ungapped_extension_seed_list) const;

	bool CheckSimilarity(AlphabetCoder::Code *query_sequence,
			uint32_t query_position, AlphabetCoder::Code *database_sequence,
			uint32_t database_position, uint32_t subsequence_length,
			std::vector<AlphabetCoder::Code> &redueced_code_map,
			Distance representation_member_similarity,
			Distance hit_similarity_threshold) const;

	bool IsSufficientSimilarityCluster(const int query_representation_distance,
			const int representation_member_distance,
			const int distance_threshold) const;

	bool used_clustering_;
	float hit_similarity_threshold_rate_;
	uint32_t subsequence_length_;
	int max_representation_member_distance_;
	Distance distance_threshold_;
	DatabaseParameters *database_parameters_ptr_;
	QueryParameters *query_parameters_ptr_;
	DatabaseParameters::SequencesIndex *sequences_index_;
	DatabaseParameters::ClusteringSequencesIndex *clustering_sequences_index_;
	uint32_t *query_seed_list_offsets_;
	std::vector<SeedSearcherGpuQueryParameters::HashPositionDataList> *hash_position_data_lists_;

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

#endif /* SEED_SEARCHER_GPU_H_ */
