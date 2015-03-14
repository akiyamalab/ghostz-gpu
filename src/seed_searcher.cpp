/*
 * seed_search.cpp
 *
 *  Created on: 2013/02/07
 *      Author: shu
 */

#include <vector>
#include <stdint.h>
#include <iostream>
#include <iomanip>
#include <limits.h>
#include <algorithm>
#include <iterator>
#include <assert.h>
#include <float.h>
#include <cmath>
#include "protein_type.h"
#include "alphabet_coder.h"
#include "seed_searcher.h"

#include <fstream>
#include <sstream>
#include "reduced_alphabet_coder.h"
#include "reduced_alphabet_file_reader.h"

#include <omp.h>

using namespace std;

SeedSearcher::SeedSearcher() :
		used_clustering_(true), sequence_delimiter_(0), hit_similarity_threshold_rate_(
				0.0f), number_of_threads_(1), subsequence_length_(0), max_representation_member_distance_(
				0), distance_threshold_(0), database_parameters_ptr_(NULL), query_parameters_ptr_(
				NULL), sequences_index_(NULL), clustering_sequences_index_(
				NULL), query_seed_list_offsets_(NULL), hash_position_data_lists_(
				NULL), query_sequences_(NULL), database_sequence_(NULL), score_matrix_(
				NULL), ungapped_extension_cutoffs_(NULL), gapped_extension_triggers_(
				NULL)
#if 1
				, hash_hits_count(0), calculate_distance_count(0), filter_out_from_distance_count(
				0), representation_hit_count(0), representation_no_similarity_check_hit_count(
				0), member_hit_count(0), member_no_similarity_check_hit_count(
				0), ungapped_extension_count(0)
#endif
{
}
SeedSearcher::~SeedSearcher() {
}

uint32_t SeedSearcher::Reset(SeedSearcherQueryParameters &query_parameters,
		SeedSearcherDatabaseParameters &database_parameters) {
	database_parameters_ptr_ = &database_parameters;
	query_parameters_ptr_ = &query_parameters;
	used_clustering_ = database_parameters_ptr_->IsClustering();
	hit_similarity_threshold_rate_ = 0.8f;
	sequences_index_ = &(database_parameters_ptr_->GetSequencesIndex());

	hash_position_data_lists_ =
			&(query_parameters_ptr_->GetHashPositionDataLists());
	subsequence_length_ = database_parameters_ptr_->GetSubsequenceLength();
	clustering_sequences_index_ =
			&(database_parameters_ptr_->GetClusteringSequencesIndex());

	max_representation_member_distance_ =
			database_parameters_ptr_->GetNumberMaxMismatchs();
	distance_threshold_ = subsequence_length_
			- hit_similarity_threshold_rate_ * subsequence_length_;
	query_sequences_ = &query_parameters.GetQuerySequences();
	database_sequence_ = database_parameters.GetSequence();
	sequence_delimiter_ = query_parameters.GetSequenceDelimiter();
	redueced_code_map_ = &database_parameters.GetReduecedCodeMap();
	score_matrix_ = &query_parameters.GetScoreMatrix();
	ungapped_extension_cutoffs_ =
			&query_parameters.GetUngappedExtensionCutoffs();
	gapped_extension_triggers_ = &query_parameters.GetGappedExtensionTriggers();
	return 0;
}

bool SeedSearcher::Search(const uint32_t hash_position_data_list_i,
		std::vector<uint32_t> &hitting_query_position_data_i_list,
		TempChainFilterSeedList *temp_ungapped_extension_seed_data) {

	const SeedSearcherQueryParameters::HashPositionDataList &hash_position_data_list =
			(*hash_position_data_lists_)[hash_position_data_list_i];

	SearchFromClusteringSubsequences(*query_sequences_, database_sequence_,
			sequence_delimiter_, hash_position_data_list,
			*clustering_sequences_index_, subsequence_length_,
			*redueced_code_map_, max_representation_member_distance_,
			distance_threshold_, *score_matrix_, *ungapped_extension_cutoffs_,
			*gapped_extension_triggers_, hitting_query_position_data_i_list,
			temp_ungapped_extension_seed_data);

	SearchFromNonClusteringSubsequences(*query_sequences_, database_sequence_,
			sequence_delimiter_, hash_position_data_list, *sequences_index_,
			*score_matrix_, *ungapped_extension_cutoffs_,
			*gapped_extension_triggers_, temp_ungapped_extension_seed_data);

	return 0;
}

int SeedSearcher::SearchFromClusteringSubsequences(
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
		TempChainFilterSeedList *temp_ungapped_extension_seed_data) {
	size_t representatives_length = 0;
	const SeedSearcherDatabaseParameters::ClusteringSequencesIndex::RepresentativeData *representatives;
	if (database_parameters_ptr_->IsClustering()
			&& !clustering_sequences_index.GetRepresentatives(
					hash_position_data_list.hash, &representatives,
					&representatives_length)) {
#if DEBUG 
		hash_hits_count += representatives_length
		* hash_position_data_list.position_data.size();
#endif

		for (size_t representatives_i = 0;
				representatives_i < representatives_length;
				++representatives_i) {
			hitting_query_position_data_i_list.clear();
			uint32_t representative_position =
					representatives[representatives_i].second;
			size_t number_query_position_data =
					hash_position_data_list.position_data.size();
			const AlphabetCoder::Code *representative_subsequence_center =
					database_sequence + representative_position;
			for (size_t query_position_data_i = 0;
					query_position_data_i < number_query_position_data;
					++query_position_data_i) {
				const SeedSearcherQueryParameters::HashPositionData &query_position_data =
						hash_position_data_list.position_data[query_position_data_i];
				int query_representation_distance = 0;
				if (query_position_data.similarity_check) {
#if DEBUG 
					++calculate_distance_count;
#endif

					query_representation_distance =
							SeedSearcherCommon::CalculateDistance(
									query_sequences[query_position_data.query_id]
											+ query_position_data.position,
									representative_subsequence_center,
									subsequence_length, redueced_code_map);

					if (IsSufficientSimilarityCluster(
							query_representation_distance,
							max_representation_member_distance,
							distance_threshold)) {
						hitting_query_position_data_i_list.push_back(
								query_position_data_i);
					}
				} else {
					hitting_query_position_data_i_list.push_back(
							query_position_data_i);
				}
				if (distance_threshold >= query_representation_distance) {
#if DEBUG 
					if (query_position_data.similarity_check) {
						++representation_hit_count;
					} else {
						++representation_no_similarity_check_hit_count;
					}
#endif
					int score =
							UngappedExtend(
									query_sequences[query_position_data.query_id],
									query_position_data.position,
									database_sequence, representative_position,
									sequence_delimiter, score_matrix,
									ungapped_extension_cutoffs[query_position_data.query_id],
									gapped_extension_triggers[query_position_data.query_id]);
					if (score
							> gapped_extension_triggers[query_position_data.query_id]) {
						AddSeed(query_position_data.query_id,
								query_position_data.position,
								representative_position,
								temp_ungapped_extension_seed_data);
					}
				}
			}

			const SeedSearcherDatabaseParameters::ClusteringSequencesIndex::Position *members;
			size_t members_length = 0;
			if (!hitting_query_position_data_i_list.empty()
					&& !clustering_sequences_index.GetMembers(
							representatives[representatives_i].first, &members,
							&members_length)) {
				for (size_t members_i = 0; members_i < members_length;
						++members_i) {
					for (size_t i = 0;
							i < hitting_query_position_data_i_list.size();
							++i) {
						const SeedSearcherQueryParameters::HashPositionData &query_position_data =
								hash_position_data_list.position_data[hitting_query_position_data_i_list[i]];
#if DEBUG 
						if (query_position_data.similarity_check) {
							++member_hit_count;
							++hash_hits_count;
						} else {
							++member_no_similarity_check_hit_count;
							++hash_hits_count;
						}
#endif
						int score =
								UngappedExtend(
										query_sequences[query_position_data.query_id],
										query_position_data.position,
										database_sequence, members[members_i],
										sequence_delimiter, score_matrix,
										ungapped_extension_cutoffs[query_position_data.query_id],
										gapped_extension_triggers[query_position_data.query_id]);
						if (score
								> gapped_extension_triggers[query_position_data.query_id]) {
							AddSeed(query_position_data.query_id,
									query_position_data.position,
									members[members_i],
									temp_ungapped_extension_seed_data);
						}
					}
				}
			}

		}

	}
	return 0;
}

int SeedSearcher::SearchFromNonClusteringSubsequences(
		const vector<AlphabetCoder::Code *> &query_sequences,
		const AlphabetCoder::Code *database_sequence,
		const AlphabetCoder::Code sequence_delimiter,
		const SeedSearcherQueryParameters::HashPositionDataList &hash_position_data_list,
		const SeedSearcherDatabaseParameters::SequencesIndex &sequences_index,
		const ScoreMatrix &score_matrix,
		const std::vector<int> &ungapped_extension_cutoffs,
		const std::vector<int> &gapped_extension_triggers,
		TempChainFilterSeedList *temp_ungapped_extension_seed_data) {

	const SeedSearcherDatabaseParameters::SequencesIndex::Value *values = NULL;
	size_t values_length = 0;
	if (!sequences_index.GetValues(hash_position_data_list.hash, &values,
			&values_length)) {
#if DEBUG
		hash_hits_count += values_length
		* hash_position_data_list.position_data.size();
#endif
		for (size_t values_i = 0; values_i < values_length; ++values_i) {
			for (std::vector<SeedSearcherQueryParameters::HashPositionData>::const_iterator queries_position_data_it =
					hash_position_data_list.position_data.begin();
					queries_position_data_it
							!= hash_position_data_list.position_data.end();
					++queries_position_data_it) {
#if 0
				Dostamce similarity = CalculateSimilarity(
						query_sequences[queries_position_data_it->query_id]
						+ queries_position_data_it->query_position,
						database_sequence + values[values_i],
						subsequence_length / 2, redueced_code_map);
				if (similarity <= hit_similarity_threshold) {
					continue;
				}
#endif
				int score =
						UngappedExtend(
								query_sequences[queries_position_data_it->query_id],
								queries_position_data_it->position,
								database_sequence, values[values_i],
								sequence_delimiter, score_matrix,
								ungapped_extension_cutoffs[queries_position_data_it->query_id],
								gapped_extension_triggers[queries_position_data_it->query_id]);
				if (score
						> gapped_extension_triggers[queries_position_data_it->query_id]) {
					AddSeed(queries_position_data_it->query_id,
							queries_position_data_it->position,
							values[values_i],
							temp_ungapped_extension_seed_data);
				}
			}
		}
	}

	return 0;
}

void SeedSearcher::AddSeed(uint32_t query_id, uint32_t query_position,
		uint32_t database_position,
		TempChainFilterSeedList *temp_ungapped_extension_seed_data) const {
	uint32_t bucket_i = query_id % number_of_threads_;
	TempChainFilterSeedBucket &bucket =
			temp_ungapped_extension_seed_data->buckets[bucket_i];
	bucket.query_ids.push_back(query_id);
	bucket.query_positions.push_back(query_position);
	bucket.database_positions.push_back(database_position);
}

bool SeedSearcher::IsSufficientSimilarityCluster(
		const int query_representation_distance,
		const int max_representation_member_distance,
		const int distance_threshold) const {

	int lower_bound = max(
			query_representation_distance - max_representation_member_distance,
			0);
	return distance_threshold >= lower_bound;
}

int SeedSearcher::DumpSearchLog() {
#if DEBUG
	cout << "hash hits count " << hash_hits_count << endl;
	cout << "calculate_distance_count " << calculate_distance_count << endl;
	cout << "representation_hit_count " << representation_hit_count << endl;
	cout << "member_hit_count " << member_hit_count << endl;
	cout << "representation_no_similarity_check_hit_count "
			<< representation_no_similarity_check_hit_count << endl;
	cout << "member_no_similarity_check_hit_count "
			<< member_no_similarity_check_hit_count << endl;
	cout << "ungapped_extension_count " << ungapped_extension_count << endl;
#endif
	return 0;
}
