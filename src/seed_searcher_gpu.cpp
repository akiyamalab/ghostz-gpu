/*
 * seed_searcher_gpu.cpp
 *
 *  Created on: Jul 10, 2014
 *      Author: shu
 */

#include "seed_searcher_gpu.h"
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
#include "distance_calculator_gpu.h"

#include <fstream>
#include <sstream>
#include "reduced_alphabet_coder.h"
#include "reduced_alphabet_file_reader.h"

using namespace std;

SeedSearcherGpu::SeedSearcherGpu() :
		used_clustering_(true), hit_similarity_threshold_rate_(0.0f), subsequence_length_(
				0), max_representation_member_distance_(0), distance_threshold_(
				0), database_parameters_ptr_(NULL), query_parameters_ptr_(NULL), sequences_index_(
				NULL), clustering_sequences_index_(NULL), query_seed_list_offsets_(
				NULL), hash_position_data_lists_(NULL)
#if 1
				, hash_hits_count(0), calculate_distance_count(0), filter_out_from_distance_count(
				0), representation_hit_count(0), representation_no_similarity_check_hit_count(
				0), member_hit_count(0), member_no_similarity_check_hit_count(
				0), ungapped_extension_count(0)
#endif
{
}
SeedSearcherGpu::~SeedSearcherGpu() {
}

uint32_t SeedSearcherGpu::Reset(QueryParameters &query_parameters,
		DatabaseParameters &database_parameters) {

	database_parameters_ptr_ = &database_parameters;
	query_parameters_ptr_ = &query_parameters;
	used_clustering_ = database_parameters_ptr_->IsClustering();
	hit_similarity_threshold_rate_ = 0.8f;
	sequences_index_ = &(database_parameters_ptr_->GetSequencesIndex());
	query_seed_list_offsets_ =
			query_parameters_ptr_->GetHashPositionDataOffsets();
	hash_position_data_lists_ =
			&(query_parameters_ptr_->GetHashPositionDataLists());
	subsequence_length_ = database_parameters_ptr_->GetSubsequenceLength();
	clustering_sequences_index_ =
			&(database_parameters_ptr_->GetClusteringSequencesIndex());

	max_representation_member_distance_ =
			database_parameters_ptr_->GetNumberMaxMismatchs();
	distance_threshold_ = subsequence_length_
			- hit_similarity_threshold_rate_ * subsequence_length_;
	return 0;
}

int SeedSearcherGpu::SearchSeedsFromRepresentativesForSimilarityFiltering(
		const uint32_t hash_position_data_list_id,
		DistanceCalculationSeedList *distance_calculation_seed_list) const {
	distance_calculation_seed_list->StartAddingSeed(hash_position_data_list_id);
#if 0
	if (3376739 == (*hash_position_data_lists_)[hash_position_data_list_id].hash) {
		cout << "r hp id " << hash_position_data_list_id << endl;
		/*
		 for (size_t i = 0;
		 i < hash_position_data_list.check_position_data.size(); ++i) {
		 cout << hash_position_data_list.check_position_data[i].query_id
		 << endl;
		 }
		 for (size_t i = 0;
		 i < hash_position_data_list.no_check_position_data.size();
		 ++i) {
		 cout << hash_position_data_list.no_check_position_data[i].query_id
		 << endl;
		 }
		 */
	}
#endif
	if (used_clustering_) {
		SearchSeedsFromRepresentativesForSimilarityFiltering(
				(*hash_position_data_lists_)[hash_position_data_list_id],
				*clustering_sequences_index_, distance_calculation_seed_list);
	}
	distance_calculation_seed_list->FinishedAddingSeed();
	return 0;
}

int SeedSearcherGpu::SearchSeedsFromClusteringSubsequence(
		const uint32_t hash_position_data_list_ids_i,
		DistanceCalculationSeedList &distance_calculation_seed_list,
		std::vector<uint32_t> &passing_distance_calculation_i_list,
		UngappedExtensionWithTriggerSeedList *ungapped_extension_seed_list) const {
	const uint32_t hash_position_data_list_id =
			distance_calculation_seed_list.GetHashPositionDataListIds()[hash_position_data_list_ids_i];
	SeedSearcherGpu::QueryParameters::HashPositionDataList &hash_position_data_list =
			(*hash_position_data_lists_)[hash_position_data_list_id];
	bool completed_hash_position_data_list_flag =
			distance_calculation_seed_list.IsCompletedHashPositionDataList(
					hash_position_data_list_ids_i);
#if 0
	if (3376739 == hash_position_data_list.hash) {
		cout << "s hp id " << hash_position_data_list_id << " c "
		<< completed_hash_position_data_list_flag << endl;
		/*
		 for (size_t i = 0;
		 i < hash_position_data_list.check_position_data.size(); ++i) {
		 cout << hash_position_data_list.check_position_data[i].query_id
		 << endl;
		 }
		 for (size_t i = 0;
		 i < hash_position_data_list.no_check_position_data.size();
		 ++i) {
		 cout << hash_position_data_list.no_check_position_data[i].query_id
		 << endl;
		 }
		 */
	}
#endif
	if (used_clustering_) {
		SearchSeedsFromClusteringSubsequences(hash_position_data_list_ids_i,
				query_seed_list_offsets_[hash_position_data_list_id],
				hash_position_data_list, completed_hash_position_data_list_flag,
				*clustering_sequences_index_, distance_calculation_seed_list,
				max_representation_member_distance_, distance_threshold_,
				passing_distance_calculation_i_list,
				ungapped_extension_seed_list);
	}
	return 0;
}

int SeedSearcherGpu::SearchSeedsFromNonClusteringSubsequence(
		const uint32_t hash_position_data_list_ids_i,
		DistanceCalculationSeedList &distance_calculation_seed_list,
		UngappedExtensionWithTriggerSeedList *ungapped_extension_seed_list) const {
	const uint32_t hash_position_data_list_id =
			distance_calculation_seed_list.GetHashPositionDataListIds()[hash_position_data_list_ids_i];
	SeedSearcherGpu::QueryParameters::HashPositionDataList &hash_position_data_list =
			(*hash_position_data_lists_)[hash_position_data_list_id];
	bool completed_hash_position_data_list_flag =
			distance_calculation_seed_list.IsCompletedHashPositionDataList(
					hash_position_data_list_ids_i);

#if 0
	if (3376739 == hash_position_data_list.hash) {
		cout << "s hp id " << hash_position_data_list_id << " c "
		<< completed_hash_position_data_list_flag << endl;
		/*
		 for (size_t i = 0;
		 i < hash_position_data_list.check_position_data.size(); ++i) {
		 cout << hash_position_data_list.check_position_data[i].query_id
		 << endl;
		 }
		 for (size_t i = 0;
		 i < hash_position_data_list.no_check_position_data.size();
		 ++i) {
		 cout << hash_position_data_list.no_check_position_data[i].query_id
		 << endl;
		 }
		 */
	}
#endif
	if (completed_hash_position_data_list_flag) {
		SearchSeedsFromNonClusteringSubsequences(
				query_seed_list_offsets_[hash_position_data_list_id],
				hash_position_data_list, *sequences_index_,
				ungapped_extension_seed_list);
	}

	return 0;
}

int SeedSearcherGpu::SearchSeedsFromRepresentativesForSimilarityFiltering(
		const SeedSearcherGpuQueryParameters::HashPositionDataList &hash_position_data_list,
		const DatabaseParameters::ClusteringSequencesIndex &clustering_sequences_index,
		DistanceCalculationSeedList *distance_calculation_seed_list) const {
#if 0
// debug ///////////////////
	uint32_t x_query_position = 18;
	uint32_t x_db_position = 30138028;
//cout << "hash : " << hash_position_data_list.hash << endl;
////////////////////////////
#endif
	size_t representatives_length = 0;
	const SeedSearcherDatabaseParameters::ClusteringSequencesIndex::RepresentativeData *representatives;
	if (!clustering_sequences_index.GetRepresentatives(
			hash_position_data_list.hash, &representatives,
			&representatives_length)) {
		for (size_t representatives_i = 0;
				representatives_i < representatives_length;
				++representatives_i) {
			uint32_t representative_position =
					representatives[representatives_i].second;
			size_t number_query_position_data =
					hash_position_data_list.check_position_data.size();
			for (size_t query_position_data_i = 0;
					query_position_data_i < number_query_position_data;
					++query_position_data_i) {
				const SeedSearcherGpuQueryParameters::HashPositionData &query_position_data =
						hash_position_data_list.check_position_data[query_position_data_i];
#if DEBUG
				++calculate_distance_count;
#endif
				distance_calculation_seed_list->AddSeed(
						query_position_data.position, representative_position);
#if 0
				if (x_db_position == representative_position) {
					cout << "dis r ok. cluster id : "
					<< representatives[representatives_i].first
					<< ", rep : "
					<< representatives[representatives_i].second
					<< "h vale : " << hash_position_data_list.hash
					<< " q " << query_position_data.query_id << ", "
					<< query_position_data.position << " d "
					<< representative_position << endl;
				}
#endif
			}
		}
	}
	return 0;
}

int SeedSearcherGpu::SearchSeedsFromClusteringSubsequences(
		const uint32_t hash_position_data_list_ids_i,
		const uint32_t query_seed_list_offset,
		const SeedSearcherGpuQueryParameters::HashPositionDataList &hash_position_data_list,
		const bool completed_hash_position_data_list_flag,
		const DatabaseParameters::ClusteringSequencesIndex &clustering_sequences_index,
		DistanceCalculationSeedList &distance_calculation_seed_list,
		const Distance max_representation_member_distance,
		const Distance distance_threshold,
		std::vector<uint32_t> &passing_distance_calculation_i_list,
		UngappedExtensionWithTriggerSeedList *ungapped_extension_seed_list) const {

#if 0
// debug ///////////////////
	uint32_t x_query_position = 18;
	uint32_t x_db_position = 623672280;
//cout << "hash : " << hash_position_data_list.hash << endl;
////////////////////////////
#endif

	size_t numer_calculated_seeds = 0;
	if (hash_position_data_list_ids_i == 0) {
		numer_calculated_seeds =
				distance_calculation_seed_list.GetNumberCalculatedSeeds();
	}
	size_t seed_offset = distance_calculation_seed_list.GetSeedOffset(
			hash_position_data_list_ids_i);
	size_t seed_end = distance_calculation_seed_list.GetSeedOffset(
			hash_position_data_list_ids_i + 1);
	if (distance_calculation_seed_list.GetDistancesLength() < seed_end) {
		seed_end = distance_calculation_seed_list.GetDistancesLength();
	}
	size_t number_seeds = seed_end - seed_offset;
	size_t current_seed_i = 0;
	size_t representatives_length = 0;

	const SeedSearcherDatabaseParameters::ClusteringSequencesIndex::RepresentativeData *representatives;
	if (!clustering_sequences_index.GetRepresentatives(
			hash_position_data_list.hash, &representatives,
			&representatives_length)) {
#if DEBUG
		hash_hits_count += representatives_length
		* hash_position_data_list.position_data.size();
#endif

		for (size_t representatives_i = 0;
				representatives_i < representatives_length;
				++representatives_i) {
			uint32_t cluster_id = representatives[representatives_i].first;
			uint32_t representative_position =
					representatives[representatives_i].second;
			passing_distance_calculation_i_list.clear();
#if 0
			if (x_db_position == representative_position) {
				cout << "r ok" << endl;
			}
#endif
			const vector<SeedSearcherGpu::QueryParameters::HashPositionData> &check_query_position_list =
					hash_position_data_list.check_position_data;
			size_t number_query_position_data =
					check_query_position_list.size();
			int offset = seed_offset + (int) current_seed_i
					- (int) numer_calculated_seeds;
			size_t end = number_seeds + numer_calculated_seeds;
			size_t query_seed_offset = query_seed_list_offset;
			for (size_t i = 0; i < number_query_position_data;
					++i, ++current_seed_i) {
				if (current_seed_i >= numer_calculated_seeds
						&& current_seed_i < end) {
#if 0
					if ((offset + i)
							>= distance_calculation_seed_list.GetDistancesLength()) {
						cout << "offset " << offset << endl;
						cout << "i " << i << endl;
						cout << "length "
						<< distance_calculation_seed_list.GetDistancesLength()
						<< endl;
						cout << "seed offset " << seed_offset << endl;
						cout << "number_seeds " << number_seeds << endl;
						cout << "current_seed_i " << current_seed_i << endl;
					}
					assert(
							(offset + i)
							< distance_calculation_seed_list.GetDistancesLength());
					assert(
							check_query_position_list[i].query_id
							== distance_calculation_seed_list.host_query_ids_[offset
							+ i]);
					assert(
							check_query_position_list[i].position
							== distance_calculation_seed_list.host_query_concatenated_positions_[offset
							+ i]);
					assert(
							representative_position
							== distance_calculation_seed_list.host_database_positions_[offset
							+ i]);
#endif
					int query_representation_distance =
							distance_calculation_seed_list.GetDistances()[offset
									+ i];
					if (distance_threshold >= query_representation_distance) {
						ungapped_extension_seed_list->AddSeed(
								query_seed_offset + i, representative_position);
#if 0
						if (x_db_position == representative_position) {
							cout << "sim r ok. cluster id : " << cluster_id
							<< ", rep : " << representative_position
							<< "h vale : " << hash_position_data_list.hash
							<< " q "
							<< distance_calculation_seed_list.GetQueryIds()[current_seed_i]
							<< ", "
							<< distance_calculation_seed_list.GetQueryConcatenatedPositions()[current_seed_i]
							<< " d " << representative_position << endl;
						}
#endif
					}

					if (IsSufficientSimilarityCluster(
							query_representation_distance,
							max_representation_member_distance,
							distance_threshold)) {
						passing_distance_calculation_i_list.push_back(i);

#if 0
						if (x_db_position == representative_position) {
							cout << "sim r pass. cluster id : " << cluster_id
							<< ", rep : " << representative_position
							<< "h vale : " << hash_position_data_list.hash
							<< " q "
							<< distance_calculation_seed_list.GetQueryIds()[current_seed_i]
							<< ", "
							<< distance_calculation_seed_list.GetQueryConcatenatedPositions()[current_seed_i]
							<< " d " << representative_position << endl;
						}
#endif
					}
				}
			}
			if (completed_hash_position_data_list_flag) {
				const vector<SeedSearcherGpu::QueryParameters::HashPositionData> &no_check_query_position_list =
						hash_position_data_list.no_check_position_data;
				size_t query_seed_offset = query_seed_list_offset
						+ hash_position_data_list.check_position_data.size();
				for (size_t no_check_query_position_data_i = 0;
						no_check_query_position_data_i
								< no_check_query_position_list.size();
						++no_check_query_position_data_i) {

					ungapped_extension_seed_list->AddSeed(
							query_seed_offset + no_check_query_position_data_i,
							representative_position);

#if 0
					if (x_db_position == representative_position) {
						cout << "no sim r ok. cluster id : " << cluster_id
						<< ", rep : " << representative_position
						<< "h vale : " << hash_position_data_list.hash
						<< " q "
						<< no_check_query_position_data.query_id << ", "
						<< no_check_query_position_data.position
						<< " d " << representative_position << endl;
					}
#endif
				}
			}

			const SeedSearcherDatabaseParameters::ClusteringSequencesIndex::Position *members;
			size_t members_length = 0;
			if ((!passing_distance_calculation_i_list.empty()
					|| completed_hash_position_data_list_flag)
					&& !clustering_sequences_index.GetMembers(cluster_id,
							&members, &members_length)) {

				size_t query_seed_offset = query_seed_list_offset;
				for (size_t pass_distance_calculation_i_list_i = 0;
						pass_distance_calculation_i_list_i
								< passing_distance_calculation_i_list.size();
						++pass_distance_calculation_i_list_i) {
					ungapped_extension_seed_list->AddSeed(
							query_seed_offset
									+ passing_distance_calculation_i_list[pass_distance_calculation_i_list_i],
							members, members_length);
#if 0
					if (x_db_position == member_position) {
						cout << "sim m ok. cluster id : " << cluster_id
						<< ", rep : " << representative_position
						<< "h vale : "
						<< hash_position_data_list.hash << " q "
						<< distance_calculation_seed_list.GetQueryIds()[seed_i]
						<< ", "
						<< distance_calculation_seed_list.GetQueryConcatenatedPositions()[seed_i]
						<< " d " << member_position << endl;
					}
#endif
				}

				if (completed_hash_position_data_list_flag) {
					const vector<
							SeedSearcherGpu::QueryParameters::HashPositionData> &no_check_query_position_list =
							hash_position_data_list.no_check_position_data;
					query_seed_offset =
							query_seed_list_offset
									+ hash_position_data_list.check_position_data.size();
					for (size_t no_check_query_position_data_i = 0;
							no_check_query_position_data_i
									< no_check_query_position_list.size();
							++no_check_query_position_data_i) {

						ungapped_extension_seed_list->AddSeed(
								query_seed_offset
										+ no_check_query_position_data_i,
								members, members_length);

#if 0
						if (x_db_position == member_position) {
							cout << "no sim m ok. cluster id : "
							<< cluster_id << ", rep : "
							<< representative_position
							<< "h vale : "
							<< hash_position_data_list.hash << " q "
							<< no_check_query_position_data.query_id
							<< ", "
							<< no_check_query_position_data.position
							<< " d " << member_position << endl;
						}

#endif

					}
				}
			}
		}
	}

	return 0;
}

int SeedSearcherGpu::SearchSeedsFromNonClusteringSubsequences(
		const uint32_t query_seed_list_offset,
		const SeedSearcherGpuQueryParameters::HashPositionDataList &hash_position_data_list,
		const DatabaseParameters::SequencesIndex &sequences_index,
		UngappedExtensionWithTriggerSeedList *ungapped_extension_seed_list) const {
#if 0
// debug ///////////////////
	uint32_t x_query_position = 18;
	uint32_t x_db_position = 623672280;
////////////////////////////
#endif

	const SeedSearcherDatabaseParameters::SequencesIndex::Value *values = NULL;
	size_t values_length = 0;
	if (!sequences_index.GetValues(hash_position_data_list.hash, &values,
			&values_length)) {
#if DEBUG
		hash_hits_count += values_length
		* hash_position_data_list.check_position_data.size();
#endif
		const vector<SeedSearcherGpu::QueryParameters::HashPositionData> &check_query_position_list =
				hash_position_data_list.check_position_data;
		size_t query_seed_offset = query_seed_list_offset;
		for (size_t check_query_position_data_i = 0;
				check_query_position_data_i < check_query_position_list.size();
				++check_query_position_data_i) {

			ungapped_extension_seed_list->AddSeed(
					query_seed_offset + check_query_position_data_i, values,
					values_length);
#if 0
			if (x_db_position == values[values_i]) {
				cout << "e ok hash vale : " << hash_position_data_list.hash
				<< " q " << check_query_position_data.query_id
				<< ", " << check_query_position_data.position
				<< " d " << values[values_i] << endl;
			}
#endif
		}

		const vector<SeedSearcherGpu::QueryParameters::HashPositionData> &no_check_query_position_list =
				hash_position_data_list.no_check_position_data;
		query_seed_offset = query_seed_list_offset
				+ check_query_position_list.size();
		for (size_t no_check_query_position_data_i = 0;
				no_check_query_position_data_i
						< no_check_query_position_list.size();
				++no_check_query_position_data_i) {

			ungapped_extension_seed_list->AddSeed(
					query_seed_offset + no_check_query_position_data_i, values,
					values_length);
#if 0
			if (x_db_position == values[values_i]) {
				cout << "e ok hash vale : " << hash_position_data_list.hash
				<< " q " << no_check_query_position_data.query_id
				<< ", " << no_check_query_position_data.position
				<< " d " << values[values_i] << endl;
			}
#endif
		}
	}
	return 0;
}

bool SeedSearcherGpu::IsSufficientSimilarityCluster(
		const int query_representation_distance,
		const int max_representation_member_distance,
		const int distance_threshold) const {

	int lower_bound = max(
			query_representation_distance - max_representation_member_distance,
			0);
	return distance_threshold >= lower_bound;
}

int SeedSearcherGpu::DumpSearchLog() {
	cout << "hash hits count " << hash_hits_count << endl;
	cout << "calculate_distance_count " << calculate_distance_count << endl;
	cout << "representation_hit_count " << representation_hit_count << endl;
	cout << "member_hit_count " << member_hit_count << endl;
	cout << "representation_no_similarity_check_hit_count "
			<< representation_no_similarity_check_hit_count << endl;
	cout << "member_no_similarity_check_hit_count "
			<< member_no_similarity_check_hit_count << endl;
	cout << "ungapped_extension_count " << ungapped_extension_count << endl;
	return 0;
}

