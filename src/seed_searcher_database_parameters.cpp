/*
 * seed_searcher_gpu_database_parameters.cpp
 *
 *  Created on: Jul 14, 2014
 *      Author: shu
 */


#include <assert.h>
#include <limits.h>
#include <float.h>
#include <algorithm>
#include <time.h>
#include <iostream>
#include <iomanip>
#include <tr1/unordered_map>
#include "k_mer_sequences_index.h"
#include "reduced_alphabet_k_mer_hash_function.h"
#include "reduced_alphabet_variable_hash_function.h"
#include "seed_searcher_common.h"
#include "seed_searcher_database_parameters.h"

#include <omp.h>

#ifdef TAU
#include <stdio.h>
#include <TAU.h>
#endif

using namespace std;

SeedSearcherDatabaseParameters::SeedSearcherDatabaseParameters() :
		clustering_(false), sequence_delimiter_(0), length_(0), subsequence_length_(
				0), cluster_threshold_rate_(0.0f), sequence_(NULL), redueced_code_map_(
				0) {
}

SeedSearcherDatabaseParameters::~SeedSearcherDatabaseParameters() {
}

int SeedSearcherDatabaseParameters::Build(AlphabetCoder::Code *sequence,
		uint32_t length, BuildParameters &parameters) {
	clustering_ = parameters.clustering;
	sequence_ = sequence;
	length_ = length;
	subsequence_length_ = parameters.subsequence_length;
	sequence_delimiter_ = parameters.sequence_delimiter;
	redueced_code_map_ = parameters.similairty_code_map;
	cluster_threshold_rate_ = 0.9f;
	int max_mismatch_count = GetNumberMaxMismatchs();
	ScoreMatrix &score_matrix = parameters.score_matrix;
	uint32_t number_letters = score_matrix.GetNumberLetters();
	int score_threshold_ = parameters.seed_threshold;
	uint32_t max_length = 8;
	vector<int> code_scores(parameters.max_indexing_code + 1);
	for (uint32_t i = 0; i < number_letters; ++i) {
		if (i < parameters.hash_code_map.size()) {
			AlphabetCoder::Code c = parameters.hash_code_map[i];
			if (c <= parameters.max_indexing_code) {
				code_scores[c] = std::max(code_scores[c],
						score_matrix.GetMatrix()[i * number_letters + i]);
			}
		}
	}

#ifdef DEBUG
	cout << "length : " << length_ << endl;
	cout << "subsequence length : " << subsequence_length_ << endl;
	cout << "code scores : " << endl;
	for (uint32_t i = 0; i < code_scores.size(); ++i) {
		cout << i << " : " << code_scores[i] << endl;
	}
#endif
	hash_function_ = SequencesIndex::HashFunction(parameters.max_indexing_code,
			score_threshold_, max_length, parameters.hash_code_map,
			code_scores);

	std::tr1::unordered_map<SequencesIndex::HashFunction::Hash, uint32_t> hash_map(
			0);
	std::tr1::unordered_map<SequencesIndex::HashFunction::Hash, uint32_t>::iterator find_it;
	vector<pair<SequencesIndex::HashKey, SequencesIndex::Value> > indexed_key_value_pairs(
			0);
	indexed_key_value_pairs.reserve(length_);
	uint32_t prev_hash_length = 0;
	for (uint32_t sequence_i = 0; sequence_i < length_; ++sequence_i) {
		SequencesIndex::HashKey h = 0;
		uint32_t hash_length = 0;
		int ret = hash_function_.CalculateHash(sequence_ + sequence_i, &h,
				&hash_length);
		if (!ret) {
			find_it = hash_map.find(h);
			if (find_it == hash_map.end()) {
				hash_map.insert(make_pair(h, 1));
			}
			uint32_t seed_position = SeedSearcherCommon::CalculateSeedPosition(
					sequence_i, hash_length);
			if ((prev_hash_length != 0)
					&& ((prev_hash_length - 1) >= hash_length)) {
				indexed_key_value_pairs.back() = make_pair(h, seed_position);
				prev_hash_length = hash_length;
			} else {
				indexed_key_value_pairs.push_back(make_pair(h, seed_position));
				prev_hash_length = hash_length;
			}
		} else {
			prev_hash_length = 0;
		}
	}

#if DEBUG
	cout << "key_value_pairs in base sequences index is "
	<< indexed_key_value_pairs.size() << endl;
#endif
	sequences_index_.Build(indexed_key_value_pairs.begin(),
			indexed_key_value_pairs.end(), hash_function_);

	if (clustering_) {
		vector<SequencesIndex::HashFunction::Hash> hashs(0);
		hashs.reserve(hash_map.size());
		for (std::tr1::unordered_map<SequencesIndex::HashFunction::Hash,
				uint32_t>::iterator hash_map_it = hash_map.begin();
				hash_map_it != hash_map.end(); ++hash_map_it) {
			hashs.push_back(hash_map_it->first);
		}
		hash_map.clear();
		const uint32_t hashs_size = hashs.size();
		ClusteringData init_clustering_data;
		init_clustering_data.clustering_state = kNonClusteringMember;
		init_clustering_data.representative_center = length;
		vector<ClusteringData> clustering_data(length, init_clustering_data);

#ifdef _OPENMP
		//omp_set_num_threads(parameters.number_threads);
#endif // _OPENMP
//#pragma omp parallel for schedule(dynamic, 10)
		for (uint32_t hashs_i = 0; hashs_i < hashs_size; ++hashs_i) {
			SequencesIndex::HashFunction::Hash hash = hashs[hashs_i];
			const KMerSequencesIndex::Value *values = NULL;
			size_t values_length = 0;
			uint32_t one_direction_length = subsequence_length_ / 2;
			if (!sequences_index_.GetValues(hash, &values, &values_length)) {
				for (size_t subsequences_i = 0; subsequences_i < values_length;
						++subsequences_i) {
					SequencesIndex::Value subsequence_center =
							values[subsequences_i];
					if (subsequence_center < one_direction_length) {
						clustering_data[subsequence_center].clustering_state =
								kNonClustering;
					} else {
						AlphabetCoder::Code *subsequence_sequence = sequence_
								+ subsequence_center - one_direction_length;
						for (size_t subsequence_i = 0;
								subsequence_i < subsequence_length_;
								++subsequence_i) {
							if (subsequence_sequence[subsequence_i]
									== sequence_delimiter_) {
								clustering_data[subsequence_center].clustering_state =
										kNonClustering;
								break;
							}
						}
					}
				}

				for (size_t subsequences_i = 0; subsequences_i < values_length;
						++subsequences_i) {
					SequencesIndex::Value representative_center =
							values[subsequences_i];
					if (clustering_data[representative_center].clustering_state
							!= kNonClusteringMember) {
						continue;
					}
					clustering_data[representative_center].clustering_state =
							kRepresentative;
					AlphabetCoder::Code *representative_sequence_center =
							sequence_ + representative_center;
					for (size_t member_i = subsequences_i + 1;
							member_i < values_length; ++member_i) {
						SequencesIndex::Value non_clustering_member_center =
								values[member_i];
						ClusteringData &non_clustering_member_data =
								clustering_data[non_clustering_member_center];
						if (non_clustering_member_data.clustering_state
								!= kNonClusteringMember) {
							continue;
						}
						AlphabetCoder::Code *non_clustering_member_sequence_center =
								sequence_ + non_clustering_member_center;

						int mismatch_count =
								SeedSearcherCommon::CalculateDistance(
										representative_sequence_center,
										non_clustering_member_sequence_center,
										subsequence_length_,
										redueced_code_map_);
						if (mismatch_count <= max_mismatch_count) {
							non_clustering_member_data.clustering_state =
									kMember;
							non_clustering_member_data.representative_center =
									representative_center;
						}
					}
				}
			}
		}

		vector<uint32_t> member_counts(clustering_data.size(), 0);
		for (size_t i = 0; i < clustering_data.size(); ++i) {
			if (clustering_data[i].clustering_state == kMember) {
#if DEBUG
				AlphabetCoder::Code *representative_sequence = sequence_
				+ clustering_data[i].representative_center;
				AlphabetCoder::Code *clustering_member_sequence = sequence_ + i;
				int mismatch_count = SeedSearcherCommon::CalculateDistance(
						representative_sequence, clustering_member_sequence,
						subsequence_length_, redueced_code_map_);

				assert(mismatch_count <= max_mismatch_count);
#endif
				++member_counts[clustering_data[i].representative_center];
			}
		}

		indexed_key_value_pairs.clear();
		vector<
				pair<ClusteringSequencesIndex::HashKey,
						ClusteringSequencesIndex::RepresentativeData> > indexed_hash_representations_pairs;
		vector<
				pair<ClusteringSequencesIndex::ClusterId,
						ClusteringSequencesIndex::Position> > indexed_cluster_id_position_pairs;
		std::tr1::unordered_map<ClusteringSequencesIndex::Position,
				ClusteringSequencesIndex::ClusterId> cluster_id_map;
		indexed_hash_representations_pairs.reserve(length_);
		indexed_cluster_id_position_pairs.reserve(length_);

		uint32_t number_clusters = 0;
		for (uint32_t hashs_i = 0; hashs_i < hashs_size; ++hashs_i) {
			SequencesIndex::HashFunction::Hash hash = hashs[hashs_i];
			const KMerSequencesIndex::Value *values = NULL;
			size_t values_length = 0;
			if (!sequences_index_.GetValues(hash, &values, &values_length)) {
				for (size_t values_i = 0; values_i < values_length;
						++values_i) {
					uint32_t sequence_center = values[values_i];
					if (clustering_data[sequence_center].clustering_state
							== kRepresentative
							&& member_counts[sequence_center] > 0) {
						indexed_hash_representations_pairs.push_back(
								make_pair(hash,
										make_pair(number_clusters,
												sequence_center)));
						cluster_id_map.insert(
								make_pair(sequence_center, number_clusters));
						++number_clusters;
					} else if (clustering_data[sequence_center].clustering_state
							== kMember) {
						std::tr1::unordered_map<
								ClusteringSequencesIndex::Position,
								ClusteringSequencesIndex::ClusterId>::iterator find_it =
								cluster_id_map.find(
										clustering_data[sequence_center].representative_center);
						assert(find_it != cluster_id_map.end());
						indexed_cluster_id_position_pairs.push_back(
								make_pair(find_it->second, sequence_center));
					} else {
						indexed_key_value_pairs.push_back(
								make_pair(hash, sequence_center));
					}
				}
			}
		}

#if DEBUG
		size_t representative_count = 0;
		size_t member_count = 0;
		size_t no_member_cluster_representative_count = 0;
		size_t non_clustering_subsequence_count = 0;

		for (size_t i = 0; i < clustering_data.size(); ++i) {
			if (clustering_data[i].clustering_state == kRepresentative
					&& member_counts[i] > 0) {
				++representative_count;
			}
			if (clustering_data[i].clustering_state == kRepresentative
					&& member_counts[i] == 0) {
				++no_member_cluster_representative_count;
			}
			if (clustering_data[i].clustering_state == kMember) {
				++member_count;
			}
			if (clustering_data[i].clustering_state == kNonClustering) {
				++non_clustering_subsequence_count;
			}
		}
		cout << "database length is                        " << length << endl;
		cout << "total number of subsequences is           "
				<< representative_count + no_member_cluster_representative_count
						+ member_count + non_clustering_subsequence_count
				<< endl;
		cout << "representative_count is                   "
				<< representative_count << endl;
		cout << "member_count is                           " << member_count
				<< endl;
		cout << "no_member_cluster_representative_count is "
				<< no_member_cluster_representative_count << endl;
		cout << "non_clustering_subsequence_count is       "
				<< non_clustering_subsequence_count << endl;
#endif
		sequences_index_.Build(indexed_key_value_pairs.begin(),
				indexed_key_value_pairs.end(), hash_function_);

		clusetering_sequences_index_.Build(
				indexed_hash_representations_pairs.begin(),
				indexed_hash_representations_pairs.end(),
				indexed_cluster_id_position_pairs.begin(),
				indexed_cluster_id_position_pairs.end(), hash_function_);
	}
	return 0;
}

int SeedSearcherDatabaseParameters::Build(AlphabetCoder::Code *sequence,
		uint32_t length, std::istream &is) {
	sequence_ = sequence;
	length_ = length;
	is.read((char *) &clustering_, sizeof(clustering_));
	is.read((char *) &sequence_delimiter_, sizeof(sequence_delimiter_));
	is.read((char *) &subsequence_length_, sizeof(subsequence_length_));
	sequences_index_.Load(is);
	clusetering_sequences_index_.Load(is);
	is.read((char *) &cluster_threshold_rate_, sizeof(cluster_threshold_rate_));
	size_t redueced_code_map_size = 0;
	is.read((char *) &redueced_code_map_size, sizeof(redueced_code_map_size));
	redueced_code_map_.resize(redueced_code_map_size);
	is.read((char *) &redueced_code_map_[0],
			sizeof(redueced_code_map_[0]) * redueced_code_map_.size());
	return 0;
}

int SeedSearcherDatabaseParameters::BuildCommonParameters(
		CommonParameters* common_parameters) {
	common_parameters->sequence_delimiter = sequence_delimiter_;
	common_parameters->subsequence_length = subsequence_length_;
	common_parameters->hash_function = hash_function_;
	common_parameters->redueced_code_map = redueced_code_map_;
	return 0;
}

int SeedSearcherDatabaseParameters::BuildCommonParameters(std::istream &is,
		CommonParameters* common_parameters) {
	is.read((char *) &common_parameters->sequence_delimiter,
			sizeof(common_parameters->sequence_delimiter));
	is.read((char *) &common_parameters->subsequence_length,
			sizeof(common_parameters->subsequence_length));
	
	size_t redueced_code_map_size = 0;
	is.read((char *) &redueced_code_map_size, sizeof(redueced_code_map_size));
	common_parameters->redueced_code_map.resize(redueced_code_map_size);
	is.read((char *) &common_parameters->redueced_code_map[0],
			sizeof(common_parameters->redueced_code_map[0]) * common_parameters->redueced_code_map.size());
	common_parameters->hash_function.Load(is);
	return 0;
}

int SeedSearcherDatabaseParameters::Save(std::ostream &os) {
	os.write((char *) &clustering_, sizeof(clustering_));
	os.write((char *) &sequence_delimiter_, sizeof(sequence_delimiter_));
	os.write((char *) &subsequence_length_, sizeof(subsequence_length_));
	sequences_index_.Save(os);
	clusetering_sequences_index_.Save(os);
	os.write((char *) &cluster_threshold_rate_,
			sizeof(cluster_threshold_rate_));
	size_t redueced_code_map_size = redueced_code_map_.size();
	os.write((char *) &redueced_code_map_size, sizeof(redueced_code_map_size));
	os.write((char *) &redueced_code_map_[0],
			sizeof(redueced_code_map_[0]) * redueced_code_map_.size());
	return 0;
}

int SeedSearcherDatabaseParameters::SaveCommonParameters(std::ostream &os) {
	os.write((char *) &sequence_delimiter_, sizeof(sequence_delimiter_));
	os.write((char *) &subsequence_length_, sizeof(subsequence_length_));
	size_t redueced_code_map_size = redueced_code_map_.size();
	os.write((char *) &redueced_code_map_size, sizeof(redueced_code_map_size));
	os.write((char *) &redueced_code_map_[0],
			sizeof(redueced_code_map_[0]) * redueced_code_map_.size());
	hash_function_.Save(os);
	return 0;
}

bool SeedSearcherDatabaseParameters::IsClustering() {
	return clustering_;
}

uint32_t SeedSearcherDatabaseParameters::GetSequenceLength() {
	return length_;
}

uint32_t SeedSearcherDatabaseParameters::GetSubsequenceLength() {
	return subsequence_length_;
}

AlphabetCoder::Code * SeedSearcherDatabaseParameters::GetSequence() {
	return sequence_;
}

SeedSearcherDatabaseParameters::SequencesIndex & SeedSearcherDatabaseParameters::GetSequencesIndex() {
	return sequences_index_;
}

int SeedSearcherDatabaseParameters::GetNumberMaxMismatchs() {
	return subsequence_length_ * (1 - cluster_threshold_rate_);
}

std::vector<AlphabetCoder::Code> &SeedSearcherDatabaseParameters::GetReduecedCodeMap() {
	return redueced_code_map_;
}

SeedSearcherDatabaseParameters::ClusteringSequencesIndex & SeedSearcherDatabaseParameters::GetClusteringSequencesIndex() {
	return clusetering_sequences_index_;
}

