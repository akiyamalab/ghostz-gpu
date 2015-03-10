/*
 * seed_searcher_gpu_database_parameters.h
 *
 *  Created on: Jul 14, 2014
 *      Author: shu
 */

#ifndef SEED_SEARCHER_GPU_DATABASE_PARAMETERS_H_
#define SEED_SEARCHER_GPU_DATABASE_PARAMETERS_H_

#include <stdint.h>
#include <vector>
#include <map>
#include "k_mer_sequences_index.h"
#include "variable_hash_clustering_seuences_index.h"
#include "variable_hash_sequences_index.h"
#include "alphabet_coder.h"

class SeedSearcherDatabaseParameters {
public:
	typedef VariableHashClusteringSeuencesIndex ClusteringSequencesIndex;
	typedef VariableHashSequencesIndex SequencesIndex;

	struct BuildParameters {
		bool clustering;
		uint32_t subsequence_length;
		AlphabetCoder::Code max_indexing_code;
		AlphabetCoder::Code sequence_delimiter;
		int seed_threshold;
		float threshold;
		uint32_t number_threads;
		std::vector<AlphabetCoder::Code> hash_code_map;
		std::vector<AlphabetCoder::Code> similairty_code_map;
		ScoreMatrix score_matrix;
	};

	struct CommonParameters {
		AlphabetCoder::Code sequence_delimiter;
		uint32_t subsequence_length;
		std::vector<AlphabetCoder::Code> redueced_code_map;
		SequencesIndex::HashFunction hash_function;
	};

	struct SubsequenceData {
		uint32_t indexed_position;
	};

	SeedSearcherDatabaseParameters();
	virtual ~SeedSearcherDatabaseParameters();
	int Build(AlphabetCoder::Code *sequence, uint32_t length,
			BuildParameters &parameters);
	int Build(AlphabetCoder::Code *sequence, uint32_t length, std::istream &is);

	int BuildCommonParameters(CommonParameters* common_parameters);
	int BuildCommonParameters(std::istream &is,
			CommonParameters* common_parameters);

	int Save(std::ostream &os);
	int SaveCommonParameters(std::ostream &os);

	bool IsClustering();
	uint32_t GetSequenceLength();
	uint32_t GetSubsequenceLength();
	AlphabetCoder::Code* GetSequence();
	SequencesIndex &GetSequencesIndex();
	int GetNumberMaxMismatchs();
	std::vector<AlphabetCoder::Code> &GetReduecedCodeMap();
	ClusteringSequencesIndex &GetClusteringSequencesIndex();

private:
	enum ClusteringSubsequenceState {
		kNonClustering, kRepresentative, kMember, kNonClusteringMember
	};
	struct ClusteringData {
		ClusteringSubsequenceState clustering_state;
		uint32_t representative_center;
	};

	bool clustering_;
	AlphabetCoder::Code sequence_delimiter_;
	uint32_t length_;
	uint32_t subsequence_length_;
	float cluster_threshold_rate_;
	SequencesIndex::HashFunction hash_function_;
	AlphabetCoder::Code *sequence_;
	SequencesIndex sequences_index_;
	ClusteringSequencesIndex clusetering_sequences_index_;
	std::vector<AlphabetCoder::Code> redueced_code_map_;

};

#endif /* SEED_SEARCHER_GPU_DATABASE_PARAMETERS_H_ */
