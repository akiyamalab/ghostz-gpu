/*
 * seed_searcher_query_parameters.h
 *
 *  Created on: 2013/04/23
 *      Author: shu
 */

#ifndef SEED_SEARCHER_QUERY_PARAMETERS_H_
#define SEED_SEARCHER_QUERY_PARAMETERS_H_

#include <stdint.h>
#include "seed_searcher_database_parameters.h"
#include "alphabet_coder.h"
#include "queries.h"

class SeedSearcherQueryParameters {
public:
	struct HashPositionData {
		bool similarity_check;
		uint32_t query_id;
		uint32_t position;
	};
	struct HashPositionDataList {
		SeedSearcherDatabaseParameters::SequencesIndex::HashKey hash;
		std::vector<HashPositionData> position_data;
	};

	struct BuildParameters {
		ScoreMatrix score_matrix;
		std::vector<int>* ungapped_extension_cutoffs_ptr;
		std::vector<int>* gapped_extension_triggers_ptr;
		SeedSearcherDatabaseParameters::CommonParameters *common_parameters;
	};

	int Build(int thread_id, int number_threads,
			boost::mutex *next_query_id_mutex, boost::barrier *barrier,
			Queries &queries, BuildParameters &parameters);

	//int Build(Queries &queries, BuildParameters &parameters);

	AlphabetCoder::Code GetSequenceDelimiter();
	std::vector<int> &GetUngappedExtensionCutoffs();
	std::vector<int> &GetGappedExtensionTriggers();
	ScoreMatrix &GetScoreMatrix();
	uint32_t GetNumberOfHashPositionDataLists();
	std::vector<AlphabetCoder::Code *> &GetQuerySequences();
	std::vector<HashPositionDataList> &GetHashPositionDataLists();

private:
	static const uint32_t kQueryIdsIncrement = 1 << 6;
	struct TempHashPositionData {
		SeedSearcherDatabaseParameters::SequencesIndex::HashKey hash;
		HashPositionData data;
	};

	typedef std::vector<TempHashPositionData> TempHashPositionDataBucket;
	typedef std::vector<TempHashPositionDataBucket> TempHashPositionDataBucketList;

	int PopQueryIds(size_t number_queries, boost::mutex *next_query_id_mutex,
			std::vector<uint32_t> *query_ids);

	AlphabetCoder::Code sequence_delimiter_;
	uint32_t number_using_hash_position_data_list_;
	uint32_t next_query_id_;
	std::vector<uint32_t> hash_position_data_offsets_;
	std::vector<uint32_t> seed_query_ids_;
	std::vector<uint32_t> seed_query_positions_;
	std::vector<HashPositionDataList> hash_position_data_lists_;
	std::vector<TempHashPositionDataBucketList> temp_hash_position_data_bucket_lists_;
	std::vector<uint32_t> number_using_hash_position_data_lists_;
	std::vector<uint32_t> number_using_hash_position_data_;
	std::vector<int>* ungapped_extension_cutoffs_ptr_;
	std::vector<int>* gapped_extension_triggers_ptr_;
	ScoreMatrix score_matrix_;
	std::vector<AlphabetCoder::Code *> query_sequences_;
};

#endif /* SEED_SEARCHER_QUERY_PARAMETERS_H_ */
