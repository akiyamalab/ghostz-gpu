/*
 * seed_searcher_gpu_query_parameters.h
 *
 *  Created on: Jul 14, 2014
 *      Author: shu
 */

#ifndef SEED_SEARCHER_GPU_QUERY_PARAMETERS_H_
#define SEED_SEARCHER_GPU_QUERY_PARAMETERS_H_

#include <vector>
#include <boost/thread.hpp>
#include "alphabet_coder.h"
#include "queries.h"
#include "seed_searcher_database_parameters.h"

class SeedSearcherGpuQueryParameters {
public:

	struct HashPositionData {
		uint32_t query_id;
		uint32_t position;
	};

	struct HashPositionDataList {
		SeedSearcherDatabaseParameters::SequencesIndex::HashKey hash;
		std::vector<HashPositionData> check_position_data;
		std::vector<HashPositionData> no_check_position_data;
	};

	struct BuildParameters {
		uint32_t *query_sequence_offsets;
		SeedSearcherDatabaseParameters::CommonParameters *common_parameters;
	};

	int Build(int thread_id, int number_threads,
			boost::mutex *next_query_id_mutex, boost::barrier *barrier,
			Queries &queries, BuildParameters &parameters);

	AlphabetCoder::Code GetSequenceDelimiter();
	uint32_t GetNumberOfHashPositionDataLists();
	std::vector<HashPositionDataList> &GetHashPositionDataLists();
	uint32_t *GetHashPositionDataOffsets();
	uint32_t *GetSeedQueryIds() {
		return &seed_query_ids_[0];
	}

	uint32_t *GetSeedQueryPositions() {
		return &seed_query_positions_[0];
	}

private:
	static const uint32_t kQueryIdsIncrement = 1 << 6;
	struct TempHashPositionData {
		bool similarity_check;
		SeedSearcherDatabaseParameters::SequencesIndex::HashKey hash;
		HashPositionData data;
	};

	typedef std::vector<TempHashPositionData> TempHashPositionDataBucket;
	typedef std::vector<TempHashPositionDataBucket> TempHashPositionDataBucketList;

	int PopQueryIds(size_t number_queries, boost::mutex *next_query_id_mutex,
			std::vector<uint32_t> *query_ids);
	int BuildSubsequences(AlphabetCoder::Code *sequence, const uint32_t length,
			const AlphabetCoder::Code &sequence_delimiter,
			const uint32_t min_length);

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
};

#endif /* SEED_SEARCHER_GPU_QUERY_PARAMETERS_H_ */
