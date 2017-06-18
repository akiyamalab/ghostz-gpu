/*
 * aligner_presearch_thread.h
 *
 *  Created on: Sep 19, 2014
 *      Author: shu
 */

#ifndef ALIGNER_PRESEARCH_THREAD_H_
#define ALIGNER_PRESEARCH_THREAD_H_

class GappedExtender;

class AlignerPresearchThread {
public:
	typedef Database<SeedSearcher> DatabaseType;
	typedef AlignerCommon::AligningCommonParameters AligningParameters;
	typedef AlignerCommon::PresearchedResult PresearchedResult;
	typedef SeedSearcher::TempChainFilterSeedBucket TempChainFilterSeedBucket;
	typedef SeedSearcher::TempChainFilterSeedList TempChainFilterSeedList;

	struct ThreadSharedParameters {
		boost::mutex next_query_id_mutex;
		boost::mutex next_hash_position_data_list_id_mutex;
		boost::mutex next_chain_filter_query_id_mutex;
		uint32_t next_hash_position_data_list_id;
		uint32_t next_chain_filter_query_id;
		size_t number_hash_position_data_lists;
		size_t queries_concatenated_sequence_length;
		Queries *queries;
		std::vector<int> *ungapped_extension_cutoffs;
		std::vector<int> *gapped_extension_triggers;
		DatabaseType *database;
		AlphabetCoder::Code *database_concatenated_sequence;
		uint32_t database_concatenated_sequence_length;
		AligningParameters *parameters;
		SeedSearcher *seed_searcher;
		SeedSearcher::QueryParameters *seed_searcher_query_parameters;
		std::vector<TempChainFilterSeedList> *temp_chain_filter_seed_lists;
		std::vector<std::vector<SeedSearcher::Hit> > *chain_filter_seed_lists;
		std::vector<std::vector<PresearchedResult> > *results_list;
		boost::barrier *presearch_barrier;
		boost::barrier *all_barrier;
	};

	struct ThreadParameters {
		int thread_id;
		int gpu_id;
		int gpu_data_list_id;
		int gapped_extension_cutoff;
		size_t gpu_seeds_memory_size;
		ThreadSharedParameters *shared_parameters;
	};

	AlignerPresearchThread();
	virtual ~AlignerPresearchThread();
	void Run(ThreadParameters& thread_parameters);

protected:
	typedef AlignerCommon::AlignmentPositionLessPosition AlignmentPositionLessPosition;
	typedef AlignerCommon::PresearchedResultGreaterScore PresearchedResultGreaterScore;
	typedef AlignerCommon::ResultGreaterScore ResultGreaterScore;
	typedef AlignerCommon::AlignmentPosition AlignmentPosition;
	typedef AlignerCommon::Coordinate Coordinate;
	typedef AlignerCommon::Result Result;

	static const size_t kMaxThreadExtensionSeedListSize;
	static const size_t kHashPositionDataListIdsIncrement = 1 << 8;
	static const size_t kQueryIdsIncrement = 1 << 8;
	static const size_t kNumberOfExtensionSeedLists = 2;
	static const size_t kExtensionSizeOnCpu = 1 << 10;

	enum GpuRun {
		kDistanceCalculation = 0,
		kUngappedExtensionWithTrigger = 1,
		kUngappedExtension = 2,
		kGappedExtension = 3,
		kIdle = 4,
	};

	size_t GetSeedsMemorySize(size_t distance_calculation_size,
			size_t ungapped_extension_with_trigger_size, size_t extension_size);

	int PopHashPositionDataListIds(
			std::vector<uint32_t> *query_ids);
	int PopHashPositionDataListIdsForSeedSearch(
			std::vector<uint32_t> *query_ids);
	int PopQueryIdsForChainFilter(std::vector<uint32_t> *query_ids);

	void AddResults(DatabaseType &database, AligningParameters &parameters,
			std::vector<std::vector<AlignmentPosition> > &alignment_start_positions_in_database,
			std::vector<PresearchedResult> &added_results,
			std::vector<PresearchedResult> &results);
	void UpdateDatabaseData(ThreadSharedParameters* shared_parameters);
	void ResetRepresentativeSearchQueryData();
	void ResetSeedSearchQueryData();

	bool CanRunSeedSearch();

	int Extend(UngappedExtender &ungapped_extender,
			GappedExtender &gapped_extender, uint32_t query_id, Query *query,
			const AlphabetCoder::Code *database_concatenated_sequence,
			AligningParameters &parameters,
			std::vector<int> &ungapped_extension_cutoffs,
			int gapped_extension_cutoff,
			std::vector<SeedSearcher::Hit> &hits,
			std::vector<PresearchedResult> &results);

	void SetChainFilteringSeeds();


	ThreadParameters *thread_parameters_;

	std::vector<uint32_t> hash_position_data_list_ids;
	std::vector<uint32_t> seed_search_hash_position_data_list_ids;
	std::vector<uint32_t> chain_filter_query_ids;
	std::vector<std::pair<uint32_t, uint32_t> > temp_ungapped_extension_seed_data_;
};

#endif /* ALIGNER_PRESEARCH_THREAD_H_ */
