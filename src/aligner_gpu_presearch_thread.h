/*
 * aligner_gpu_presearch_thread.h
 *
 *  Created on: Aug 8, 2014
 *      Author: shu
 */

#ifndef ALIGNER_GPU_PRESEARCH_THREAD_H_
#define ALIGNER_GPU_PRESEARCH_THREAD_H_

#include <vector>
#include <boost/thread/barrier.hpp>
#include "aligner.h"
#include "seed_searcher_gpu.h"
#include "gpu_stream_controller.h"
#include "aligner_common.h"
#include "aligner_gpu_data.h"
#include "distance_calculation_seed_list.h"
#include "ungapped_extension_with_trigger_seed_list.h"
#include "extension_seed_list.h"

class UngappedExtenderGpu;
class UngappedExtender;
class GappedExtender;

class AlignerGpuPresearchThread {
public:
	typedef Database<SeedSearcherGpu> DatabaseType;
	typedef AlignerCommon::AligningCommonParameters AligningParameters;
	typedef AlignerCommon::PresearchedResult PresearchedResult;
	struct TempChainFilterSeedBucket {
		std::vector<uint32_t> query_ids;
		std::vector<uint32_t> query_concatenated_positions;
		std::vector<uint32_t> database_positions;
	};

	struct TempChainFilterSeedList {
		std::vector<TempChainFilterSeedBucket> buckets;
	};

	struct ThreadSharedParameters {
		boost::mutex next_query_id_mutex;
		boost::mutex next_hash_position_data_list_id_mutex;
		boost::mutex next_chain_filter_query_id_mutex;
		uint32_t next_hash_position_data_list_id;
		uint32_t next_chain_filter_query_id;
		size_t number_hash_position_data_lists;
		size_t queries_concatenated_sequence_length;
		Queries *queries;
		AlphabetCoder::Code *queries_concatenated_sequence;
		std::vector<uint32_t> *query_sequence_offsets;
		std::vector<int> *ungapped_extension_cutoffs;
		std::vector<int> *gapped_extension_triggers;
		DatabaseType *database;
		AlphabetCoder::Code *database_concatenated_sequence;
		uint32_t database_concatenated_sequence_length;
		AligningParameters *parameters;
		SeedSearcherGpu *seed_searcher;
		SeedSearcherGpu::QueryParameters *seed_searcher_query_parameters;
		std::vector<AlignerGpuPresearchThread::TempChainFilterSeedList> *temp_chain_filter_seed_lists;
		std::vector<std::vector<SeedSearcherGpu::Hit> > *chain_filter_seed_lists;
		std::vector<AlignerGpuData> *gpu_data_list;
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
		GpuStreamController *gpu_stream_controller;
		ThreadSharedParameters *shared_parameters;
	};

	AlignerGpuPresearchThread();
	virtual ~AlignerGpuPresearchThread();
	void Run(ThreadParameters& thread_parameters);

private:
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
	static const uint32_t kNumberOfGpuResources =
			GpuStreamController::kMaxNumberOfConcurrentExecutions;
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

	int PopHashPositionDataListIdsForRepresentativeSearch(
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

	int RunDistanceCalculation(GpuStreamController &controller,
			HostSeedsMemory &host_seeds_memory);
	int RunUngappedExtensionWithTrigger(GpuStreamController &controller,
			HostSeedsMemory &host_seeds_memory);
	int RunUngappedExtension(size_t seed_lists_id,
			GpuStreamController &controller,
			HostSeedsMemory &host_seeds_memory);
	int RunGappedExtension(size_t seed_lists_id,
			GpuStreamController &controller,
			HostSeedsMemory &host_seeds_memory);
	int MoveDistanceCalculationResults(int gpu_controller_resource_id,
			GpuStreamController &controller,
			DistanceCalculationSeedList *distance_calculation_seed_list);
	int MoveUngappedExtensionWithTriggerResults(int gpu_controller_resource_id,
			GpuStreamController &controller,
			TempChainFilterSeedList *ungapped_extension_finished_seeds_list);

	int MoveGappedExtensionResults(int gpu_controller_resource_id,
			GpuStreamController &controller, size_t seed_list_id,
			std::vector<std::vector<AlignmentPosition> > &alignment_start_positions_in_database);

	int FinalizeUngappedExtension(int gpu_controller_resource_id,
			GpuStreamController &controller, size_t seed_list_id);

	bool CanRunSeedSearch();
	bool CanRunRepresentativeSearch();
	bool CanRunSearch();
	bool CanRunExtension(GpuStreamController& gpu_stream_controller);

	int ExtendOnCpu(UngappedExtender &ungapped_extender,
			GappedExtender &gapped_extender, uint32_t query_id, Query *query,
			const AlphabetCoder::Code *database_concatenated_sequence,
			AligningParameters &parameters,
			std::vector<int> &ungapped_extension_cutoffs,
			int gapped_extension_cutoff,
			std::vector<SeedSearcherGpu::Hit> &hits,
			std::vector<PresearchedResult> &results);

	int UngappedExtendBackSeedsOnCpu(size_t seed_list_id, size_t extension_size,
			UngappedExtender &ungapped_extender, Queries &queries,
			const AlphabetCoder::Code *queries_concatenated_sequence,
			const AlphabetCoder::Code *database_concatenated_sequence,
			AligningParameters &parameters,
			std::vector<int> &ungapped_extension_cutoffs);

	int GappedExtendBackSeedsOnCpu(size_t seed_list_id, size_t extension_size,
			GappedExtender &gapped_extender, Queries &queries,
			const AlphabetCoder::Code *queries_concatenated_sequence,
			const AlphabetCoder::Code *database_concatenated_sequence,
			AligningParameters &parameters, int gapped_extension_cutoff);

	void SetChainFilteringSeeds();

	HostSeedsMemory* GetSeedsMemory() {
		size_t id = available_seeds_memory_ids_.back();
		available_seeds_memory_ids_.pop_back();
		assert(!host_seeds_memories_[id].IsAllocated());
		return &host_seeds_memories_[id];
	}

	void ReleaseSeedsMemory(HostSeedsMemory &memory) {
		available_seeds_memory_ids_.push_back(memory.GetId());
		memory.Release();
	}

	void ReleaseSeedsMemory(DistanceCalculationSeedList &list) {
		HostSeedsMemory* memory = list.MoveHostMemory();
		ReleaseSeedsMemory(*memory);
	}

	void ReleaseSeedsMemory(UngappedExtensionWithTriggerSeedList &list) {
		HostSeedsMemory* memory = list.MoveHostMemory();
		ReleaseSeedsMemory(*memory);
	}

	void ReleaseSeedsMemory(ExtensionSeedList &list) {
		HostSeedsMemory* memory = list.GetSeedsMemory();
		available_seeds_memory_ids_.push_back(memory->GetId());
		list.ReleaseHostMemory();
	}

	size_t GetPreviousExtensionSeedListId(size_t id) {
		if (id == 0) {
			return 1;
		} else {
			return 0;
		}
	}

	GpuRun gpu_run_state_list_[kNumberOfGpuResources];
	int last_run_gpu_controller_resource_id_;
	size_t run_host_seeds_memory_ids_[kNumberOfGpuResources];
	HostSeedsMemory host_seeds_memories_[kNumberOfGpuResources];
	ThreadParameters *thread_parameters_;
	std::vector<uint32_t> available_seeds_memory_ids_;
	std::vector<uint32_t> representative_search_hash_position_data_list_ids;
	std::vector<uint32_t> seed_search_hash_position_data_list_ids;
	std::vector<uint32_t> chain_filter_query_ids;
	std::vector<std::pair<uint32_t, uint32_t> > temp_ungapped_extension_seed_data_;
	DistanceCalculationSeedList distance_calculation_seed_list_;
	UngappedExtensionWithTriggerSeedList ungapped_extension_with_trigger_seed_list_;
	ExtensionSeedList extension_seed_lists_[kNumberOfExtensionSeedLists];
	GpuStreamController::DistanceCalculationSeedsParameters last_distance_calculation_seeds_parameters_list_[kNumberOfGpuResources];
	GpuStreamController::UngappedExtensionWithTriggerSeedsParameters last_ungapped_extension_with_trigger_seeds_parameters_list_[kNumberOfGpuResources];
};

#endif /* ALIGNER_GPU_PRESEARCH_THREAD_H_ */
