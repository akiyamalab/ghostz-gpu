/*
 * aligner_gpu_presearch_thread.cpp
 *
 *  Created on: Aug 8, 2014
 *      Author: shu
 */

#include <string>
#include <vector>
#include <queue>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <assert.h>
#include "limits.h"
#include "logger.h"
#include "alphabet_coder.h"
#include "chain_filter.h"
#include "edit_blocks.h"
#include "score_matrix_reader.h"
#include "score_matrix.h"
#include "sequence_type.h"
#include "protein_type.h"
#include "dna_type.h"
#include "ungapped_extender.h"
#include "ungapped_extender_gpu.h"
#include "gapped_extender.h"
#include "gapped_extender_gpu.h"
#include "seed_searcher_gpu.h"
#include "seed_searcher_query_parameters.h"
#include "reduced_alphabet_coder.h"
#include "aligner_gpu_data.h"
#include "cuda_common.h"
#include "gpu_stream_controller.h"
#include "aligner_gpu_presearch_thread.h"

#include <boost/thread.hpp>

using namespace std;

const size_t AlignerGpuPresearchThread::kMaxThreadExtensionSeedListSize = 1
		<< 20;

AlignerGpuPresearchThread::AlignerGpuPresearchThread() :
		last_run_gpu_controller_resource_id_(-1), thread_parameters_(NULL) {
	for (size_t i = 0; i < kNumberOfGpuResources; ++i) {
		gpu_run_state_list_[i] = kIdle;
		run_host_seeds_memory_ids_[i] = -1;
	}
}

AlignerGpuPresearchThread::~AlignerGpuPresearchThread() {
// TODO Auto-generated destructor stub
}

bool AlignerGpuPresearchThread::CanRunRepresentativeSearch() {
	return !distance_calculation_seed_list_.FilledBuffer()
			&& !PopHashPositionDataListIdsForRepresentativeSearch(
					&representative_search_hash_position_data_list_ids);
}

void AlignerGpuPresearchThread::Run(ThreadParameters& thread_parameters) {
	thread_parameters_ = &thread_parameters;
	int thread_id = thread_parameters_->thread_id;
	CUDA_CHECK_RETURN(cudaSetDevice(thread_parameters.gpu_id));
	ThreadSharedParameters *shared_parameters =
			thread_parameters_->shared_parameters;
	Queries &queries = *(shared_parameters->queries);
	AligningParameters& parameters = *(shared_parameters->parameters);
	if (thread_id < 0) {
		DatabaseType::PreloadTarget target = DatabaseType::kSequence
				| DatabaseType::kOffset
				| DatabaseType::kSeedSearcherCommonParameters
				| DatabaseType::kSeedSearcherParameters;
		thread_parameters.shared_parameters->database->Preload(0, target);
	} else {
		SeedSearcherGpu::QueryParameters::BuildParameters seed_search_query_parameters_build_parameter;
		seed_search_query_parameters_build_parameter.query_sequence_offsets =
				&(*shared_parameters->query_sequence_offsets)[0];
		seed_search_query_parameters_build_parameter.common_parameters =
				&shared_parameters->database->GetSeedSearcherCommonParameters();

		shared_parameters->seed_searcher_query_parameters->Build(thread_id,
				shared_parameters->parameters->number_threads,
				&shared_parameters->next_query_id_mutex,
				shared_parameters->presearch_barrier, queries,
				seed_search_query_parameters_build_parameter);
		if (thread_id == 0) {
			shared_parameters->number_hash_position_data_lists =
					shared_parameters->seed_searcher_query_parameters->GetNumberOfHashPositionDataLists();
		}
		shared_parameters->presearch_barrier->wait();

		if (thread_id < (int) shared_parameters->gpu_data_list->size()) {
			AlignerGpuData &gpu_data =
					(*shared_parameters->gpu_data_list)[thread_parameters_->gpu_data_list_id];
			gpu_data.SetGpuQueriesSequence(
					shared_parameters->queries_concatenated_sequence,
					shared_parameters->queries_concatenated_sequence_length);
			gpu_data.SetGpuUngappedExtensionCutoffs(
					&(*shared_parameters->ungapped_extension_cutoffs)[0],
					shared_parameters->queries->GetNumberOfSequences());

			gpu_data.SetGpuGappedExtensionTriggers(
					&(*shared_parameters->gapped_extension_triggers)[0],
					shared_parameters->queries->GetNumberOfSequences());
			std::vector<AlphabetCoder::Code> &reduced_code_map =
					shared_parameters->database->GetSeedSearcherCommonParameters().redueced_code_map;
			gpu_data.SetGpuReducedCodeMap(&reduced_code_map[0],
					reduced_code_map.size());
			gpu_data.SetGpuScoreMatrix(
					shared_parameters->parameters->score_matrix.GetMatrix(),
					shared_parameters->parameters->score_matrix.GetNumberLetters());

			size_t hash_position_data_length =
					shared_parameters->seed_searcher_query_parameters->GetHashPositionDataOffsets()[shared_parameters->number_hash_position_data_lists];
			uint32_t *query_seed_query_ids =
					shared_parameters->seed_searcher_query_parameters->GetSeedQueryIds();
			uint32_t *query_seed_positions =
					shared_parameters->seed_searcher_query_parameters->GetSeedQueryPositions();

			gpu_data.SetGpuQueryIds(&query_seed_query_ids[0],
					hash_position_data_length);

			gpu_data.SetGpuConcatenatedQuerySequencePositions(
					&query_seed_positions[0], hash_position_data_length);

		}
		shared_parameters->presearch_barrier->wait();
	}

	size_t distance_calculation_size = thread_parameters_->gpu_seeds_memory_size;
	size_t ungapped_extension_with_trigger_size =
			thread_parameters_->gpu_seeds_memory_size;
	size_t extension_size = std::min(kMaxThreadExtensionSeedListSize,
			thread_parameters_->gpu_seeds_memory_size / 2);
	size_t gpu_seeds_memory_size = GetSeedsMemorySize(distance_calculation_size,
			ungapped_extension_with_trigger_size, extension_size);
	for (size_t i = 0; i < kNumberOfGpuResources; ++i) {
		host_seeds_memories_[i].Init(i, gpu_seeds_memory_size);
		available_seeds_memory_ids_.push_back(i);
	}
	if (thread_id == 0) {
		cout << "gpu seed size : " << gpu_seeds_memory_size << endl;
	}
	distance_calculation_size = std::min(distance_calculation_size,
			gpu_seeds_memory_size);
	ungapped_extension_with_trigger_size = std::min(
			ungapped_extension_with_trigger_size, gpu_seeds_memory_size);
	extension_size = std::min(extension_size, gpu_seeds_memory_size / 2);

	temp_ungapped_extension_seed_data_.reserve(gpu_seeds_memory_size);
	last_run_gpu_controller_resource_id_ = -1;
	GpuStreamController gpu_stream_controller;
	gpu_stream_controller.Init(gpu_seeds_memory_size);
	thread_parameters_->gpu_stream_controller = &gpu_stream_controller;
	ChainFilter chain_filter(queries.GetQuery(0)->GetSequenceDelimiter(),
			parameters.score_matrix);
	TempChainFilterSeedList & temp_chain_filter_seed_list =
			thread_id >= 0 ?
					(*shared_parameters->temp_chain_filter_seed_lists)[thread_id] :
					(*shared_parameters->temp_chain_filter_seed_lists)[0];
	if (thread_id >= 0) {
		temp_chain_filter_seed_list.buckets.resize(
				shared_parameters->parameters->number_threads);
	}

	GappedExtender gapped_extender;
	vector<PresearchedResult> tmp_cpu_presearch_results;
	vector<std::vector<AlignmentPosition> > alignment_start_positions_in_database;
	uint32_t number_letters = parameters.score_matrix.GetNumberLetters();
	DistanceCalculatorGpu distance_calculator_gpu;
	assert(
			(*shared_parameters->gpu_data_list)[thread_parameters.gpu_data_list_id].GetGpuId()
					== thread_parameters.gpu_id);
	distance_calculator_gpu.SetQueries(
			(*shared_parameters->gpu_data_list)[thread_parameters.gpu_data_list_id].GetGpuQueriesSequence());
	distance_calculator_gpu.SetReducedCodeMap(
			(*shared_parameters->gpu_data_list)[thread_parameters.gpu_data_list_id].GetGpuReducedCodeMap());

	distance_calculator_gpu.SetSubsequenceLength(
			shared_parameters->database->GetSeedSearcherCommonParameters().subsequence_length);
	UngappedExtenderGpu ungapped_extender_gpu;
	ungapped_extender_gpu.SetQueries(
			shared_parameters->database->GetSequenceDelimiter(),
			(*shared_parameters->gpu_data_list)[thread_parameters.gpu_data_list_id].GetGpuQueriesSequence(),
			(*shared_parameters->gpu_data_list)[thread_parameters.gpu_data_list_id].GetGpuUngappedExtensionCutoffs(),
			(*shared_parameters->gpu_data_list)[thread_parameters.gpu_data_list_id].GetGpuGappedExtensionTriggers());
	ungapped_extender_gpu.SetScoreMatrix(
			(*shared_parameters->gpu_data_list)[thread_parameters.gpu_data_list_id].GetGpuScoreMatrix(),
			number_letters);

	ungapped_extender_gpu.SetQuerySeedDataList(
			(*shared_parameters->gpu_data_list)[thread_parameters.gpu_data_list_id].GetGpuQuerySeedQueryIds(),
			(*shared_parameters->gpu_data_list)[thread_parameters.gpu_data_list_id].GetGpuQuerySeedQueryPositions());

	GappedExtenderGpu gapped_extender_gpu;
	gapped_extender_gpu.SetQueries(
			shared_parameters->database->GetSequenceDelimiter(),
			(*shared_parameters->gpu_data_list)[thread_parameters.gpu_data_list_id].GetGpuQueriesSequence());
	gapped_extender_gpu.SetScoreParameters(
			(*shared_parameters->gpu_data_list)[thread_parameters.gpu_data_list_id].GetGpuScoreMatrix(),

			number_letters, shared_parameters->parameters->gap_open,
			shared_parameters->parameters->gap_extension,
			thread_parameters_->gapped_extension_cutoff);

	std::vector<uint32_t> passing_distance_calculation_i_list;
	UngappedExtender ungapped_extender;
	bool representative_search_finished = false;
	bool seed_search_finished = false;
// debug ///////
	bool gapped_extension_gpu_flag = true;
///////////////
	shared_parameters->all_barrier->wait();
	while (shared_parameters->database->GetChunkId()
			< (int) shared_parameters->database->GetNumberChunks()) {
		representative_search_finished = false;
		seed_search_finished = false;

		assert(available_seeds_memory_ids_.size() == kNumberOfGpuResources);
		assert(
				gpu_stream_controller.GetNumberOfAvailableResources()
						== kNumberOfGpuResources);
		for (size_t i = 0; i < kNumberOfGpuResources; ++i) {
			assert(host_seeds_memories_[i].IsAllocated() == false);
		}

		for (size_t i = 0; i < 2; ++i) {
			assert(
					distance_calculation_seed_list_.GetNumberHashPositionDataListIds()
							== 0);
		}

		for (size_t i = 0;
				i
						< thread_parameters.shared_parameters->parameters->number_threads;
				++i) {
			TempChainFilterSeedBucket &seed_bucket =
					temp_chain_filter_seed_list.buckets[i];
			assert(seed_bucket.query_ids.empty());
			assert(seed_bucket.query_concatenated_positions.empty());
			assert(seed_bucket.database_positions.empty());
		}

		if (thread_id == 0) {
			cout << "start to presearch against database chunk "
					<< shared_parameters->database->GetChunkId() << endl;
			shared_parameters->next_hash_position_data_list_id = 0;
			shared_parameters->next_chain_filter_query_id = 0;
			UpdateDatabaseData(shared_parameters);
		}
		shared_parameters->all_barrier->wait();

		if (0 <= thread_id
				&& thread_id < (int) shared_parameters->gpu_data_list->size()) {
			(*shared_parameters->gpu_data_list)[thread_parameters_->gpu_data_list_id].SetGpuDatabaseSequence(
					shared_parameters->database_concatenated_sequence,
					shared_parameters->database_concatenated_sequence_length);

		}
		shared_parameters->all_barrier->wait();
		if (thread_id < 0) {
			int next_chunk_id = shared_parameters->database->GetChunkId() + 1;
			if (next_chunk_id
					< (int) shared_parameters->database->GetNumberChunks()) {
				DatabaseType::PreloadTarget target = DatabaseType::kSequence
						| DatabaseType::kOffset
						| DatabaseType::kSeedSearcherParameters;
				thread_parameters.shared_parameters->database->Preload(
						next_chunk_id, target);
			} else {
				DatabaseType::PreloadTarget target = DatabaseType::kName
						| DatabaseType::kSequence | DatabaseType::kOffset;
				thread_parameters.shared_parameters->database->Preload(0,
						target);
			}
		} else {
			distance_calculator_gpu.SetDatabase(
					(*shared_parameters->gpu_data_list)[thread_parameters.gpu_data_list_id].GetGpuDatabaseSequence());
			ungapped_extender_gpu.SetDatabase(
					(*shared_parameters->gpu_data_list)[thread_parameters.gpu_data_list_id].GetGpuDatabaseSequence());
			gapped_extender_gpu.SetDatabase(
					(*shared_parameters->gpu_data_list)[thread_parameters.gpu_data_list_id].GetGpuDatabaseSequence());
			gpu_stream_controller.SetDistanceCalculator(
					&distance_calculator_gpu);
			gpu_stream_controller.SetUngappedExtender(&ungapped_extender_gpu);
			gpu_stream_controller.SetGappedExtender(&gapped_extender_gpu);

			if (alignment_start_positions_in_database.size()
					< shared_parameters->database->GetNumberSequencesInChunk()) {
				alignment_start_positions_in_database.resize(
						shared_parameters->database->GetNumberSequencesInChunk());
			}

			for (size_t i = 0; i < kNumberOfGpuResources; ++i) {
				gpu_run_state_list_[i] = kIdle;
				last_distance_calculation_seeds_parameters_list_[i].size = 0;
				last_ungapped_extension_with_trigger_seeds_parameters_list_[i].size =
						0;
			}

			//gpu_stream_controller.FinishedRun(); # for profile
			distance_calculation_seed_list_.SetUp(distance_calculation_size,
					GetSeedsMemory());
			while (!PopHashPositionDataListIdsForRepresentativeSearch(
					&representative_search_hash_position_data_list_ids)
					|| distance_calculation_seed_list_.GetNumberHashPositionDataListIds()
							> 0) {
				// rep search
				while (CanRunRepresentativeSearch()) {
					uint32_t hash_position_data_list_id =
							representative_search_hash_position_data_list_ids.back();
					representative_search_hash_position_data_list_ids.pop_back();
					shared_parameters->seed_searcher->SearchSeedsFromRepresentativesForSimilarityFiltering(
							hash_position_data_list_id,
							&distance_calculation_seed_list_);
				}

				if (distance_calculation_seed_list_.GetNumberSeeds() > 0) {
					int resource_id = last_run_gpu_controller_resource_id_;
					last_run_gpu_controller_resource_id_ = -1;
					RunDistanceCalculation(gpu_stream_controller,
							*(distance_calculation_seed_list_.MoveHostMemory()));
					gpu_stream_controller.WaitRun(resource_id);
					MoveUngappedExtensionWithTriggerResults(resource_id,
							gpu_stream_controller,
							&temp_chain_filter_seed_list);
					//last_run_gpu_controller_resource_id_ = -1;
					//RunDistanceCalculation(gpu_stream_controller, *(distance_calculation_seed_list_.MoveHostMemory()));
				} else {
					ReleaseSeedsMemory(distance_calculation_seed_list_);
				}
				if (distance_calculation_seed_list_.GetNumberHashPositionDataListIds()
						> 0) {
					// seed search
					size_t number_hash_position_data_list_ids =
							distance_calculation_seed_list_.GetNumberHashPositionDataListIds();
					size_t hash_position_data_list_ids_i = 0;
					ungapped_extension_with_trigger_seed_list_.SetUp(
							ungapped_extension_with_trigger_size,
							GetSeedsMemory());

					while (hash_position_data_list_ids_i
							< number_hash_position_data_list_ids) {
						if (ungapped_extension_with_trigger_seed_list_.Filled()) {
							int resource_id =
									last_run_gpu_controller_resource_id_;
							last_run_gpu_controller_resource_id_ = -1;
							RunUngappedExtensionWithTrigger(
									gpu_stream_controller,
									*(ungapped_extension_with_trigger_seed_list_.MoveHostMemory()));
							gpu_stream_controller.WaitRun(resource_id);
							if (gpu_run_state_list_[resource_id]
									== kDistanceCalculation) {
								MoveDistanceCalculationResults(resource_id,
										gpu_stream_controller,
										&distance_calculation_seed_list_);
								//last_run_gpu_controller_resource_id_ = -1;
							} else {
								MoveUngappedExtensionWithTriggerResults(
										resource_id, gpu_stream_controller,
										&temp_chain_filter_seed_list);
								//last_run_gpu_controller_resource_id_ = -1;
							}
							ungapped_extension_with_trigger_seed_list_.SetUp(
									ungapped_extension_with_trigger_size,
									GetSeedsMemory());
						}
						if (!ungapped_extension_with_trigger_seed_list_.Filled()) {
							shared_parameters->seed_searcher->SearchSeedsFromNonClusteringSubsequence(
									hash_position_data_list_ids_i,
									distance_calculation_seed_list_,
									&ungapped_extension_with_trigger_seed_list_);
							++hash_position_data_list_ids_i;
						}
					}

					if (gpu_run_state_list_[last_run_gpu_controller_resource_id_]
							== kDistanceCalculation) {
						gpu_stream_controller.WaitRun(
								last_run_gpu_controller_resource_id_);
						MoveDistanceCalculationResults(
								last_run_gpu_controller_resource_id_,
								gpu_stream_controller,
								&distance_calculation_seed_list_);
						last_run_gpu_controller_resource_id_ = -1;
					}

					hash_position_data_list_ids_i = 0;
					while (hash_position_data_list_ids_i
							< number_hash_position_data_list_ids) {
						if (ungapped_extension_with_trigger_seed_list_.Filled()) {
							int resource_id =
									last_run_gpu_controller_resource_id_;
							last_run_gpu_controller_resource_id_ = -1;
							RunUngappedExtensionWithTrigger(
									gpu_stream_controller,
									*(ungapped_extension_with_trigger_seed_list_.MoveHostMemory()));
							gpu_stream_controller.WaitRun(resource_id);
							MoveUngappedExtensionWithTriggerResults(resource_id,
									gpu_stream_controller,
									&temp_chain_filter_seed_list);

							ungapped_extension_with_trigger_seed_list_.SetUp(
									ungapped_extension_with_trigger_size,
									GetSeedsMemory());
						}
						if (!ungapped_extension_with_trigger_seed_list_.Filled()) {
							shared_parameters->seed_searcher->SearchSeedsFromClusteringSubsequence(
									hash_position_data_list_ids_i,
									distance_calculation_seed_list_,
									passing_distance_calculation_i_list,
									&ungapped_extension_with_trigger_seed_list_);
							++hash_position_data_list_ids_i;
						}
					}
					while (ungapped_extension_with_trigger_seed_list_.GetSeedsLength()
							> 0) {
						int resource_id = last_run_gpu_controller_resource_id_;
						last_run_gpu_controller_resource_id_ = -1;
						RunUngappedExtensionWithTrigger(gpu_stream_controller,
								*(ungapped_extension_with_trigger_seed_list_.MoveHostMemory()));
						gpu_stream_controller.WaitRun(resource_id);
						MoveUngappedExtensionWithTriggerResults(resource_id,
								gpu_stream_controller,
								&temp_chain_filter_seed_list);
						ungapped_extension_with_trigger_seed_list_.SetUp(
								ungapped_extension_with_trigger_size,
								GetSeedsMemory());
					}
					ReleaseSeedsMemory(
							ungapped_extension_with_trigger_seed_list_);
					distance_calculation_seed_list_.SetUp(
							distance_calculation_size, GetSeedsMemory());
				}
				assert(
						distance_calculation_seed_list_.GetNumberHashPositionDataListIds()
								<= 1);
			}
			gpu_stream_controller.WaitRun(last_run_gpu_controller_resource_id_);
			MoveUngappedExtensionWithTriggerResults(
					last_run_gpu_controller_resource_id_, gpu_stream_controller,
					&temp_chain_filter_seed_list);
			last_run_gpu_controller_resource_id_ = -1;
			// finalize
			ReleaseSeedsMemory(distance_calculation_seed_list_);
			assert(
					distance_calculation_seed_list_.GetNumberHashPositionDataListIds()
							== 0);
			assert(
					distance_calculation_seed_list_.GetNumberHashPositionDataListIds()
							== 0);
			assert(
					ungapped_extension_with_trigger_seed_list_.GetSeedsLength()
							== 0);

			cout << thread_id << " thread finished seed search " << endl;

			shared_parameters->presearch_barrier->wait();
			SetChainFilteringSeeds();
			shared_parameters->presearch_barrier->wait();

			assert(available_seeds_memory_ids_.size() == kNumberOfGpuResources);
			assert(
					gpu_stream_controller.GetNumberOfAvailableResources()
							== kNumberOfGpuResources);
			cout << thread_id << " thread starts extension " << endl;
			size_t next_chain_filter_query_id = UINT_MAX;
			vector<SeedSearcher::Hit> *next_hits = NULL;
			for (size_t extension_seed_lists_i = 0;
					extension_seed_lists_i < kNumberOfExtensionSeedLists;
					++extension_seed_lists_i) {
				ExtensionSeedList &seed_list =
						extension_seed_lists_[extension_seed_lists_i];
				seed_list.SetUp(extension_size, GetSeedsMemory());
			}
			while (next_chain_filter_query_id != UINT_MAX
					|| !PopQueryIdsForChainFilter(&chain_filter_query_ids)) {
				for (size_t extension_seed_lists_i = 0;
						extension_seed_lists_i < kNumberOfExtensionSeedLists;
						++extension_seed_lists_i) {
					ExtensionSeedList &seed_list =
							extension_seed_lists_[extension_seed_lists_i];
					if (seed_list.GetState() == ExtensionSeedList::kRun) {
						gpu_stream_controller.WaitRun(
								last_run_gpu_controller_resource_id_);
						seed_list.FinishRun();
					}
					if (seed_list.GetState() == ExtensionSeedList::kEnd) {
						MoveGappedExtensionResults(
								last_run_gpu_controller_resource_id_,
								gpu_stream_controller, extension_seed_lists_i,
								alignment_start_positions_in_database);
						last_run_gpu_controller_resource_id_ = -1;
					}
					assert(seed_list.GetState() == ExtensionSeedList::kIdle);
					while (next_chain_filter_query_id != UINT_MAX
							|| !PopQueryIdsForChainFilter(
									&chain_filter_query_ids)) {
						while (next_chain_filter_query_id != UINT_MAX
								|| !chain_filter_query_ids.empty()) {
							if (next_chain_filter_query_id == UINT_MAX) {
								next_chain_filter_query_id =
										chain_filter_query_ids.back();
								chain_filter_query_ids.pop_back();
								Query *extended_query = queries.GetQuery(
										next_chain_filter_query_id);

								next_hits =
										&(*shared_parameters->chain_filter_seed_lists)[next_chain_filter_query_id];
#if 0
								for (size_t i = 0; i < next_hits.size(); ++i) {
									assert(0 < next_hits[i].query_sequence_position);
									assert(0 < next_hits[i].database_sequence_position);
									assert(
											next_hits[i].query_sequence_position
											< extended_query->GetSequenceLength());
									assert(
											next_hits[i].database_sequence_position
											< shared_parameters->database_concatenated_sequence_length);
								}
#endif
								chain_filter.Filter(
										extended_query->GetSequence(),
										extended_query->GetSequenceLength(),
										(*shared_parameters->ungapped_extension_cutoffs)[next_chain_filter_query_id],
										shared_parameters->database_concatenated_sequence,
										*next_hits);

								if (!gapped_extension_gpu_flag
										|| seed_list.GetState()
												== ExtensionSeedList::kRun
										|| next_hits->size()
												>= extension_size) {
									tmp_cpu_presearch_results.clear();
									ExtendOnCpu(ungapped_extender,
											gapped_extender,
											next_chain_filter_query_id,
											extended_query,
											shared_parameters->database_concatenated_sequence,
											parameters,
											(*shared_parameters->ungapped_extension_cutoffs),
											thread_parameters_->gapped_extension_cutoff,
											*next_hits,
											tmp_cpu_presearch_results);

									vector<PresearchedResult> &results_list =
											(*shared_parameters->results_list)[next_chain_filter_query_id];

									AddResults((*shared_parameters->database),
											parameters,
											alignment_start_positions_in_database,
											tmp_cpu_presearch_results,
											results_list);
#if 0
									// debug ///////////////////////////////////////////
									cout << "query id = " << extended_query_id << endl;
									for (size_t i = 0; i < results_list.size(); ++i) {
										cout << "q = " << results_list[i].hit.query_position
										<< ", d = "
										<< results_list[i].hit.database_position
										<< " score = " << results_list[i].score << endl;
									}
									////////////////////////////////////////////////////
#endif
									next_hits->clear();
									next_hits = NULL;
									next_chain_filter_query_id = UINT_MAX;
								}
							}
							if (next_chain_filter_query_id != UINT_MAX) {
								if (!next_hits->empty()) {
									if (!seed_list.PushExtensionSeeds(
											next_chain_filter_query_id,
											(*shared_parameters->query_sequence_offsets)[next_chain_filter_query_id],
											*next_hits)) {
										next_hits->clear();
									} else {
										break;
									}
								}
								if (next_hits->empty()) {
									next_hits = NULL;
									next_chain_filter_query_id = UINT_MAX;
								}
							}
						}
						if (next_chain_filter_query_id != UINT_MAX) {
							break;
						}
					}

					if (seed_list.GetSeedsLength() > 0) {
						int resource_id = last_run_gpu_controller_resource_id_;
						last_run_gpu_controller_resource_id_ = -1;

						seed_list.ConvertGappedExtensionSeeds(queries,
								(*shared_parameters->query_sequence_offsets));
						RunGappedExtension(extension_seed_lists_i,
								gpu_stream_controller,
								*(seed_list.GetSeedsMemory()));

						if (resource_id >= 0
								&& gpu_run_state_list_[resource_id] != kIdle) {

							gpu_stream_controller.WaitRun(resource_id);
							size_t prev_id = GetPreviousExtensionSeedListId(
									extension_seed_lists_i);
							assert(
									resource_id == -1
											|| gpu_run_state_list_[resource_id]
													== kGappedExtension);
							MoveGappedExtensionResults(resource_id,
									gpu_stream_controller, prev_id,
									alignment_start_positions_in_database);
						}
					}
				}
			}
			for (size_t extension_seed_lists_i = 0;
					extension_seed_lists_i < kNumberOfExtensionSeedLists;
					++extension_seed_lists_i) {
				ExtensionSeedList &seed_list =
						extension_seed_lists_[extension_seed_lists_i];
				if (seed_list.GetState() != ExtensionSeedList::kIdle) {
					gpu_stream_controller.WaitRun(
							last_run_gpu_controller_resource_id_);
					MoveGappedExtensionResults(
							last_run_gpu_controller_resource_id_,
							gpu_stream_controller, extension_seed_lists_i,
							alignment_start_positions_in_database);
					last_run_gpu_controller_resource_id_ = -1;
				}
				assert(seed_list.GetState() == ExtensionSeedList::kIdle);
				ReleaseSeedsMemory(seed_list);
			}
			cout << thread_id << " thread finished presearch " << endl;
		}
		shared_parameters->all_barrier->wait();
		if (thread_id == 0) {
			shared_parameters->database->NextChunk();
			//cudaDeviceReset();
			//exit(0);
		}
		shared_parameters->all_barrier->wait();
	}

}

size_t AlignerGpuPresearchThread::GetSeedsMemorySize(
		size_t distance_calculation_size,
		size_t ungapped_extension_with_trigger_size, size_t extension_size) {
	size_t max_size = distance_calculation_size;
	max_size = max(max_size, ungapped_extension_with_trigger_size);
	max_size = max(max_size, extension_size * 2);
	return max_size;
}
int AlignerGpuPresearchThread::PopHashPositionDataListIdsForRepresentativeSearch(
		std::vector<uint32_t> *query_ids) {
	uint32_t start_id = 0;
	uint32_t end_id = 0;
	if (query_ids->empty()) {
		boost::unique_lock<boost::mutex> lock(
				thread_parameters_->shared_parameters->next_hash_position_data_list_id_mutex);
		start_id =
				thread_parameters_->shared_parameters->next_hash_position_data_list_id;
		end_id =
				std::min(
						thread_parameters_->shared_parameters->number_hash_position_data_lists,
						start_id + kHashPositionDataListIdsIncrement);
		thread_parameters_->shared_parameters->next_hash_position_data_list_id =
				end_id;
	}
	for (uint32_t i = start_id; i < end_id; ++i) {
		query_ids->push_back(i);
	}
	return query_ids->empty() ? 1 : 0;
}

int AlignerGpuPresearchThread::PopQueryIdsForChainFilter(
		std::vector<uint32_t> *query_ids) {
	uint32_t start_id = 0;
	uint32_t end_id = 0;
	if (query_ids->empty()) {
		boost::unique_lock<boost::mutex> lock(
				thread_parameters_->shared_parameters->next_chain_filter_query_id_mutex);
		start_id =
				thread_parameters_->shared_parameters->next_chain_filter_query_id;
		end_id =
				std::min(
						thread_parameters_->shared_parameters->queries->GetNumberOfSequences(),
						start_id + kQueryIdsIncrement);
		thread_parameters_->shared_parameters->next_chain_filter_query_id =
				end_id;
	}
	for (uint32_t i = start_id; i < end_id; ++i) {
		query_ids->push_back(i);
	}
	return query_ids->empty() ? 1 : 0;
}

void AlignerGpuPresearchThread::AddResults(DatabaseType &database,
		AligningParameters &parameters,
		std::vector<std::vector<AlignmentPosition> > &alignment_start_positions_in_database,
		std::vector<PresearchedResult> &added_results,
		std::vector<PresearchedResult> &results) {
	AlignerCommon::AddResults(database, parameters,
			alignment_start_positions_in_database, added_results, results);
}

void AlignerGpuPresearchThread::UpdateDatabaseData(
		ThreadSharedParameters* shared_parameters) {
	shared_parameters->database_concatenated_sequence =
			shared_parameters->database->GetConcatenatedSequence();
	shared_parameters->database_concatenated_sequence_length =
			shared_parameters->database->GetConcatenatedSequenceLength();
	shared_parameters->seed_searcher->Reset(
			*(shared_parameters->seed_searcher_query_parameters),
			shared_parameters->database->GetSeedSearcherParameters());
}

int AlignerGpuPresearchThread::RunDistanceCalculation(
		GpuStreamController &controller, HostSeedsMemory &host_seeds_memory) {
	GpuStreamController::DistanceCalculationSeedsParameters p;
	p.size = std::min(distance_calculation_seed_list_.GetNumberSeeds(),
			distance_calculation_seed_list_.GetNumberSeedsThreshold());
	last_run_gpu_controller_resource_id_ = controller.Run(p, host_seeds_memory);
	run_host_seeds_memory_ids_[last_run_gpu_controller_resource_id_] =
			host_seeds_memory.GetId();
	last_distance_calculation_seeds_parameters_list_[last_run_gpu_controller_resource_id_] =
			p;
	gpu_run_state_list_[last_run_gpu_controller_resource_id_] =
			kDistanceCalculation;
	return 0;
}

int AlignerGpuPresearchThread::RunUngappedExtensionWithTrigger(
		GpuStreamController &controller, HostSeedsMemory &host_seeds_memory) {
	GpuStreamController::UngappedExtensionWithTriggerSeedsParameters p;
	p.query_seed_size =
			ungapped_extension_with_trigger_seed_list_.GetQuerySeedsLength();
	p.size =
			std::min(
					ungapped_extension_with_trigger_seed_list_.GetSeedsLength(),
					ungapped_extension_with_trigger_seed_list_.GetSeedsLengthThreshold());

	uint32_t* query_seed_starts = host_seeds_memory.GetQuerySeedStarts();
	query_seed_starts[p.query_seed_size] = UINT_MAX;
	last_run_gpu_controller_resource_id_ = controller.Run(p, host_seeds_memory);
	run_host_seeds_memory_ids_[last_run_gpu_controller_resource_id_] =
			host_seeds_memory.GetId();
	last_ungapped_extension_with_trigger_seeds_parameters_list_[last_run_gpu_controller_resource_id_] =
			p;
	gpu_run_state_list_[last_run_gpu_controller_resource_id_] =
			kUngappedExtensionWithTrigger;
	return 0;
}

int AlignerGpuPresearchThread::RunGappedExtension(size_t seed_lists_id,
		GpuStreamController &controller, HostSeedsMemory &host_seeds_memory) {
#if 0
	cout << "before run " << endl;
	cout << "num push q = "
	<< extension_seed_lists_[gapped_extension_seed_lists_ids_].query_ids_.size_()
	<< endl;
	cout << "num run q = "
	<< extension_seed_lists_[gapped_extension_gpu_results_lists_run_id_].query_ids_.size_()
	<< endl;
#endif
	ExtensionSeedList &seed_list = extension_seed_lists_[seed_lists_id];
	GpuStreamController::GappedExtensionSeesParameters p;
	p.forward_length = seed_list.GetForwardGpuSeedsLength();
	p.reverse_length = seed_list.GetReverseGpuSeedsLength();
	p.forward_offset = seed_list.GetForwardGpuSeedsOffset();
	p.reverse_offset = seed_list.GetReverseGpuSeedsOffset();
	last_run_gpu_controller_resource_id_ = controller.Run(p, host_seeds_memory);
	run_host_seeds_memory_ids_[last_run_gpu_controller_resource_id_] =
			host_seeds_memory.GetId();
	gpu_run_state_list_[last_run_gpu_controller_resource_id_] =
			kGappedExtension;
	seed_list.StartRun();
	return 0;
}

int AlignerGpuPresearchThread::MoveDistanceCalculationResults(
		int gpu_controller_resource_id, GpuStreamController &controller,
		DistanceCalculationSeedList *distance_calculation_seed_list) {
	if (gpu_run_state_list_[gpu_controller_resource_id]
			== kDistanceCalculation) {
		HostSeedsMemory &host_seeds_memory =
				host_seeds_memories_[run_host_seeds_memory_ids_[gpu_controller_resource_id]];
		size_t number_seeds =
				last_distance_calculation_seeds_parameters_list_[gpu_controller_resource_id].size;
		const DistanceCalculatorGpu::Distance* result_values =
				(DistanceCalculatorGpu::Distance*) host_seeds_memory.GetResultValues();
		memcpy(distance_calculation_seed_list->GetDistances(), result_values,
				sizeof(result_values[0]) * number_seeds);
		//std::copy(result_values, result_values + host_buffer_size, distance_calculation_seed_list->GetDistances());
#if 0
		uint32_t* query_concatenated_positions =
		host_buffer.GetQueryConcatenatedPositions();
		uint32_t* database_positions = host_buffer.GetDatabasePositions();
		assert(
				distance_calculation_seed_list->GetSeedsLength()
				== host_buffer.GetSize());
		for (uint32_t i = 0; i < number_seeds; ++i) {
			cout << i <<": q "
			<< distance_calculation_seed_list->GetQueryConcatenatedPositions()[i]
			<< ", d "
			<< distance_calculation_seed_list->GetDatabasePositions()[i]
			<< ", dis " << result_values[i] << endl;
			assert(
					distance_calculation_seed_list->GetQueryConcatenatedPositions()[i]
					== query_concatenated_positions[i]);
			assert(
					distance_calculation_seed_list->GetDatabasePositions()[i]
					== database_positions[i]);
		}
#endif
		ReleaseSeedsMemory(host_seeds_memory);
		last_distance_calculation_seeds_parameters_list_[gpu_controller_resource_id].size =
				0;
		gpu_run_state_list_[gpu_controller_resource_id] = kIdle;
		controller.ReleaseResource(gpu_controller_resource_id);
	}
	return 0;
}
int AlignerGpuPresearchThread::MoveUngappedExtensionWithTriggerResults(
		int gpu_controller_resource_id, GpuStreamController &controller,
		TempChainFilterSeedList *ungapped_extension_finished_seeds_list) {
	if (gpu_run_state_list_[gpu_controller_resource_id]
			== kUngappedExtensionWithTrigger) {
		HostSeedsMemory &host_seeds_memory =
				host_seeds_memories_[run_host_seeds_memory_ids_[gpu_controller_resource_id]];
		size_t number_seeds =
				last_ungapped_extension_with_trigger_seeds_parameters_list_[gpu_controller_resource_id].size;
		const uint32_t* query_seed_query_ids =
				thread_parameters_->shared_parameters->seed_searcher_query_parameters->GetSeedQueryIds();
		const uint32_t* query_seed_query_positions =
				thread_parameters_->shared_parameters->seed_searcher_query_parameters->GetSeedQueryPositions();
		const uint32_t* query_seed_ids = host_seeds_memory.GetQuerySeedIds();
		const uint32_t* query_seed_starts =
				host_seeds_memory.GetQuerySeedStarts();
		const uint32_t* database_positions =
				host_seeds_memory.GetDatabasePositions();
		const char* next_step_flags =
				(char*) host_seeds_memory.GetResultValues();
		temp_ungapped_extension_seed_data_.clear();
		uint32_t next_query_seed_starts_i = 1;
		for (uint32_t i = 0; i < number_seeds; ++i) {
			if (next_step_flags[i] > 0) {
				while (query_seed_starts[next_query_seed_starts_i] <= i) {
					++next_query_seed_starts_i;
				}
				temp_ungapped_extension_seed_data_.push_back(
						make_pair(query_seed_ids[next_query_seed_starts_i - 1],
								i));
			}
#if 0
			// debug ///////////////////
			uint32_t x_query_position = 18;
			uint32_t x_db_position = 623672280;
			if (x_db_position == database_positions[i]) {
				cout << "q " << query_ids[i] << ", "
				<< query_concatenated_positions[i] << " d "
				<< database_positions[i] << " score "
				<< next_step_flags[i] << endl;
			}
			cout << query_ids[i] << ", " << query_concatenated_positions[i]
			<< " : " << database_positions[i] << endl;
			////////////////////////////
#endif
		}

		size_t buckets_size =
				ungapped_extension_finished_seeds_list->buckets.size();
		for (uint32_t i = 0; i < temp_ungapped_extension_seed_data_.size();
				++i) {
			uint32_t query_seed_i = temp_ungapped_extension_seed_data_[i].first;
			uint32_t seed_i = temp_ungapped_extension_seed_data_[i].second;
			const uint32_t query_id = query_seed_query_ids[query_seed_i];
			uint32_t bucket_i = query_id % buckets_size;
			TempChainFilterSeedBucket &bucket =
					ungapped_extension_finished_seeds_list->buckets[bucket_i];
			bucket.query_ids.push_back(query_id);
			bucket.query_concatenated_positions.push_back(
					query_seed_query_positions[query_seed_i]);
			bucket.database_positions.push_back(database_positions[seed_i]);
#if 0
			cout << query_ids[i] << ", " << query_concatenated_positions[i]
			<< " : " << database_positions[i] << endl;
#endif

		}
		ReleaseSeedsMemory(host_seeds_memory);
		last_ungapped_extension_with_trigger_seeds_parameters_list_[gpu_controller_resource_id].size =
				0;
		gpu_run_state_list_[gpu_controller_resource_id] = kIdle;
		controller.ReleaseResource(gpu_controller_resource_id);
	}
	return 0;
}

int AlignerGpuPresearchThread::MoveGappedExtensionResults(
		int gpu_controller_resource_id, GpuStreamController &controller,
		size_t seed_list_id,
		std::vector<std::vector<AlignmentPosition> > &alignment_start_positions_in_database) {
	ExtensionSeedList &seed_list = extension_seed_lists_[seed_list_id];

#if 0
	cout << "move resutls " << endl;
	cout << "num push q = "
	<< extension_seed_lists_[gapped_extension_seed_lists_ids_].query_ids_.size_()
	<< endl;
	cout << "num run q = "
	<< extension_seed_lists_[gapped_extension_gpu_results_lists_run_id_].query_ids_.size_()
	<< endl;
#endif

	std::vector<uint32_t> &query_offsets =
			(*thread_parameters_->shared_parameters->query_sequence_offsets);
	if (seed_list.GetSeedsLength() != 0) {
		ExtensionSeedList::SeedData *seed_data_list =
				seed_list.GetSeedDataList();
		HostSeedsMemory *host_memory = seed_list.GetSeedsMemory();
		uint32_t* query_concatenated_positions =
				host_memory->GetQueryConcatenatedPositions();
		uint32_t* database_positions = host_memory->GetDatabasePositions();
		int* result_values = host_memory->GetResultValues();
		uint32_t query_id = seed_data_list[0].query_id;
		uint32_t query_offset = query_offsets[query_id];
		size_t number_seeds = seed_list.GetSeedsLength();
		vector<AlignerCommon::PresearchedResult> temp_results;
		AlignerCommon::PresearchedResult temp_result;
		temp_result.database_chunk_id = UINT_MAX;
		temp_result.subject_id_in_chunk = UINT_MAX;
		temp_result.hit_count = 1;

		DatabaseType &database =
				*(thread_parameters_->shared_parameters->database);
		AligningParameters &parameters =
				*(thread_parameters_->shared_parameters->parameters);
		for (size_t seed_i = 0; seed_i < number_seeds; ++seed_i) {
			ExtensionSeedList::SeedData &seed_data = seed_data_list[seed_i];
			if (query_id != seed_data.query_id) {
				std::vector<PresearchedResult> &results =
						(*thread_parameters_->shared_parameters->results_list)[query_id];
#if 0
				// debug ///////////////////////////////////////////
				cout << "query id = " << query_id << endl;
				for (size_t i = 0; i < results.size(); ++i) {
					cout << "q = " << results[i].hit.query_position << ", d = "
					<< results[i].hit.database_position << " score = "
					<< results[i].score << endl;
				}
				////////////////////////////////////////////////////
#endif
				AddResults(database, parameters,
						alignment_start_positions_in_database, temp_results,
						results);
				query_id = seed_data.query_id;
				query_offset = query_offsets[query_id];
				temp_results.clear();
			}
			temp_result.score = 0;
			temp_result.hit = seed_data.hit;

			temp_result.start.query_position =
					query_concatenated_positions[seed_data.reverse_seed_id]
							- query_offset;
			temp_result.start.database_position =
					database_positions[seed_data.reverse_seed_id];
			temp_result.score += result_values[seed_data.reverse_seed_id];

			temp_result.end.query_position =
					query_concatenated_positions[seed_data.foward_seed_id]
							- query_offset;
			temp_result.end.database_position =
					database_positions[seed_data.foward_seed_id];
			temp_result.score += result_values[seed_data.foward_seed_id];

			temp_results.push_back(temp_result);

#if 0
			if (query_id == 10) {
				cout << "reverse gapped score : "
				<< result_values[seed_data.reverse_seed_id] << endl;
				cout << "forward gapped score : "
				<< result_values[seed_data.foward_seed_id] << endl;
				cout << "gapped extension gpu" << endl;
				cout << "seed_i " << seed_i << endl;
				cout << "s " << temp_result.score << endl;
				cout << "q " << temp_result.start.query_position << " - "
				<< temp_result.end.query_position << "( "
				<< temp_result.hit.query_position << ")" << endl;
				cout << "d " << temp_result.start.database_position << " - "
				<< temp_result.end.database_position << "( "
				<< temp_result.hit.database_position << ")" << endl;
			}
#endif
		}
		if (!temp_results.empty()) {
			std::vector<PresearchedResult> &results =
					(*thread_parameters_->shared_parameters->results_list)[query_id];
			AddResults(database, parameters,
					alignment_start_positions_in_database, temp_results,
					results);
		}
	}
	seed_list.Clear();
	gpu_run_state_list_[gpu_controller_resource_id] = kIdle;
	controller.ReleaseResource(gpu_controller_resource_id);
	return 0;
}

int AlignerGpuPresearchThread::FinalizeUngappedExtension(
		int gpu_controller_resource_id, GpuStreamController &controller,
		size_t seed_list_id) {
	ExtensionSeedList &seed_list = extension_seed_lists_[seed_list_id];
	seed_list.FinishRun();
	gpu_run_state_list_[gpu_controller_resource_id] = kIdle;
	controller.ReleaseResource(gpu_controller_resource_id);
	return 0;
}

int AlignerGpuPresearchThread::ExtendOnCpu(UngappedExtender &ungapped_extender,
		GappedExtender &gapped_extender, uint32_t query_id, Query *query,
		const AlphabetCoder::Code *database_concatenated_sequence,
		AligningParameters &parameters,
		std::vector<int> &ungapped_extension_cutoffs,
		int gapped_extension_cutoff, std::vector<SeedSearcherGpu::Hit> &hits,
		std::vector<PresearchedResult> &results) {
	PresearchedResult tmp_result;
	for (size_t hit_id = 0; hit_id < hits.size(); ++hit_id) {
#if DEBUG
		++total_gapped_extention_count;
#endif
		int score = 0;
		int query_position;
		int database_position;
		tmp_result.score = 0;
		tmp_result.hit.query_position = hits[hit_id].query_sequence_position;
		tmp_result.hit.database_position =
				hits[hit_id].database_sequence_position;

		gapped_extender.ExtendOneSide(
				query->GetSequence() + tmp_result.hit.query_position - 1,
				tmp_result.hit.query_position - 1,

				database_concatenated_sequence
						+ tmp_result.hit.database_position - 1,
				query->GetSequenceDelimiter(), true, parameters.score_matrix,
				parameters.gap_open, parameters.gap_extension,
				gapped_extension_cutoff, &score, &query_position,
				&database_position, NULL);

		tmp_result.score += score;
		tmp_result.start.query_position = tmp_result.hit.query_position - 1
				+ query_position;
		tmp_result.start.database_position = tmp_result.hit.database_position
				- 1 + database_position;

		gapped_extender.ExtendOneSide(
				query->GetSequence() + tmp_result.hit.query_position,
				query->GetSequenceLength() - tmp_result.hit.query_position,
				database_concatenated_sequence
						+ tmp_result.hit.database_position,
				query->GetSequenceDelimiter(), false, parameters.score_matrix,
				parameters.gap_open, parameters.gap_extension,
				gapped_extension_cutoff, &score, &query_position,
				&database_position, NULL);
		tmp_result.score += score;
		tmp_result.end.query_position += query_position + 1;
		tmp_result.end.database_position += database_position + 1;

		tmp_result.hit_count = hits[hit_id].k_mer_count;
		results.push_back(tmp_result);
#if 0
		if (query_id == 10) {
			cout << "gapped extension cpu" << endl;
			cout << "hit id " << hit_id << endl;
			cout << "s " << tmp_result.score << endl;
			cout << "q " << tmp_result.start.query_position << " - "
			<< tmp_result.end.query_position << "( "
			<< tmp_result.hit.query_position << ")" << endl;
			cout << "d " << tmp_result.start.database_position << " - "
			<< tmp_result.end.database_position << "( "
			<< tmp_result.hit.database_position << ")" << endl;
		}
#endif
	}
	return 0;
}

int AlignerGpuPresearchThread::GappedExtendBackSeedsOnCpu(size_t seed_list_id,
		size_t default_extension_size, GappedExtender &gapped_extender,
		Queries &queries,
		const AlphabetCoder::Code *queries_concatenated_sequence,
		const AlphabetCoder::Code *database_concatenated_sequence,
		AligningParameters &parameters, int gapped_extension_cutoff) {

	AlphabetCoder::Code sequence_delimiter =
			queries.GetQuery(0)->GetSequenceDelimiter();
	ExtensionSeedList &seed_list = extension_seed_lists_[seed_list_id];
	HostSeedsMemory *host_memory = seed_list.GetSeedsMemory();
	uint32_t* query_concatenated_positions =
			host_memory->GetQueryConcatenatedPositions();
	uint32_t* database_positions = host_memory->GetDatabasePositions();
	int* result_values = host_memory->GetResultValues();
	size_t max_query_length = queries.GetMaxQuerySequenceLength();
	size_t reverse_offset = seed_list.GetReverseGpuSeedsOffset();
	size_t reverse_end = reverse_offset + seed_list.GetReverseGpuSeedsLength();
	size_t reverse_start = reverse_offset;
	if (reverse_end > (reverse_offset + default_extension_size)) {
		reverse_start = reverse_end - default_extension_size;
	}
	for (size_t seed_i = reverse_start; seed_i < reverse_end; ++seed_i) {
#if DEBUG
		++total_gapped_extention_count;
#endif
		int query_position;
		int database_position;
#if 0
		size_t q_len = 0;
		for (int i = query_concatenated_positions[seed_i]; i > 0 && queries_concatenated_sequence[i] != sequence_delimiter; --i) {
			++q_len;
		}
		assert(q_len < max_query_length);
#endif
		gapped_extender.ExtendOneSide(
				queries_concatenated_sequence
						+ query_concatenated_positions[seed_i],
				max_query_length,
				database_concatenated_sequence + database_positions[seed_i],
				sequence_delimiter, true, parameters.score_matrix,
				parameters.gap_open, parameters.gap_extension,
				gapped_extension_cutoff, &result_values[seed_i],
				&query_position, &database_position, NULL);

		query_concatenated_positions[seed_i] += query_position;
		database_positions[seed_i] += database_position;
	}

	size_t foward_offset = seed_list.GetForwardGpuSeedsOffset();
	size_t foward_end = foward_offset + seed_list.GetForwardGpuSeedsLength();
	size_t foward_start = foward_offset;
	if (foward_end > (foward_offset + default_extension_size)) {
		foward_start = foward_end - default_extension_size;
	}
	for (size_t seed_i = foward_start; seed_i < foward_end; ++seed_i) {
#if DEBUG
		++total_gapped_extention_count;
#endif
		int query_position;
		int database_position;
		gapped_extender.ExtendOneSide(
				queries_concatenated_sequence
						+ query_concatenated_positions[seed_i],
				max_query_length,
				database_concatenated_sequence + database_positions[seed_i],
				sequence_delimiter, false, parameters.score_matrix,
				parameters.gap_open, parameters.gap_extension,
				gapped_extension_cutoff, &result_values[seed_i],
				&query_position, &database_position, NULL);

		query_concatenated_positions[seed_i] += query_position;
		database_positions[seed_i] += database_position;

	}
	size_t extension_size = foward_end - foward_start;
	if (extension_size > 0) {
		seed_list.SetComputedFlagForBackSeeds(extension_size);
		return 0;
	} else {
		return 1;
	}
}

void AlignerGpuPresearchThread::SetChainFilteringSeeds() {
	int number_threads =
			thread_parameters_->shared_parameters->parameters->number_threads;
	int thread_id = thread_parameters_->thread_id;
	std::vector<AlignerGpuPresearchThread::TempChainFilterSeedList> &temp_chain_filter_seed_lists =
			*(thread_parameters_->shared_parameters->temp_chain_filter_seed_lists);
	std::vector<std::vector<SeedSearcherGpu::Hit> > &chain_filter_seed_lists =
			*(thread_parameters_->shared_parameters->chain_filter_seed_lists);
	std::vector<uint32_t> &query_sequence_offsets =
			*(thread_parameters_->shared_parameters->query_sequence_offsets);
	for (int i = 0; i < number_threads; ++i) {
		TempChainFilterSeedBucket &seed_bucket =
				temp_chain_filter_seed_lists[i].buckets[thread_id];
		size_t number_seeds = seed_bucket.query_ids.size();
		for (size_t seed_i = 0; seed_i < number_seeds; ++seed_i) {
			SeedSearcherGpu::Hit h;
			h.k_mer_count = 1;
			uint32_t query_id = seed_bucket.query_ids[seed_i];
			h.query_sequence_position =
					seed_bucket.query_concatenated_positions[seed_i]
							- query_sequence_offsets[query_id];
			h.database_sequence_position =
					seed_bucket.database_positions[seed_i];
			chain_filter_seed_lists[query_id].push_back(h);
		}
		seed_bucket.query_ids.clear();
		seed_bucket.query_concatenated_positions.clear();
		seed_bucket.database_positions.clear();
	}
}
