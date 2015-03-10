/*
 * aligner_presearch_thread.cpp
 *
 *  Created on: Sep 19, 2014
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
#include "aligner_presearch_thread.h"

#include <boost/thread.hpp>

using namespace std;

const size_t AlignerPresearchThread::kMaxThreadExtensionSeedListSize = 1 << 20;

AlignerPresearchThread::AlignerPresearchThread() :
		thread_parameters_(NULL) {
}

AlignerPresearchThread::~AlignerPresearchThread() {
// TODO Auto-generated destructor stub
}

bool AlignerPresearchThread::CanRunSeedSearch() {
	return !PopHashPositionDataListIds(&hash_position_data_list_ids);
}

void AlignerPresearchThread::Run(ThreadParameters& thread_parameters) {
	thread_parameters_ = &thread_parameters;
	int thread_id = thread_parameters_->thread_id;
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
		SeedSearcher::QueryParameters::BuildParameters seed_search_query_parameters_build_parameter;
		seed_search_query_parameters_build_parameter.score_matrix =
				parameters.score_matrix;
		seed_search_query_parameters_build_parameter.ungapped_extension_cutoffs_ptr =
				shared_parameters->ungapped_extension_cutoffs;
		seed_search_query_parameters_build_parameter.gapped_extension_triggers_ptr =
				shared_parameters->gapped_extension_triggers;
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

	}

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

	std::vector<uint32_t> hitting_query_position_data_i_list;
	UngappedExtender ungapped_extender;
	bool representative_search_finished = false;
	bool seed_search_finished = false;

	shared_parameters->all_barrier->wait();
	while (shared_parameters->database->GetChunkId()
			< (int) shared_parameters->database->GetNumberChunks()) {
		representative_search_finished = false;
		seed_search_finished = false;

		for (size_t i = 0;
				i
						< thread_parameters.shared_parameters->parameters->number_threads;
				++i) {
			TempChainFilterSeedBucket &seed_bucket =
					temp_chain_filter_seed_list.buckets[i];
			assert(seed_bucket.query_ids.empty());
			assert(seed_bucket.query_positions.empty());
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
			if (alignment_start_positions_in_database.size()
					< shared_parameters->database->GetNumberSequencesInChunk()) {
				alignment_start_positions_in_database.resize(
						shared_parameters->database->GetNumberSequencesInChunk());
			}

			while (!PopHashPositionDataListIds(&hash_position_data_list_ids)) {
				while (CanRunSeedSearch()) {
					uint32_t hash_position_data_list_id =
							hash_position_data_list_ids.back();
					hash_position_data_list_ids.pop_back();
					shared_parameters->seed_searcher->Search(
							hash_position_data_list_id,
							hitting_query_position_data_i_list,
							&temp_chain_filter_seed_list);
				}
			}

			cout << thread_id << " thread finished seed search " << endl;

			shared_parameters->presearch_barrier->wait();
			SetChainFilteringSeeds();
			shared_parameters->presearch_barrier->wait();

			cout << thread_id << " thread starts extension " << endl;

			while (!PopQueryIdsForChainFilter(&chain_filter_query_ids)) {
				size_t next_chain_filter_query_id =
						chain_filter_query_ids.back();
				chain_filter_query_ids.pop_back();
				Query *extended_query = queries.GetQuery(
						next_chain_filter_query_id);

				vector<SeedSearcher::Hit>* next_hits =
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
				chain_filter.Filter(extended_query->GetSequence(),
						extended_query->GetSequenceLength(),
						(*shared_parameters->ungapped_extension_cutoffs)[next_chain_filter_query_id],
						shared_parameters->database_concatenated_sequence,
						*next_hits);

				tmp_cpu_presearch_results.clear();
				Extend(ungapped_extender, gapped_extender,
						next_chain_filter_query_id, extended_query,
						shared_parameters->database_concatenated_sequence,
						parameters,
						(*shared_parameters->ungapped_extension_cutoffs),
						thread_parameters_->gapped_extension_cutoff, *next_hits,
						tmp_cpu_presearch_results);

				vector<PresearchedResult> &results_list =
						(*shared_parameters->results_list)[next_chain_filter_query_id];

				AddResults((*shared_parameters->database), parameters,
						alignment_start_positions_in_database,
						tmp_cpu_presearch_results, results_list);
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

size_t AlignerPresearchThread::GetSeedsMemorySize(
		size_t distance_calculation_size,
		size_t ungapped_extension_with_trigger_size, size_t extension_size) {
	size_t max_size = distance_calculation_size;
	max_size = max(max_size, ungapped_extension_with_trigger_size);
	max_size = max(max_size, extension_size * 2);
	return max_size;
}
int AlignerPresearchThread::PopHashPositionDataListIds(
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

int AlignerPresearchThread::PopQueryIdsForChainFilter(
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

void AlignerPresearchThread::AddResults(DatabaseType &database,
		AligningParameters &parameters,
		std::vector<std::vector<AlignmentPosition> > &alignment_start_positions_in_database,
		std::vector<PresearchedResult> &added_results,
		std::vector<PresearchedResult> &results) {
	AlignerCommon::AddResults(database, parameters,
			alignment_start_positions_in_database, added_results, results);
}

void AlignerPresearchThread::UpdateDatabaseData(
		ThreadSharedParameters* shared_parameters) {
	shared_parameters->database_concatenated_sequence =
			shared_parameters->database->GetConcatenatedSequence();
	shared_parameters->database_concatenated_sequence_length =
			shared_parameters->database->GetConcatenatedSequenceLength();
	shared_parameters->seed_searcher->Reset(
			*(shared_parameters->seed_searcher_query_parameters),
			shared_parameters->database->GetSeedSearcherParameters());
}

int AlignerPresearchThread::Extend(UngappedExtender &ungapped_extender,
		GappedExtender &gapped_extender, uint32_t query_id, Query *query,
		const AlphabetCoder::Code *database_concatenated_sequence,
		AligningParameters &parameters,
		std::vector<int> &ungapped_extension_cutoffs,
		int gapped_extension_cutoff, std::vector<SeedSearcher::Hit> &hits,
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

void AlignerPresearchThread::SetChainFilteringSeeds() {
	int number_threads =
			thread_parameters_->shared_parameters->parameters->number_threads;
	int thread_id = thread_parameters_->thread_id;
	std::vector<AlignerPresearchThread::TempChainFilterSeedList> &temp_chain_filter_seed_lists =
			*(thread_parameters_->shared_parameters->temp_chain_filter_seed_lists);
	std::vector<std::vector<SeedSearcher::Hit> > &chain_filter_seed_lists =
			*(thread_parameters_->shared_parameters->chain_filter_seed_lists);
	for (int i = 0; i < number_threads; ++i) {
		TempChainFilterSeedBucket &seed_bucket =
				temp_chain_filter_seed_lists[i].buckets[thread_id];
		size_t number_seeds = seed_bucket.query_ids.size();
		for (size_t seed_i = 0; seed_i < number_seeds; ++seed_i) {
			SeedSearcher::Hit h;
			h.k_mer_count = 1;
			uint32_t query_id = seed_bucket.query_ids[seed_i];
			h.query_sequence_position = seed_bucket.query_positions[seed_i];
			h.database_sequence_position =
					seed_bucket.database_positions[seed_i];
			chain_filter_seed_lists[query_id].push_back(h);
		}
		seed_bucket.query_ids.clear();
		seed_bucket.query_positions.clear();
		seed_bucket.database_positions.clear();
	}
}
