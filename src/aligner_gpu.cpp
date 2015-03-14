/*
 * aligner_gpu.cpp
 *
 *  Created on: Jul 10, 2014
 *      Author: shu
 */

#include "aligner_gpu.h"
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
#include "seed_searcher_gpu.h"
#include "seed_searcher_query_parameters.h"
#include "reduced_alphabet_coder.h"
#include "aligner_gpu_data.h"
#include "cuda_common.h"
#include "gpu_stream_controller.h"
#include "aligner_gpu_presearch_thread.h"
#include "aligner_build_results_thread.h"

#include <boost/thread.hpp>

using namespace std;

AlignerGpu::AlignerGpu() {

}

AlignerGpu::~AlignerGpu() {

}

void AlignerGpu::Align(string &queries_filename, string &database_filename,
		string &output_filename, AligningParameters &parameters) {
	Logger *logger = Logger::GetInstance();
	Queries::Parameters queries_parameters;
	ifstream queries_is(queries_filename.c_str());
	AlignerCommon::BuildQueriesParameters(parameters, queries_parameters);
	DatabaseType database(database_filename);
	ofstream os(output_filename.c_str());
	stringstream ss;
	for (Queries queries(queries_is, queries_parameters);
			queries.GetNumberOfSequences() != 0; queries.Next()) {
		ss.str("");
		ss << "number queries is " << queries.GetNumberOfSequences() << endl;
		logger->Log(ss.str());
		vector<vector<Aligner::PresearchedResult> > presearch_results_list(
				queries.GetNumberOfSequences());
		vector<vector<Aligner::Result> > results_list(
				queries.GetNumberOfSequences());
		ss.str("");
		ss << "start presearch " << endl;
		logger->Log(ss.str());
		Presearch(queries, database, parameters, presearch_results_list);
		ss.str("");
		ss << "start build results" << endl;
		logger->Log(ss.str());
		BuildResults(queries, database, parameters, presearch_results_list,
				results_list);
		ss.str("");
		cout << "write results" << endl;
		logger->Log(ss.str());
		AlignerCommon::WriteOutput(os, queries, database, parameters, results_list);
	}
	queries_is.close();
	os.close();
}

void AlignerGpu::Presearch(Queries &queries, DatabaseType &database,
		AligningParameters &parameters,
		std::vector<std::vector<PresearchedResult> > &results_list) {
	Logger *logger = Logger::GetInstance();
	int device_count = 0;
	cudaGetDeviceCount(&device_count);
	stringstream ss;
	ss << device_count << " GPUs are available." << endl;
	logger->Log(ss.str());
	if (parameters.number_gpus == -1 || parameters.number_gpus > device_count) {
		parameters.number_gpus = device_count;
	}
	stringstream ss;
	ss << "use " << parameters.number_gpus << " GPUs." << endl;
	logger->Log(ss.str());
	vector<int> gpu_ids;
	for (int device_i = 0; device_i < parameters.number_gpus; ++device_i) {
		gpu_ids.push_back(device_i);
	}
	bool database_preload = true;
	Statistics statistics(*(parameters.aligning_sequence_type_ptr));
	Statistics::KarlinParameters gapped_karlin_parameters;
	statistics.CalculateGappedKarlinParameters(parameters.score_matrix,
			parameters.gap_open, parameters.gap_extension,
			&gapped_karlin_parameters);
	int gapped_extension_cutoff = Statistics::Normalized2Nominal(
			parameters.normalized_presearched_gapped_extension_cutoff,
			gapped_karlin_parameters);
	SeedSearcherGpu seed_searcher;
	vector<AlignerGpuPresearchThread::TempChainFilterSeedList> temp_chain_filter_seed_lists(
			parameters.number_threads);
	vector<vector<SeedSearcherGpu::Hit> > chain_filter_seed_lists(
			queries.GetNumberOfSequences());
	vector<vector<std::vector<AlignmentPosition> > > alignment_start_positions_in_database_list(
			parameters.number_threads);
	std::vector<AlphabetCoder::Code> queries_concatenated_sequence;
	std::vector<uint32_t> query_sequence_offsets;
	std::vector<int> ungapped_extension_cutoffs(queries.GetNumberOfSequences());
	std::vector<int> gapped_extension_triggers(queries.GetNumberOfSequences());

	uint32_t total_length_query_sequences = 0;
	for (uint32_t query_id = 0; query_id < queries.GetNumberOfSequences();
			++query_id) {
		Query *query = queries.GetQuery(query_id);
		uint32_t query_sequence_length = query->GetSequenceLength();
		total_length_query_sequences += query_sequence_length;
	}

	queries_concatenated_sequence.resize(
			total_length_query_sequences - queries.GetNumberOfSequences() + 1);
	query_sequence_offsets.resize(queries.GetNumberOfSequences());
	SetQueriesData(queries, database.GetSequenceDelimiter(), parameters,
			queries_concatenated_sequence, query_sequence_offsets,
			ungapped_extension_cutoffs, gapped_extension_triggers);
	vector<AlignerGpuData> gpu_data_list(gpu_ids.size());
	for (size_t i = 0; i < gpu_ids.size(); ++i) {
		int gpu_id = gpu_ids[i];
		AlignerGpuData &gpu_data = gpu_data_list[i];
		gpu_data.SetGpuId(gpu_id);
	}
	SeedSearcherGpu::QueryParameters seed_searcher_query_parameter;
#if DEBUG
	uint32_t total_gapped_extention_count = 0;
#endif
	database.ResetChunk();
	boost::thread_group threads;
	size_t number_presearch_threads = parameters.number_threads;
	size_t number_all_threads = number_presearch_threads;
	if (database_preload) {
		++number_all_threads;
	}
	boost::barrier presearch_barrier(number_presearch_threads);
	boost::barrier all_barrier(number_all_threads);
	AlignerGpuPresearchThread::ThreadSharedParameters shared_parameters;
	shared_parameters.next_hash_position_data_list_id = 0;
	shared_parameters.next_chain_filter_query_id = 0;
	shared_parameters.queries_concatenated_sequence_length =
			queries_concatenated_sequence.size();
	shared_parameters.queries = &queries;
	shared_parameters.queries_concatenated_sequence =
			&queries_concatenated_sequence[0];
	shared_parameters.query_sequence_offsets = &query_sequence_offsets;
	shared_parameters.ungapped_extension_cutoffs = &ungapped_extension_cutoffs;
	shared_parameters.gapped_extension_triggers = &gapped_extension_triggers;

	shared_parameters.database = &database;
	shared_parameters.database_concatenated_sequence = NULL;
	shared_parameters.database_concatenated_sequence_length = 0;
	shared_parameters.parameters = &parameters;
	shared_parameters.seed_searcher = &seed_searcher;
	shared_parameters.seed_searcher_query_parameters =
			&seed_searcher_query_parameter;
	shared_parameters.temp_chain_filter_seed_lists =
			&temp_chain_filter_seed_lists;
	shared_parameters.chain_filter_seed_lists = &chain_filter_seed_lists;
	shared_parameters.gpu_data_list = &gpu_data_list;
	shared_parameters.results_list = &results_list;
	shared_parameters.presearch_barrier = &presearch_barrier;
	shared_parameters.all_barrier = &all_barrier;

	size_t gpu_stream_controller_buffer_size =
			kMaxSeedBufferSizePerGPU
					/ ((number_presearch_threads + gpu_ids.size() - 1)
							/ gpu_ids.size());
	for (size_t i = 0; i < number_presearch_threads; ++i) {
		AlignerGpuPresearchThread t;
		AlignerGpuPresearchThread::ThreadParameters p;
		p.thread_id = i;
		p.gpu_data_list_id = i % gpu_data_list.size();
		p.gpu_id = gpu_ids[p.gpu_data_list_id];
		p.gapped_extension_cutoff = gapped_extension_cutoff;
		p.gpu_seeds_memory_size = gpu_stream_controller_buffer_size;
		p.shared_parameters = &shared_parameters;
		threads.create_thread(
				boost::bind(&AlignerGpuPresearchThread::Run, t, p));
	}
	if (database_preload) {
		AlignerGpuPresearchThread t;
		AlignerGpuPresearchThread::ThreadParameters p;
		p.thread_id = -1;
		p.gpu_data_list_id = 0;
		p.gpu_id = gpu_ids[p.gpu_data_list_id];
		p.gapped_extension_cutoff = 0;
		p.gpu_seeds_memory_size = 0;
		p.shared_parameters = &shared_parameters;
		threads.create_thread(
				boost::bind(&AlignerGpuPresearchThread::Run, t, p));
	}
	threads.join_all();
	return;
}

void AlignerGpu::SetQueriesData(Queries &queries,
		AlphabetCoder::Code sequence_delimiter, AligningParameters &parameters,
		std::vector<AlphabetCoder::Code> &queries_concatenated_sequence,
		std::vector<uint32_t> &query_sequence_offsets,
		std::vector<int> &ungapped_extension_cutoffs,
		std::vector<int> &gapped_extension_triggers) {
	uint32_t offset = 0;
	for (uint32_t query_id = 0; query_id < queries.GetNumberOfSequences();
			++query_id) {
		Query *query = queries.GetQuery(query_id);
		const AlphabetCoder::Code *query_sequence = query->GetSequence();
		uint32_t query_sequence_length = query->GetSequenceLength();
		query_sequence_offsets[query_id] = offset;

		memcpy(&queries_concatenated_sequence[offset], query_sequence,
				sizeof(AlphabetCoder::Code) * (query_sequence_length - 1));

		offset += query_sequence_length - 1;

		ungapped_extension_cutoffs[query_id] =
				Statistics::NormalizedCutoff2NominalCutoff(
						parameters.normalized_presearched_ungapped_extension_cutoff,
						query->GetUngappedKarlinParameters());
		gapped_extension_triggers[query_id] = Statistics::Normalized2Nominal(
				parameters.normalized_presearched_gapped_extension_trigger,
				query->GetUngappedKarlinParameters());

	}
	queries_concatenated_sequence[offset] = sequence_delimiter;
}

void AlignerGpu::BuildResults(Queries &queries, DatabaseType &database,
		AligningParameters &parameters,
		std::vector<std::vector<PresearchedResult> > &presearch_results_list,
		std::vector<std::vector<Result> > &results_list) {
	bool database_preload = true;
	Statistics statistics(*(parameters.aligning_sequence_type_ptr));
	Statistics::KarlinParameters gapped_karlin_parameters;
	statistics.CalculateGappedKarlinParameters(parameters.score_matrix,
			parameters.gap_open, parameters.gap_extension,
			&gapped_karlin_parameters);
	const int gapped_extension_cutoff = Statistics::Normalized2Nominal(
			parameters.normalized_result_gapped_extension_cutoff,
			gapped_karlin_parameters);

	vector<vector<pair<uint32_t, uint32_t> > > result_ids_list(
			database.GetNumberChunks());
	if (presearch_results_list.size() < queries.GetNumberOfSequences()) {
		presearch_results_list.resize(queries.GetNumberOfSequences());
	}

	size_t max_result = 0;
	for (uint32_t query_id = 0; query_id < queries.GetNumberOfSequences();
			++query_id) {
		results_list[query_id].resize(presearch_results_list[query_id].size());
		max_result = std::max(max_result,
				presearch_results_list[query_id].size());
		for (uint32_t result_id = 0;
				result_id < presearch_results_list[query_id].size();
				++result_id) {
			result_ids_list[presearch_results_list[query_id][result_id].database_chunk_id].push_back(
					make_pair(query_id, result_id));
		}
	}

	boost::thread_group threads;
	boost::mutex next_result_ids_list_i_mutex;
	boost::mutex next_query_id_mutex;
	int number_threads = parameters.number_threads;
	if (database_preload) {
		++number_threads;
	}
	boost::barrier all_barrier(number_threads);
	uint32_t next_result_ids_list_i = 0;
	uint32_t next_query_id = 0;
	for (size_t i = 0; i < parameters.number_threads; ++i) {
		AlignerBuildResultsThread<DatabaseType> t;
		AlignerBuildResultsThread<DatabaseType>::Parameters p;
		p.thread_id = i;
		p.next_result_ids_list_i_mutex = &next_result_ids_list_i_mutex;
		p.next_query_id_mutex = &next_query_id_mutex;
		p.next_result_ids_list_i = &next_result_ids_list_i;
		p.next_query_id = &next_query_id;
		p.queries = &queries;
		p.database = &database;
		p.presearch_results_list = &presearch_results_list[0];
		p.results_list = &results_list[0];
		p.result_ids_list = &result_ids_list[0];
		p.gapped_extension_cutoff = gapped_extension_cutoff;
		p.aligningParameters = &parameters;
		p.all_barrier = &all_barrier;
		threads.create_thread(
				boost::bind(&AlignerBuildResultsThread<DatabaseType>::Run, t,
						p));
	}
	if (database_preload) {
		AlignerBuildResultsThread<DatabaseType> t;
		AlignerBuildResultsThread<DatabaseType>::Parameters p;
		p.thread_id = -1;
		p.next_result_ids_list_i_mutex = &next_result_ids_list_i_mutex;
		p.next_query_id_mutex = &next_query_id_mutex;
		p.next_result_ids_list_i = &next_result_ids_list_i;
		p.next_query_id = &next_query_id;
		p.queries = &queries;
		p.database = &database;
		p.presearch_results_list = &presearch_results_list[0];
		p.results_list = &results_list[0];
		p.result_ids_list = &result_ids_list[0];
		p.gapped_extension_cutoff = gapped_extension_cutoff;
		p.aligningParameters = &parameters;
		p.all_barrier = &all_barrier;
		threads.create_thread(
				boost::bind(&AlignerBuildResultsThread<DatabaseType>::Run, t,
						p));
	}
	threads.join_all();
}

