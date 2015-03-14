/*
 * aligner.cpp
 *
 *  Created on: 2012/10/05
 *      Author: shu
 */

#include "aligner.h"
#include <string>
#include <vector>
#include <queue>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <assert.h>
#include <limits.h>

#include <boost/thread.hpp>

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
#include "gapped_extender.h"
#include "seed_searcher.h"
#include "seed_searcher_query_parameters.h"
#include "reduced_alphabet_coder.h"
#include "aligner_build_results_thread.h"

using namespace std;

void Aligner::BuildDatabase(string &input_filename, string &database_filename,
		DatabaseParameters paramters) {
	Logger *logger = Logger::GetInstance();
	DatabaseType::Parameters database_parameters;
	BuildDatabaseParameters(paramters, database_parameters);
	ifstream is(input_filename.c_str());
	if (!is) {
		logger->ErrorLog("invalid input filename");
	}
	DatabaseType database(is, database_parameters);
	if (database.GetChunkId() < 0) {
		logger->ErrorLog("invalid input file or parameters");
	} else {
		stringstream ss;
		database.Save(database_filename);
		ss.str("");
		ss << "number chunks is " << database.GetNumberChunks() << endl;
		ss << "number of sequences is " << database.GetNumberTotalSequences()
				<< endl;
		ss << "total size is " << database.GetDatabaseTotalLenght() << endl;
		logger->Log(ss.str());
	}
	is.close();
}

void Aligner::Align(string &queries_filename, string &database_filename,
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
		ss << "write results" << endl;
		logger->Log(ss.str());
		WriteOutput(os, queries, database, parameters, results_list);
	}
	queries_is.close();
	os.close();
}

void Aligner::BuildDatabaseParameters(DatabaseParameters &parameters,
		DatabaseType::Parameters &database_parameters) {
	database_parameters.chunk_size = parameters.chunk_size;
	database_parameters.chunk_build_option = parameters.chunk_build_option;
	database_parameters.sequence_type_ptr = parameters.sequence_type_ptr;
	database_parameters.seed_search_parameters_build_parameters.number_threads =
			parameters.number_threads;
	database_parameters.seed_search_parameters_build_parameters.seed_threshold =
			parameters.seed_threshold;
	database_parameters.seed_search_parameters_build_parameters.score_matrix =
			parameters.score_matrix;
	database_parameters.sequence_delimiter =
			AlignerCommon::GetSequenceDelimiter(*parameters.sequence_type_ptr);
	database_parameters.seed_search_parameters_build_parameters.sequence_delimiter =
			database_parameters.sequence_delimiter;
	database_parameters.seed_search_parameters_build_parameters.clustering =
			parameters.clustering;
	database_parameters.seed_search_parameters_build_parameters.subsequence_length =
			parameters.clustering_subsequence_length;
	AlphabetCoder coder(*parameters.sequence_type_ptr);
	if (parameters.hash_alphabet_sets.empty()) {
		database_parameters.seed_search_parameters_build_parameters.max_indexing_code =
				coder.GetMaxRegularLetterCode();
	} else {
		ReducedAlphabetCoder hash_reduced_alphabet_coder(
				*(parameters.sequence_type_ptr), parameters.hash_alphabet_sets);
		AlphabetCoder::Code max_code = coder.GetMaxCode();
		vector<AlphabetCoder::Code> &hash_reduced_code_map =
				database_parameters.seed_search_parameters_build_parameters.hash_code_map;
		hash_reduced_code_map.resize(coder.GetMaxCode() + 1);
		for (AlphabetCoder::Code code = coder.GetMinCode(); code <= max_code;
				++code) {
			char c = coder.Decode(code);
			AlphabetCoder::Code reduced_code =
					hash_reduced_alphabet_coder.Encode(c);
			hash_reduced_code_map[code] = reduced_code;
		}
		database_parameters.seed_search_parameters_build_parameters.max_indexing_code =
				hash_reduced_alphabet_coder.GetMaxRegularLetterCode();

		ReducedAlphabetCoder similarity_reduced_alphabet_coder(
				*(parameters.sequence_type_ptr),
				parameters.similarity_alphabet_sets);
		vector<AlphabetCoder::Code> &similairty_reduced_code_map =
				database_parameters.seed_search_parameters_build_parameters.similairty_code_map;
		similairty_reduced_code_map.resize(coder.GetMaxCode() + 1);
		for (AlphabetCoder::Code code = coder.GetMinCode(); code <= max_code;
				++code) {
			char c = coder.Decode(code);
			AlphabetCoder::Code reduced_code =
					similarity_reduced_alphabet_coder.Encode(c);
			similairty_reduced_code_map[code] = reduced_code;
		}
	}
}

void Aligner::Presearch(Queries &queries, DatabaseType &database,
		AligningParameters &parameters,
		std::vector<std::vector<PresearchedResult> > &results_list) {
	bool database_preload = true;
	Statistics statistics(*(parameters.aligning_sequence_type_ptr));
	Statistics::KarlinParameters gapped_karlin_parameters;
	statistics.CalculateGappedKarlinParameters(parameters.score_matrix,
			parameters.gap_open, parameters.gap_extension,
			&gapped_karlin_parameters);
	int gapped_extension_cutoff = Statistics::Normalized2Nominal(
			parameters.normalized_presearched_gapped_extension_cutoff,
			gapped_karlin_parameters);
	SeedSearcher seed_searcher;
	vector<AlignerPresearchThread::TempChainFilterSeedList> temp_chain_filter_seed_lists(
			parameters.number_threads);
	vector<vector<SeedSearcherGpu::Hit> > chain_filter_seed_lists(
			queries.GetNumberOfSequences());
	vector<vector<std::vector<AlignmentPosition> > > alignment_start_positions_in_database_list(
			parameters.number_threads);
	std::vector<int> ungapped_extension_cutoffs(queries.GetNumberOfSequences());
	std::vector<int> gapped_extension_triggers(queries.GetNumberOfSequences());
	SetQueriesData(queries, parameters, ungapped_extension_cutoffs,
			gapped_extension_triggers);
	SeedSearcher::QueryParameters seed_searcher_query_parameter;
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
	AlignerPresearchThread::ThreadSharedParameters shared_parameters;
	shared_parameters.next_hash_position_data_list_id = 0;
	shared_parameters.next_chain_filter_query_id = 0;
	shared_parameters.queries = &queries;
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
	shared_parameters.results_list = &results_list;
	shared_parameters.presearch_barrier = &presearch_barrier;
	shared_parameters.all_barrier = &all_barrier;

	for (size_t i = 0; i < number_presearch_threads; ++i) {
		AlignerPresearchThread t;
		AlignerPresearchThread::ThreadParameters p;
		p.thread_id = i;
		p.gapped_extension_cutoff = gapped_extension_cutoff;
		p.shared_parameters = &shared_parameters;
		threads.create_thread(boost::bind(&AlignerPresearchThread::Run, t, p));
	}
	if (database_preload) {
		AlignerPresearchThread t;
		AlignerPresearchThread::ThreadParameters p;
		p.thread_id = -1;
		p.gpu_data_list_id = 0;
		p.gapped_extension_cutoff = 0;
		p.gpu_seeds_memory_size = 0;
		p.shared_parameters = &shared_parameters;
		threads.create_thread(boost::bind(&AlignerPresearchThread::Run, t, p));
	}
	threads.join_all();
}

void Aligner::SetQueriesData(Queries &queries, AligningParameters &parameters,
		std::vector<int> &ungapped_extension_cutoffs,
		std::vector<int> &gapped_extension_triggers) {
	for (uint32_t query_id = 0; query_id < queries.GetNumberOfSequences();
			++query_id) {
		Query *query = queries.GetQuery(query_id);
		ungapped_extension_cutoffs[query_id] =
				Statistics::NormalizedCutoff2NominalCutoff(
						parameters.normalized_presearched_ungapped_extension_cutoff,
						query->GetUngappedKarlinParameters());
		gapped_extension_triggers[query_id] = Statistics::Normalized2Nominal(
				parameters.normalized_presearched_gapped_extension_trigger,
				query->GetUngappedKarlinParameters());

	}
}

void Aligner::BuildResults(Queries &queries, DatabaseType &database,
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

void Aligner::WriteOutput(std::ostream &os, Queries &queries,
		DatabaseType &database, AligningParameters &parameters,
		std::vector<std::vector<Result> > &results) {
	AlignerCommon::WriteOutput(os, queries, database, parameters, results);
}

void Aligner::AddResults(DatabaseType &database, AligningParameters &parameters,
		std::vector<std::vector<AlignmentPosition> > &alignment_start_positions_in_database,
		std::vector<PresearchedResult> &added_results,
		std::vector<PresearchedResult> &results) {
	AlignerCommon::AddResults(database, parameters,
			alignment_start_positions_in_database, added_results, results);
}

