/*
 * aligner_build_results_thread.h
 *
 *  Created on: Aug 11, 2014
 *      Author: shu
 */

#ifndef ALIGNER_BUILD_RESULTS_THREAD_H_
#define ALIGNER_BUILD_RESULTS_THREAD_H_

#include <boost/thread.hpp>
#include <boost/thread/barrier.hpp>

#include "queries.h"
#include "database.h"
#include "aligner_common.h"
#include "gapped_extender.h"

template<typename TDatabase>
class AlignerBuildResultsThread {
public:
	typedef TDatabase Database;
	struct Parameters {
		boost::mutex *next_result_ids_list_i_mutex;
		boost::mutex *next_query_id_mutex;
		int thread_id;
		uint32_t *next_result_ids_list_i;
		uint32_t *next_query_id;
		Queries *queries;
		Database *database;
		std::vector<AlignerCommon::PresearchedResult> *presearch_results_list;
		std::vector<AlignerCommon::Result> *results_list;
		std::vector<std::pair<uint32_t, uint32_t> > *result_ids_list;
		int gapped_extension_cutoff;
		AlignerCommon::AligningCommonParameters *aligningParameters;
		boost::barrier *all_barrier;
	};
	AlignerBuildResultsThread();
	virtual ~AlignerBuildResultsThread();

	void Run(Parameters &parameters);

private:
	static const uint32_t kResultIdsIIncrement = 1 << 5;
	static const uint32_t kQueryIdsIncrement = 1 << 5;
	int PopResultIdsI(size_t result_ids_list_size,
			std::vector<uint32_t> *result_ids_i_stack);
	int PopQueryIds(size_t number_queries, std::vector<uint32_t> *query_ids);

	Parameters *parameters_;
};

template<typename TDatabase>
AlignerBuildResultsThread<TDatabase>::AlignerBuildResultsThread() :
		parameters_(NULL) {

}

template<typename TDatabase>
AlignerBuildResultsThread<TDatabase>::~AlignerBuildResultsThread() {

}

template<typename TDatabase>
void AlignerBuildResultsThread<TDatabase>::Run(Parameters &parameters) {
	parameters_ = &parameters;
	if (parameters_->thread_id == 0) {
		parameters_->database->ResetChunk();
		*(parameters_->next_query_id) = 0;
	}

	std::vector<uint32_t> result_ids_i_stack;
	parameters_->all_barrier->wait();
	while (parameters_->database->GetChunkId()
			< (int) parameters_->database->GetNumberChunks()) {
		if (parameters_->thread_id == 0) {
			std::cout << "trace back against database chunk "
					<< parameters_->database->GetChunkId() << std::endl;
		}
		if (parameters_->thread_id < 0) {
			typename Database::PreloadTarget target = Database::kName
					| Database::kSequence | Database::kOffset;
			std::cout << "preload thread begin to preload chunk "
					<< parameters_->database->GetChunkId() + 1 << std::endl;
			parameters_->database->Preload(
					parameters_->database->GetChunkId() + 1, target);
		} else {
			const uint32_t database_chunk_id =
					parameters_->database->GetChunkId();
			const AlphabetCoder::Code *database_concatenated_sequence =
					parameters_->database->GetConcatenatedSequence();
			std::vector<std::pair<uint32_t, uint32_t> > &result_ids =
					parameters_->result_ids_list[database_chunk_id];
			const size_t result_ids_list_size = result_ids.size();

			EditBlocks edit_blocks;
			EditBlocks tmp_edit_blocks;
			GappedExtender gapped_extender;
			while (!PopResultIdsI(result_ids_list_size, &result_ids_i_stack)) {
				while (!result_ids_i_stack.empty()) {
					uint32_t result_ids_i = result_ids_i_stack.back();
					result_ids_i_stack.pop_back();
					edit_blocks.Clear();
					assert(result_ids_i < result_ids_list_size);
					uint32_t query_id = result_ids[result_ids_i].first;
					uint32_t result_id = result_ids[result_ids_i].second;
					assert(
							query_id
									< parameters_->queries->GetNumberOfSequences());
					assert(result_id < 10);
					assert(
							result_id
									< (parameters_->results_list)[query_id].size());
					Query *query = parameters_->queries->GetQuery(query_id);
					int sum_score = 0;
					AlignerCommon::Coordinate hit;
					AlignerCommon::Coordinate start;
					AlignerCommon::Coordinate end;
					AlignerCommon::PresearchedResult *presearched_result_ptr =
							&parameters_->presearch_results_list[query_id][result_id];
					hit.query_position =
							presearched_result_ptr->hit.query_position;
					hit.database_position =
							presearched_result_ptr->hit.database_position;
					int score;
					int query_position;
					int database_position;
					gapped_extender.ExtendOneSide(
							query->GetSequence() + hit.query_position - 1,
							hit.query_position - 1,
							database_concatenated_sequence
									+ hit.database_position - 1,
							query->GetSequenceDelimiter(), true,
							parameters_->aligningParameters->score_matrix,
							parameters_->aligningParameters->gap_open,
							parameters_->aligningParameters->gap_extension,
							parameters_->gapped_extension_cutoff, &score,
							&query_position, &database_position,
							&tmp_edit_blocks);
					sum_score += score;
					start.query_position = hit.query_position - 1
							+ query_position;
					start.database_position = hit.database_position - 1
							+ database_position;
					tmp_edit_blocks.Reverse();
					edit_blocks.Add(tmp_edit_blocks);
					gapped_extender.ExtendOneSide(
							query->GetSequence() + hit.query_position,
							query->GetSequenceLength() - hit.query_position,
							database_concatenated_sequence
									+ hit.database_position,
							query->GetSequenceDelimiter(), false,
							parameters_->aligningParameters->score_matrix,
							parameters_->aligningParameters->gap_open,
							parameters_->aligningParameters->gap_extension,
							parameters_->gapped_extension_cutoff, &score,
							&query_position, &database_position,
							&tmp_edit_blocks);
					sum_score += score;
					end.query_position = hit.query_position + query_position;
					end.database_position = hit.database_position
							+ database_position;
					edit_blocks.Add(tmp_edit_blocks);
					std::vector<EditBlocks::EditOpType> edits =
							edit_blocks.ToVector();
					AlignerCommon::BuildResult(*query, *(parameters_->database),
							presearched_result_ptr->database_chunk_id,
							presearched_result_ptr->subject_id_in_chunk,
							sum_score, hit, start, end, edits,
							parameters_->results_list[query_id][result_id]);
				}
			}
		}
		parameters_->all_barrier->wait();
		if (parameters_->thread_id == 0) {
			parameters_->database->NextChunk();
			*(parameters_->next_result_ids_list_i) = 0;
		}
		parameters_->all_barrier->wait();
	}

	parameters_->all_barrier->wait();
	std::vector<uint32_t> query_ids;
	size_t number_queries = parameters_->queries->GetNumberOfSequences();
	while (!PopQueryIds(number_queries, &query_ids)) {
		while (!query_ids.empty()) {
			uint32_t query_id = query_ids.back();
			query_ids.pop_back();
			assert(parameters_->results_list[query_id].size() <= 10);
			sort(parameters_->results_list[query_id].begin(),
					parameters_->results_list[query_id].end(),
					AlignerCommon::ResultGreaterScore());
		}
	}
}

template<typename TDatabase>
int AlignerBuildResultsThread<TDatabase>::PopResultIdsI(size_t result_ids_list_size,
		std::vector<uint32_t> *result_ids_i_stack) {
	uint32_t start_id = 0;
	uint32_t end_id = 0;
	assert(result_ids_i_stack->empty());
	{
		boost::unique_lock<boost::mutex> lock(*(parameters_->next_result_ids_list_i_mutex));
		start_id = *(parameters_->next_result_ids_list_i);
		end_id = std::min(result_ids_list_size,
				size_t(start_id + kResultIdsIIncrement));
		*(parameters_->next_result_ids_list_i) = end_id;
	}
	for (uint32_t i = start_id; i < end_id; ++i) {
		result_ids_i_stack->push_back(i);
	}
	return result_ids_i_stack->empty() ? 1 : 0;
}

template<typename TDatabase>
int AlignerBuildResultsThread<TDatabase>::PopQueryIds(size_t number_queries,
		std::vector<uint32_t> *query_ids) {
	uint32_t start_id = 0;
	uint32_t end_id = 0;
	assert(query_ids->empty());
	{
		boost::unique_lock<boost::mutex> lock(*(parameters_->next_query_id_mutex));
		start_id = *(parameters_->next_query_id);
		end_id = std::min(number_queries,
				size_t(start_id + kQueryIdsIncrement));
		*(parameters_->next_query_id) = end_id;
	}
	for (uint32_t i = start_id; i < end_id; ++i) {
		query_ids->push_back(i);
	}
	return query_ids->empty() ? 1 : 0;
}


#endif /* ALIGNER_BUILD_RESULTS_THREAD_H_ */
