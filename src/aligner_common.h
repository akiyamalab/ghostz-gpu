/*
 * aligner_common.h
 *
 *  Created on: Aug 8, 2014
 *      Author: shu
 */

#ifndef ALIGNER_COMMON_H_
#define ALIGNER_COMMON_H_

#include <string>
#include <vector>
#include <tr1/memory>
#include "alphabet_coder.h"
#include "edit_blocks.h"
#include "sequence_type.h"
#include "score_matrix.h"
#include "queries.h"
#include "database.h"
#include "cuda_common.h"

class AlignerCommon {
public:
	template<typename TDatabase>
	struct DatabaseCommonParameters {
		bool clustering;
		uint32_t clustering_subsequence_length;
		uint32_t number_threads;
		std::string filename_prefix;
		unsigned int chunk_size;
		typename TDatabase::ChunkBuildOption chunk_build_option;
		std::tr1::shared_ptr<SequenceType> sequence_type_ptr;
		unsigned int seed_threshold;
		ScoreMatrix score_matrix;
		std::vector<std::string> hash_alphabet_sets;
		std::vector<std::string> similarity_alphabet_sets;
	};

	struct AligningCommonParameters {
		bool filter;
		int gap_open;
		int gap_extension;
		int number_gpus;
		float normalized_presearched_ungapped_extension_cutoff;
		float normalized_presearched_gapped_extension_trigger;
		float normalized_presearched_gapped_extension_cutoff;
		float normalized_result_gapped_extension_cutoff;
		uint32_t number_threads;
		unsigned int queries_chunk_size;
		unsigned int max_number_results;
		unsigned int max_number_one_subject_results;
		ScoreMatrix score_matrix;
		std::tr1::shared_ptr<SequenceType> queries_file_sequence_type_ptr;
		std::tr1::shared_ptr<SequenceType> aligning_sequence_type_ptr;
	};

	typedef cuda_common::Coordinate Coordinate;
#if 0
	typedef struct {
		uint32_t query_position;
		uint32_t database_position;
	}Coordinate;
#endif
	typedef struct {
		int database_chunk_id;
		uint32_t subject_id_in_chunk;
		int score;
		Coordinate hit;
		Coordinate start;
		Coordinate end;
		uint32_t hit_count;
	} PresearchedResult;

	typedef struct {
		uint32_t result_id;
		uint32_t position;
	} AlignmentPosition;

	typedef struct {
		uint32_t database_chunk_id;
		uint32_t subject_id_in_chunk;
		std::string subject_name;
		int score;
		uint32_t alignment_length;
		uint32_t mismatches;
		uint32_t gap_openings;
		uint32_t gap_extensions;
		Coordinate hit;
		Coordinate start;
		Coordinate end;
		uint32_t hit_count;
	} Result;

	static int AlignmentComp(const int a_score, const Coordinate &a_start,
			const Coordinate &a_end, const int b_score,
			const Coordinate &b_start, const Coordinate &b_end);

	static bool IsContained(Coordinate &a_start, Coordinate &a_end,
			Coordinate &b_start, Coordinate &b_end);

	struct AlignmentPositionLessPosition {
		bool operator()(const AlignmentPosition &a1,
				const AlignmentPosition &a2) const {
			return a1.position < a2.position;
		}
	};

	struct PresearchedResultGreaterScore {
		bool operator ()(const PresearchedResult &r1,
				const PresearchedResult &r2) {
			return AlignmentComp(r1.score, r1.start, r1.end, r2.score, r2.start,
					r2.end) > 0;
		}
	};

	struct ResultGreaterScore {
		bool operator ()(const Result &r1, const Result &r2) {
			return AlignmentComp(r1.score, r1.start, r1.end, r2.score, r2.start,
					r2.end) > 0;
		}
	};

	static AlphabetCoder::Code GetSequenceDelimiter(SequenceType &type);

	static void BuildQueriesParameters(AligningCommonParameters &parameters,
				Queries::Parameters &queries_parameters);

	template<typename TDatabase>
	static void AddResults(TDatabase &database,
			AligningCommonParameters &parameters,
			std::vector<std::vector<AlignmentPosition> > &alignment_start_positions_in_database,
			std::vector<PresearchedResult> &added_results,
			std::vector<PresearchedResult> &results);

	template<typename TDatabase>
	static void BuildResult(Query &query, TDatabase &database,
			uint32_t database_chunk_id, uint32_t subject_id_in_chunk, int score,
			Coordinate &hit, Coordinate &start, Coordinate &end,
			std::vector<EditBlocks::EditOpType> &edits, Result &result);

	template<typename TDatabase>
	static void WriteOutput(std::ostream &os, Queries &queries, TDatabase &database,
			AligningCommonParameters &parameters,
			std::vector<std::vector<Result> > &results);
};

template<typename TDatabase>
void AlignerCommon::AddResults(TDatabase &database, AligningCommonParameters &parameters,
		std::vector<std::vector<AlignmentPosition> > &alignment_start_positions_in_database,
		std::vector<PresearchedResult> &added_results,
		std::vector<PresearchedResult> &results) {
	std::vector<PresearchedResult> new_results(0);
	if (alignment_start_positions_in_database.size()
			< database.GetNumberSequencesInChunk()) {
		alignment_start_positions_in_database.resize(
				database.GetNumberSequencesInChunk());
	}
	//sort(added_results.begin(), added_results.end(), PresearchedResultGreaterScore());
	stable_sort(added_results.begin(), added_results.end(), PresearchedResultGreaterScore());

#if 0
// debug //////////////////////////////////////////////////
	cout << "sorted hits" << endl;
	for (uint32_t x = 0; x < hits.size(); ++x) {
		cout << hits[x].score << " "<<hits[x].start.query_position << " " << hits[x].start.database_position << " " << hits[x].end.query_position << " " << hits[x].end.database_position<< " " << hits[x].hit.query_position << " " << hits[x].hit.database_position << endl;
	}
	cout << endl;
///////////////////////////////////////////////////////////
#endif

	std::vector<PresearchedResult>::iterator results_iter = results.begin();
	std::vector<PresearchedResult>::iterator results_end = results.end();
	std::vector<PresearchedResult>::iterator added_results_iter =
			added_results.begin();
	std::vector<PresearchedResult>::iterator added_results_end = added_results.end();
	while (new_results.size() < parameters.max_number_results
			&& (added_results_iter != added_results_end
					|| results_iter != results_end)) {
		PresearchedResult added_alignment;
		if (added_results_iter == added_results_end) {
			added_alignment = *results_iter;
			assert(
					added_alignment.database_chunk_id
							< (int)database.GetNumberChunks());
			++results_iter;
		} else if (results_iter == results_end
				|| AlignmentComp(results_iter->score, results_iter->start,
						results_iter->end, added_results_iter->score,
						added_results_iter->start, added_results_iter->end)
						< 0) {
			added_alignment.database_chunk_id = database.GetChunkId();
			added_alignment.subject_id_in_chunk = database.GetId(
					added_results_iter->start.database_position);
			assert(added_alignment.subject_id_in_chunk < database.GetNumberSequencesInChunk());
			added_alignment.score = added_results_iter->score;
			added_alignment.start = added_results_iter->start;
			added_alignment.end = added_results_iter->end;
			added_alignment.hit = added_results_iter->hit;
			added_alignment.hit_count = added_results_iter->hit_count;
			assert(
					added_alignment.database_chunk_id
							< (int)database.GetNumberChunks());
			++added_results_iter;
		} else {
			added_alignment = *results_iter;
			assert(
					added_alignment.database_chunk_id
							< (int)database.GetNumberChunks());
			++results_iter;
		}
		if (added_alignment.database_chunk_id == database.GetChunkId()) {
			if (alignment_start_positions_in_database[added_alignment.subject_id_in_chunk].size()
					< parameters.max_number_one_subject_results) {
				std::vector<AlignmentPosition> *start_positions =
						&alignment_start_positions_in_database[added_alignment.subject_id_in_chunk];
				AlignmentPosition searched_position;
				searched_position.result_id = 0;
				searched_position.position =
						added_alignment.start.database_position;
				std::vector<AlignmentPosition>::iterator start_positions_iter =
						lower_bound(start_positions->begin(),
								start_positions->end(), searched_position,
								AlignmentPositionLessPosition());
				std::vector<AlignmentPosition>::iterator start_positions_insert_iter =
						start_positions_iter;
				searched_position.position =
						added_alignment.end.database_position;
				std::vector<AlignmentPosition>::iterator start_positions_iter_end =
						lower_bound(start_positions_iter,
								start_positions->end(), searched_position,
								AlignmentPositionLessPosition());

				bool added = true;
				for (; start_positions_iter != start_positions_iter_end;
						++start_positions_iter) {
					PresearchedResult *r_ptr =
							&new_results[start_positions_iter->result_id];
					if (IsContained(r_ptr->start, r_ptr->end,
							added_alignment.start, added_alignment.end)) {
						added = false;
						break;
					}
				}

				if (added == true) {
					AlignmentPosition p;
					p.position = added_alignment.start.query_position;
					p.result_id = new_results.size();
					alignment_start_positions_in_database[added_alignment.subject_id_in_chunk].insert(
							start_positions_insert_iter, p);
					new_results.push_back(added_alignment);
				}
			}
		} else {
			new_results.push_back(added_alignment);
		}
	}
	results.clear();
	results.insert(results.begin(), new_results.begin(), new_results.end());
	for (std::vector<PresearchedResult>::iterator iter = new_results.begin();
			iter != new_results.end(); ++iter) {
		//if (iter->database_chunk_id == database.GetChunkId()) {
		alignment_start_positions_in_database[iter->subject_id_in_chunk].clear();
		//}
	}
#if 0
// debug //////////////////////////////////
	cout << "results" << endl;
	for (uint32_t x = 0; x < results.size(); ++x) {
		//cout << results[x].score << " "<<results[x].start.query_position << " "<< results[x].start.database_position << " "<< results[x].end.query_position << " " << results[x].end.database_position<<endl;
	}
//////////////////////////////////////////
#endif
}

template<typename TDatabase>
void AlignerCommon::BuildResult(Query &query, TDatabase &database,
		uint32_t database_chunk_id, uint32_t subject_id_in_chunk, int score,
		Coordinate &hit, Coordinate &start, Coordinate &end,
		std::vector<EditBlocks::EditOpType> &edits, Result &result) {
	AlphabetCoder::Code *query_sequence = query.GetSequence();
	AlphabetCoder::Code *database_sequence = database.GetConcatenatedSequence();
	uint32_t query_position = start.query_position;
	uint32_t database_position = start.database_position;
	uint32_t m = 0;
	uint32_t o = 0;
	uint32_t e = 0;
	EditBlocks::EditOpType prev_op = EditBlocks::kSubstitution;
	for (std::vector<EditBlocks::EditOpType>::iterator it = edits.begin();
			it != edits.end(); ++it) {
		switch (*it) {
		case EditBlocks::kSubstitution:
			if (query_sequence[query_position]
					!= database_sequence[database_position]) {
				++m;
			}
			++query_position;
			++database_position;
			break;
		case EditBlocks::kGapInSeq0:
			if (prev_op != EditBlocks::kGapInSeq0) {
				++o;
			} else {
				++e;
			}
			++database_position;
			break;
		case EditBlocks::kGapInSeq1:
			if (prev_op != EditBlocks::kGapInSeq1) {
				++o;
			} else {
				++e;
			}
			++query_position;
			break;
		default:
			abort();
			break;
		}
		prev_op = *it;
	}
#if 0
/// debug //////////////////////////////////////////////////
	cout << result.score << " "<<result.start.query_position << " "<< result.start.database_position << " "<< result.end.query_position << " " << result.end.database_position<<endl;
/////////////////////////////////////////////////////////////
#endif
	result.database_chunk_id = database_chunk_id;
	result.subject_id_in_chunk = subject_id_in_chunk;
	result.subject_name = database.GetName(subject_id_in_chunk);
	result.score = score;
	result.alignment_length = edits.size();
	result.mismatches = m;
	result.gap_openings = o;
	result.gap_extensions = e;
	result.hit = hit;
	result.start.query_position = query.GetRealStart(start.query_position);
	result.end.query_position = query.GetRealEnd(end.query_position);
	result.start.database_position = start.database_position
			- (database.GetOffset(subject_id_in_chunk) - 1);
	result.end.database_position = end.database_position
			- (database.GetOffset(subject_id_in_chunk) - 1);
}


template<typename TDatabase>
void AlignerCommon::WriteOutput(std::ostream &os, Queries &queries, TDatabase &database,
		AligningCommonParameters &parameters,
		std::vector<std::vector<Result> > &results) {
	Statistics::KarlinParameters gapped_karlin_parameters;
	double alpha;
	double beta;

	Statistics statistics(*(parameters.aligning_sequence_type_ptr));
	statistics.CalculateGappedKarlinParameters(parameters.score_matrix,
			parameters.gap_open, parameters.gap_extension,
			&gapped_karlin_parameters);
	statistics.CalculateAlphaBeta(parameters.score_matrix, parameters.gap_open,
			parameters.gap_extension, &alpha, &beta);
	for (uint32_t i = 0; i < queries.GetNumberOfSequences(); ++i) {
		Query *query = queries.GetQuery(i);
		std::string query_name = query->GetName();
		//cout << i << " : " << query_name << endl;
		uint64_t search_space = statistics.CalculateSearchSpace(
				query->GetRealSequenceLength(),
				database.GetDatabaseTotalLenght(),
				database.GetNumberTotalSequences(), gapped_karlin_parameters,
				alpha, beta);
		for (std::vector<Result>::iterator it = results[i].begin();
				it != results[i].end(); ++it) {
			os << query_name << "\t";
			os << it->subject_name << "\t";
			os
					<< (1.0
							- (static_cast<float>(it->mismatches
									+ it->gap_openings + it->gap_extensions)
									/ static_cast<float>(it->alignment_length)))
							* 100 << "\t";
			os << it->alignment_length << "\t";
			os << it->mismatches << "\t";
			os << it->gap_openings << "\t";
			os << it->start.query_position << "\t";
			os << it->end.query_position << "\t";
			os << it->start.database_position << "\t";
			os << it->end.database_position << "\t";
			os
					<< Statistics::Nominal2EValue(it->score, search_space,
							gapped_karlin_parameters) << "\t";
			os
					<< Statistics::Nominal2Normalized(it->score,
							gapped_karlin_parameters);
			//os << "\t" << it->score << "\t" << it->hit_count << "\t";
			//os << "\t(" << it->hit.query_position << ", "
			//		<< it->hit.database_position << ")";
			os << '\n';
		}
	}
}

#endif /* ALIGNER_COMMON_H_ */
