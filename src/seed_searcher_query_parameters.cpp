/*
 * seed_searcher_query_parameters.cpp
 *
 *  Created on: 2013/04/23
 *      Author: shu
 */

#include "seed_searcher_common.h"
#include "seed_searcher_query_parameters.h"

using namespace std;

int SeedSearcherQueryParameters::Build(int thread_id, int number_threads,
		boost::mutex *next_query_id_mutex, boost::barrier *barrier,
		Queries &queries, BuildParameters &parameters) {
	if (thread_id == 0) {
		sequence_delimiter_ = parameters.common_parameters->sequence_delimiter;
		number_using_hash_position_data_list_ = 0;
		temp_hash_position_data_bucket_lists_.resize(number_threads);
		number_using_hash_position_data_lists_.resize(number_threads, 0);
		number_using_hash_position_data_.resize(number_threads + 1, 0);
		next_query_id_ = 0;
		hash_position_data_offsets_.clear();
		ungapped_extension_cutoffs_ptr_ =
				parameters.ungapped_extension_cutoffs_ptr;
		gapped_extension_triggers_ptr_ =
				parameters.gapped_extension_triggers_ptr;
		score_matrix_ = parameters.score_matrix;
		query_sequences_.resize(queries.GetNumberOfSequences());
	}

	barrier->wait();

	TempHashPositionDataBucketList &pushed_temp_hash_position_data_bucket_list =
			temp_hash_position_data_bucket_lists_[thread_id];
	pushed_temp_hash_position_data_bucket_list.resize(number_threads);
	uint32_t subsequence_length =
			parameters.common_parameters->subsequence_length;
	SeedSearcherDatabaseParameters::SequencesIndex::HashFunction &hash_function =
			parameters.common_parameters->hash_function;
	std::vector<uint32_t> query_ids;
	size_t number_queries = queries.GetNumberOfSequences();
	while (!PopQueryIds(number_queries, next_query_id_mutex, &query_ids)) {
		while (!query_ids.empty()) {
			uint32_t query_id = query_ids.back();
			query_ids.pop_back();
			Query *query = queries.GetQuery(query_id);
			AlphabetCoder::Code *query_sequence = query->GetSequence();
			query_sequences_[query_id] = query_sequence;
			uint32_t query_sequence_length = query->GetSequenceLength();
			uint32_t subsequence_offset = 0;
			while (subsequence_offset < query_sequence_length) {
				uint32_t subsequence_end = subsequence_offset;
				for (;
						(subsequence_end < query_sequence_length)
								&& (query_sequence[subsequence_end]
										!= sequence_delimiter_);
						++subsequence_end) {
				}
				if ((subsequence_end - subsequence_offset)
						>= subsequence_length) {
					uint32_t similarity_check_start = subsequence_offset
							+ subsequence_length / 2;
					uint32_t similarity_check_end = subsequence_end
							- subsequence_length / 2;
					for (uint32_t query_sequence_i = subsequence_offset;
							query_sequence_i < subsequence_end;
							++query_sequence_i) {
						SeedSearcherDatabaseParameters::SequencesIndex::HashKey hash =
								0;
						uint32_t hashed_sequence_length = 0;
#if 0
						SeedSearcherDatabaseParameters::SequencesIndex::HashKey ans_hash =
						0;
						uint32_t ans_hashed_sequence_length = 0;
						int ret = parameters.sequences_index->GetHashKey(
								query_sequence + query_sequence_i, &ans_hash,
								&ans_hashed_sequence_length);
						assert(
								ret
								== hash_function.CalculateHash(
										query_sequence
										+ query_sequence_i,
										&hash,
										&hashed_sequence_length));
						assert(hash == ans_hash);
						assert(ans_hashed_sequence_length ==
								hashed_sequence_length);
#endif

						if (!hash_function.CalculateHash(
								query_sequence + query_sequence_i, &hash,
								&hashed_sequence_length)) {

							uint32_t seed_position =
									SeedSearcherCommon::CalculateSeedPosition(
											query_sequence_i,
											hashed_sequence_length);

							TempHashPositionData temp_hash_position_data;
							temp_hash_position_data.hash = hash;
							temp_hash_position_data.data.query_id = query_id;
							temp_hash_position_data.data.position =
									seed_position;
							if (similarity_check_start <= seed_position
									&& query_sequence_i
											< similarity_check_end) {
								temp_hash_position_data.data.similarity_check =
										true;
							} else {
								temp_hash_position_data.data.similarity_check =
										false;
							}
							pushed_temp_hash_position_data_bucket_list[hash
									% number_threads].push_back(
									temp_hash_position_data);
						}
					}
				}
				subsequence_offset = subsequence_end + 1;
			}
		}
	}
	barrier->wait();
	typedef std::tr1::unordered_map<
			SeedSearcherDatabaseParameters::SequencesIndex::HashKey, uint32_t> HashMap;
	HashMap hash_map;
	std::vector<HashPositionDataList> temp_hash_position_data_lists;
	uint32_t &temp_number_using_hash_position_data_list =
			number_using_hash_position_data_lists_[thread_id];

	for (int list_i = 0; list_i < number_threads; ++list_i) {
		TempHashPositionDataBucket &bucket =
				temp_hash_position_data_bucket_lists_[list_i][thread_id];
		for (size_t data_i = 0; data_i < bucket.size(); ++data_i) {
			TempHashPositionData &data = bucket[data_i];
			HashMap::iterator find_it = hash_map.find(data.hash);
			uint32_t hash_position_data_list_id = 0;
			if (find_it == hash_map.end()) {
				hash_map.insert(
						std::make_pair(data.hash,
								temp_number_using_hash_position_data_list));
				++temp_number_using_hash_position_data_list;
				HashPositionDataList hash_position_data_list;
				hash_position_data_list.hash = 0;
				if (temp_hash_position_data_lists.size()
						<= temp_number_using_hash_position_data_list) {
					temp_hash_position_data_lists.push_back(
							hash_position_data_list);
				}
				hash_position_data_list_id =
						temp_number_using_hash_position_data_list - 1;
				temp_hash_position_data_lists[hash_position_data_list_id].hash =
						data.hash;
				temp_hash_position_data_lists[hash_position_data_list_id].position_data.clear();
			} else {
				hash_position_data_list_id = find_it->second;
			}
			assert(
					temp_hash_position_data_lists[hash_position_data_list_id].hash
							== data.hash);
			temp_hash_position_data_lists[hash_position_data_list_id].position_data.push_back(
					data.data);
			++number_using_hash_position_data_[thread_id];
		}
	}
	barrier->wait();
	if (thread_id == 0) {
		temp_hash_position_data_bucket_lists_.clear();
		number_using_hash_position_data_list_ = 0;
		uint32_t temp_number_of_lists = 0;
		uint32_t temp_number_of_data = 0;
		uint32_t prefix_sum_of_data = 0;
		for (int i = 0; i < number_threads; ++i) {
			temp_number_of_lists = number_using_hash_position_data_lists_[i];
			number_using_hash_position_data_lists_[i] =
					number_using_hash_position_data_list_;
			number_using_hash_position_data_list_ += temp_number_of_lists;

			temp_number_of_data = number_using_hash_position_data_[i];
			number_using_hash_position_data_[i] = prefix_sum_of_data;
			prefix_sum_of_data += temp_number_of_data;
		}
		seed_query_ids_.resize(prefix_sum_of_data);
		seed_query_positions_.resize(prefix_sum_of_data);
		hash_position_data_lists_.resize(number_using_hash_position_data_list_);
		hash_position_data_offsets_.resize(
				number_using_hash_position_data_list_ + 1);
		hash_position_data_offsets_[0] = 0;
		for (int i = 0; i < number_threads; ++i) {
			hash_position_data_offsets_[number_using_hash_position_data_lists_[i]] =
					number_using_hash_position_data_[i];
		}
		hash_position_data_offsets_[number_using_hash_position_data_list_] =
				prefix_sum_of_data;
	}
	barrier->wait();
	uint32_t list_offset = number_using_hash_position_data_lists_[thread_id];
	for (uint32_t i = 0; i < temp_hash_position_data_lists.size(); ++i) {
		swap(hash_position_data_lists_[list_offset + i],
				temp_hash_position_data_lists[i]);
	}

#if 0
	for (uint32_t i = 0; i < number_using_hash_position_data_list_; ++i) {
		HashPositionDataList &hash_position_data_list =
		hash_position_data_lists_[i];
		for (uint32_t j = 0; j < hash_position_data_list.check_position_data.size();
				++j) {
			Query *query = queries.GetQuery(
					hash_position_data_list.check_position_data[j].query_id);
			AlphabetCoder::Code *query_sequence = query->GetSequence();
			uint32_t query_sequence_length = query->GetSequenceLength();
			assert(
					hash_position_data_list.check_position_data[j].position
					< query_sequence_length);
			SeedSearcherDatabaseParameters::SequencesIndex::HashKey hash = 0;
			if (!parameters.sequences_index->GetHashKey(
							query_sequence
							+ hash_position_data_list.check_position_data[j].position,
							&hash)) {
				assert(hash == hash_position_data_list.hash);
			}
		}
	}
#endif
	return 0;
}

/*
 int SeedSearcherQueryParameters::Build(Queries &queries,
 BuildParameters &parameters) {
 sequence_delimiter_ = parameters.sequence_delimiter;
 number_using_hash_position_data_list_ = 0;
 score_matrix_ = *parameters.score_matrix;
 ungapped_extension_cutoffs_ptr_ = parameters.ungapped_extension_cutoffs;
 gapped_extension_triggers_ptr_ = parameters.gapped_extension_triggers;
 query_sequences_.clear();
 uint32_t subsequence_length = parameters.subsequence_length;
 typedef std::tr1::unordered_map<
 SeedSearcherDatabaseParameters::SequencesIndex::HashKey, uint32_t> HashMap;
 HashMap hash_map;
 for (uint32_t query_id = 0; query_id < queries.GetNumberSequences();
 ++query_id) {
 Query *query = queries.GetQuery(query_id);
 AlphabetCoder::Code *query_sequence = query->GetSequence();
 uint32_t query_sequence_length = query->GetSequenceLength();
 query_sequences_.push_back(query_sequence);
 uint32_t subsequence_offset = 0;
 while (subsequence_offset < query_sequence_length) {
 uint32_t subsequence_end = subsequence_offset;
 for (;
 (subsequence_end < query_sequence_length)
 && (query_sequence[subsequence_end]
 != sequence_delimiter_);
 ++subsequence_end) {
 }
 if ((subsequence_end - subsequence_offset) >= subsequence_length) {
 uint32_t similarity_check_start = subsequence_offset
 + subsequence_length / 2;
 uint32_t similarity_check_end = subsequence_end
 - subsequence_length / 2;
 for (uint32_t query_sequence_i = subsequence_offset;
 query_sequence_i < subsequence_end;
 ++query_sequence_i) {
 SeedSearcherDatabaseParameters::SequencesIndex::HashKey hash =
 0;
 uint32_t hashed_sequence_length = 0;
 if (!parameters.sequences_index->GetHashKey(
 query_sequence + query_sequence_i, &hash,
 &hashed_sequence_length)) {
 uint32_t seed_position =
 SeedSearcherCommon::CalculateSeedPosition(
 query_sequence_i,
 hashed_sequence_length);
 HashMap::iterator find_it = hash_map.find(hash);
 HashPositionData hash_position_data;
 if (similarity_check_start <= seed_position
 && query_sequence_i < similarity_check_end) {
 hash_position_data.similarity_check = true;
 } else {
 hash_position_data.similarity_check = false;
 }
 hash_position_data.query_id = query_id;
 hash_position_data.position = seed_position;
 if (find_it == hash_map.end()) {
 hash_map.insert(
 std::make_pair(hash,
 number_using_hash_position_data_list_));
 ++number_using_hash_position_data_list_;
 HashPositionDataList hash_position_data_list;
 hash_position_data_list.hash = 0;
 if (hash_position_data_lists_.size()
 <= number_using_hash_position_data_list_) {
 hash_position_data_lists_.push_back(
 hash_position_data_list);
 }
 hash_position_data_lists_[number_using_hash_position_data_list_
 - 1].hash = hash;
 hash_position_data_lists_[number_using_hash_position_data_list_
 - 1].position_data.clear();
 hash_position_data_lists_[number_using_hash_position_data_list_
 - 1].position_data.push_back(
 hash_position_data);

 } else {
 hash_position_data_lists_[find_it->second].position_data.push_back(
 hash_position_data);
 }
 }
 }
 }
 subsequence_offset = subsequence_end + 1;
 }
 }
 #if 0
 for (uint32_t i = 0; i < number_using_hash_position_data_list_; ++i) {
 HashPositionDataList &hash_position_data_list =
 hash_position_data_lists_[i];
 for (uint32_t j = 0; j < hash_position_data_list.position_data.size();
 ++j) {
 Query *query = queries.GetQuery(
 hash_position_data_list.position_data[j].query_id);
 AlphabetCoder::Code *query_sequence = query->GetSequence();
 uint32_t query_sequence_length = query->GetSequenceLength();
 assert(
 hash_position_data_list.position_data[j].position
 < query_sequence_length);
 SeedSearcherDatabaseParameters::SequencesIndex::HashKey hash = 0;
 if (!parameters.sequences_index->GetHashKey(
 query_sequence
 + hash_position_data_list.position_data[j].position,
 &hash)) {
 assert(hash == hash_position_data_list.hash);
 }
 }
 }
 #endif

 return 0;
 }
 */
AlphabetCoder::Code SeedSearcherQueryParameters::GetSequenceDelimiter() {
	return sequence_delimiter_;
}
std::vector<int> &SeedSearcherQueryParameters::GetUngappedExtensionCutoffs() {
	return *ungapped_extension_cutoffs_ptr_;
}
std::vector<int> &SeedSearcherQueryParameters::GetGappedExtensionTriggers() {
	return *gapped_extension_triggers_ptr_;
}
ScoreMatrix &SeedSearcherQueryParameters::GetScoreMatrix() {
	return score_matrix_;
}

uint32_t SeedSearcherQueryParameters::GetNumberOfHashPositionDataLists() {
	return number_using_hash_position_data_list_;
}

std::vector<AlphabetCoder::Code *> &SeedSearcherQueryParameters::GetQuerySequences() {
	return query_sequences_;
}

std::vector<SeedSearcherQueryParameters::HashPositionDataList> &SeedSearcherQueryParameters::GetHashPositionDataLists() {
	return hash_position_data_lists_;
}

int SeedSearcherQueryParameters::PopQueryIds(size_t number_queries,
		boost::mutex *next_query_id_mutex, std::vector<uint32_t> *query_ids) {
	uint32_t start_id = 0;
	uint32_t end_id = 0;
	assert(query_ids->empty());
	{
		boost::unique_lock<boost::mutex> lock(*next_query_id_mutex);
		start_id = next_query_id_;
		end_id = std::min(number_queries,
				size_t(start_id + kQueryIdsIncrement));
		next_query_id_ = end_id;
	}
	for (uint32_t i = start_id; i < end_id; ++i) {
		query_ids->push_back(i);
	}
	return query_ids->empty() ? 1 : 0;
}

