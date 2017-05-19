/*
 * query.cpp
 *
 *  Created on: 2009/10/09
 *      Author: shu
 */

#include "fasta_sequence_reader.h"
#include "sequence.h"
#include "alphabet_coder.h"
#include "dna_sequence.h"
#include "dna_type.h"
#include "sequence_no_filter.h"
#include "sequence_seg_filter.h"
#include "protein_type.h"
#include "statistics.h"
#include "translator.h"
#include "translated_dna_query.h"
#include "protein_query.h"
#include "protein_type.h"
#include "queries.h"
#include "logger.h"
#include <cmath>
#include <cstddef>
#include <algorithm>
#include <numeric>
#include <typeinfo>
#include <string>
#include <stdint.h>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <functional>
#include <boost/thread.hpp>
#include <fstream>

using namespace std;
Queries::Queries(const Parameters &parameters) :
		parameters_(parameters), max_query_sequence_length_(0), next_query_(QueryPtr()) {

	
}

Queries::Queries(istream &in, const Parameters &parameters) :
		parameters_(parameters), max_query_sequence_length_(0), reader_(in), next_query_(
				QueryPtr()) {
	Next();
}

void Queries::Next() {
	max_query_sequence_length_ = 0;
	SetQueries();
}

Queries::FilterPtr Queries::GetFilter(bool filter,
		std::tr1::shared_ptr<SequenceType> &aligning_sequence_type_ptr) {
	if (filter) {
		if (typeid(*(aligning_sequence_type_ptr)) == typeid(ProteinType)) {
			return FilterPtr(new SequenceSegFilter());
		} else {
			Logger *logger = Logger::GetInstance();
			logger->ErrorLog("can't use a filter for dna sequences");
			exit(1);
		}
	} else {
		return FilterPtr(new SequenceNoFilter());
	}
	return FilterPtr();
}

uint32_t Queries::GetSequenceLength(
		std::tr1::shared_ptr<SequenceType> &file_sequence_type_ptr,
		std::tr1::shared_ptr<SequenceType> &aligning_sequence_type_ptr,
		uint32_t length) {
	if (typeid(*(file_sequence_type_ptr)) == typeid(DnaType)) {
		if (typeid(*(aligning_sequence_type_ptr)) != typeid(ProteinType)) {
			Logger *logger = Logger::GetInstance();
			logger->ErrorLog("can't use aligner for query dna");
			exit(1);
		} else {
			return TranslatedDnaQuery::GetSequenceLength(length);
		}
	} else {
		return ProteinQuery::GetSequenceLength(length);
	}
}




Queries::QueryPtr Queries::BuildQuery(Sequence &sequence,
		Parameters &parameters, Translator &translator,
		SequenceFilterInterface *sequence_filter) {
	if (typeid(*(parameters.file_sequence_type_ptr)) == typeid(DnaType)) {
		if (typeid(*(parameters.aligning_sequence_type_ptr))
				!= typeid(ProteinType)) {
			Logger *logger = Logger::GetInstance();
			logger->ErrorLog("can't use aligner for query dna");
			exit(1);
		} else {
			return Queries::QueryPtr(
					new TranslatedDnaQuery(sequence,
							parameters.sequence_delimiter, translator,
							sequence_filter, parameters.score_matrix,
							parameters.ungapped_karlin_parameters));
		}
	} else {
		return Queries::QueryPtr(
				new ProteinQuery(sequence, parameters.sequence_delimiter,
						sequence_filter, parameters.score_matrix,
						parameters.ungapped_karlin_parameters));
	}
	return QueryPtr();
}

std::ifstream::pos_type Queries::GetNextChunkPosition(std::istream &is){
	string name;
	string sequence;
	
	bool reader_ret;
	FastaSequenceReader reader(is);
	unsigned int loaded_chunk_size;
	cout<<"max_size:"<<parameters_.chunk_size<<endl;
	cout<<"start:"<<is.tellg()<<endl;
	while(1){
		size_t number_queries = 0;
		size_t number_read_sequences = 0;
		size_t queries_offset = 0;
		if (loaded_chunk_size >= parameters_.chunk_size) {
			break;
		}
		//queries_offset = parameters_.number_queries;
		for (number_read_sequences = 0;
			 number_read_sequences < kNumberReadSequences;
			 ++number_read_sequences) {
			reader_ret = reader.Read(name,sequence);
			if (!reader_ret) {
				break;
			}
			uint32_t new_sequence_length = 
				GetSequenceLength(parameters_.file_sequence_type_ptr,
								  parameters_.aligning_sequence_type_ptr,
								  sequence.size());
			loaded_chunk_size += new_sequence_length;
			if (loaded_chunk_size
				>= parameters_.chunk_size) {
				break;
			}
			
		}
		if (loaded_chunk_size >= parameters_.chunk_size) {
			number_queries = number_read_sequences;
			++number_read_sequences;
		} else {
			number_queries = number_read_sequences;
		}

		if (number_read_sequences != kNumberReadSequences) {
			break;
		}
		

	}
	if(is.tellg()==ios_base::end){
		cout<<"eof"<<endl;
	}
	cout<<"current:"<<is.tellg()<<endl;
	return is.tellg();


}

void Queries::RunThreadForSetingQueries(int thread_id,
		ThreadParameters &parameters) {
	vector<string> names(kNumberReadSequences);
	vector<string> sequences(kNumberReadSequences);
	FilterPtr filter = GetFilter(parameters.parameters->filter,
			parameters.parameters->aligning_sequence_type_ptr);
	bool reader_ret;
	vector<pair<uint32_t, QueryPtr> > temp_queries;
	while (1) {
		size_t number_queries = 0;
		size_t number_read_sequences = 0;
		size_t queries_offset = 0;
		{
			boost::unique_lock<boost::mutex> lock(*(parameters.read_mutex));
			if (*(parameters.chunk_size) >= parameters.parameters->chunk_size) {
				break;
			}
			queries_offset = *parameters.number_queries;
			for (number_read_sequences = 0;
					number_read_sequences < kNumberReadSequences;
					++number_read_sequences) {
				reader_ret = parameters.reader->Read(
						names[number_read_sequences],
						sequences[number_read_sequences]);
				if (!reader_ret) {
					break;
				}
				uint32_t new_sequence_length = GetSequenceLength(
						parameters.parameters->file_sequence_type_ptr,
						parameters.parameters->aligning_sequence_type_ptr,
						sequences[number_read_sequences].size());
				*(parameters.chunk_size) += new_sequence_length;
				if (*(parameters.chunk_size)
						>= parameters.parameters->chunk_size) {
					break;
				}
				*parameters.max_query_sequence_length = std::max(
						new_sequence_length,
						*parameters.max_query_sequence_length);
			}
			if (*(parameters.chunk_size) >= parameters.parameters->chunk_size) {
				number_queries = number_read_sequences;
				++number_read_sequences;
			} else {
				number_queries = number_read_sequences;
			}
			*parameters.number_queries += number_queries;
		}

		for (size_t i = 0; i < number_read_sequences; ++i) {
			Sequence s(names[i], sequences[i]);
			QueryPtr new_query_ptr = BuildQuery(s, *(parameters.parameters),
					*(parameters.translator), filter.get());
			assert(new_query_ptr);
			if (i < number_queries) {
				temp_queries.push_back(
						make_pair(queries_offset + i, new_query_ptr));
			} else {
				(*parameters.next_query) = new_query_ptr;
			}
		}
		if (number_read_sequences != kNumberReadSequences) {
			break;
		}
	}
	parameters.barrier->wait();
	if (thread_id == 0) {
		parameters.queries->resize(*parameters.number_queries);
	}
	parameters.barrier->wait();
	for (size_t i = 0; i < temp_queries.size(); ++i) {
		(*parameters.queries)[temp_queries[i].first] = temp_queries[i].second;
	}
}

bool Queries::SetQueries() {
	queries_.clear();
	max_query_sequence_length_ = 0;
	if (parameters_.number_threads <= 1) {
		unsigned int chunk_size = 0;
		string name;
		string sequence;
		bool reader_ret;
		FilterPtr filter = GetFilter(parameters_.filter,
				parameters_.aligning_sequence_type_ptr);
		if (!next_query_) {
			reader_ret = reader_.Read(name, sequence);
			if (reader_ret) {
				Sequence s(name, sequence);
				next_query_ = BuildQuery(s, parameters_, translator_,
						filter.get());
			}
		}
		while (1) {
			if (!next_query_) {
				break;
			}
			max_query_sequence_length_ = max(max_query_sequence_length_,
					next_query_->GetSequenceLength());
			unsigned int new_chunk_size = chunk_size
					+ next_query_->GetSequenceLength();
			if (new_chunk_size >= parameters_.chunk_size) {
				break;
			}
			queries_.push_back(next_query_);
			chunk_size = new_chunk_size;
			reader_ret = reader_.Read(name, sequence);
			if (reader_ret) {
				Sequence s(name, sequence);
				next_query_ = BuildQuery(s, parameters_, translator_,
						filter.get());
			} else {
				next_query_ = QueryPtr();
				break;
			}
		}
		return chunk_size;
	} else {
		unsigned int chunk_size = 0;
		uint32_t number_queries = 0;
		if (next_query_) {
			uint32_t sequence_length = next_query_->GetSequenceLength();
			chunk_size = sequence_length;
			max_query_sequence_length_ = sequence_length;
			queries_.push_back(next_query_);
			number_queries = 1;
		}
		next_query_ = QueryPtr();
		boost::mutex read_mutex;
		boost::barrier barrier(parameters_.number_threads);
		boost::thread_group threads;
		ThreadParameters thread_parameters;
		thread_parameters.chunk_size = &chunk_size;
		thread_parameters.max_query_sequence_length =
				&max_query_sequence_length_;
		thread_parameters.number_queries = &number_queries;
		thread_parameters.next_query = &next_query_;
		thread_parameters.parameters = &parameters_;
		thread_parameters.queries = &queries_;
		thread_parameters.read_mutex = &read_mutex;
		thread_parameters.barrier = &barrier;
		thread_parameters.reader = &reader_;
		thread_parameters.translator = &translator_;
		for (int i = 0; i < parameters_.number_threads; ++i) {
			threads.create_thread(
					boost::bind(&RunThreadForSetingQueries, i,
							thread_parameters));
		}
		threads.join_all();
		return chunk_size;
	}
}

