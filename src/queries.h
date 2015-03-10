/*
 * queries.h
 *
 *  Created on: 2009/10/09
 *      Author: shu
 */

#ifndef QUERIES_H_
#define QUERIES_H_

#include "sequence.h"
#include "fasta_sequence_reader.h"
#include "sequence_filter_interface.h"
#include "query.h"
#include "statistics.h"
#include "translator.h"
#include <vector>
#include <string>
#include <stdint.h>
#include <tr1/memory>
#include <boost/thread.hpp>

class Queries {
public:
	typedef struct {
		bool filter;
		std::tr1::shared_ptr<SequenceType> file_sequence_type_ptr;
		std::tr1::shared_ptr<SequenceType> aligning_sequence_type_ptr;
		Statistics::KarlinParameters ungapped_karlin_parameters;
		unsigned int chunk_size;
		ScoreMatrix score_matrix;
		AlphabetCoder::Code sequence_delimiter;
		int number_threads;
	} Parameters;

	Queries(std::istream &in, const Parameters &parameters);

	virtual ~Queries() {
	}

	void Next();

	Query *GetQuery(uint32_t id) {
		return queries_[id].get();
	}
	size_t GetNumberOfSequences() {
		return queries_.size();
	}

	uint32_t GetMaxQuerySequenceLength() {
		return max_query_sequence_length_;
	}

private:
	typedef std::tr1::shared_ptr<Query> QueryPtr;
	typedef std::tr1::shared_ptr<SequenceFilterInterface> FilterPtr;
	struct ThreadParameters {
		uint32_t *number_queries;
		unsigned int *chunk_size;
		Parameters *parameters;
		std::vector<QueryPtr> *queries;
		uint32_t *max_query_sequence_length;
		Translator *translator;
		FastaSequenceReader *reader;
		QueryPtr *next_query;
		boost::mutex *read_mutex;
		boost::barrier *barrier;
	};

	static const uint32_t kNumberReadSequences = 1 << 5;
	static FilterPtr GetFilter(bool filter,
			std::tr1::shared_ptr<SequenceType> &aligning_sequence_type_ptr);
	static uint32_t GetSequenceLength(
			std::tr1::shared_ptr<SequenceType> &file_sequence_type_ptr,
			std::tr1::shared_ptr<SequenceType> &aligning_sequence_type_ptr,
			uint32_t length);
	static QueryPtr BuildQuery(Sequence &sequence, Parameters &parameters,
			Translator &translator, SequenceFilterInterface *sequence_filter);
	static void RunThreadForSetingQueries(int thread_id,
			ThreadParameters &parameters);
	bool SetQueries();

	Parameters parameters_;
	std::vector<QueryPtr> queries_;
	uint32_t max_query_sequence_length_;
	Translator translator_;
	FastaSequenceReader reader_;
	QueryPtr next_query_;
};

#endif /* QUERIES_H_ */
