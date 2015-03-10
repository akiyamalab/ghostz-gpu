#include "../src/dna_type.h"
#include "../src/alphabet_coder.h"
#include "../src/score_matrix.h"
#include "../src/edit_blocks.h"
#include "../src/chain_filter.h"
#include "../src/seed_searcher_common.h"
#include <gtest/gtest.h>

using namespace std;

TEST(ChainFilterTest, ChainFilterConnected) {
	DnaType type;
	AlphabetCoder coder(type);
	AlphabetCoder::Code delimiter_code = coder.GetMaxCode() + 1;
	ScoreMatrix score_matrix("test_matrix", 4, 2, -1);

	string seq0("AAAAAAC");
	vector<AlphabetCoder::Code> encoded_seq0(seq0.size() + 2);
	encoded_seq0[0] = delimiter_code;
	coder.Encode(&seq0[0], seq0.size(), &encoded_seq0[1]);
	encoded_seq0[encoded_seq0.size() - 1] = delimiter_code;

	string seq1("AAACAAG");
	vector<AlphabetCoder::Code> encoded_seq1(seq1.size() + 2);
	encoded_seq1[0] = delimiter_code;
	coder.Encode(&seq1[0], seq1.size(), &encoded_seq1[1]);
	encoded_seq1[encoded_seq1.size() - 1] = delimiter_code;

	vector<SeedSearcherCommon::Hit> hits;
	SeedSearcherCommon::Hit h;
	h.query_sequence_position = 1;
	h.database_sequence_position = 1;
	hits.push_back(h);
	h.query_sequence_position = 3;
	h.database_sequence_position = 3;
	hits.push_back(h);

	ChainFilter chain_filter(delimiter_code, score_matrix);
	chain_filter.Filter(&encoded_seq0[0], encoded_seq0.size(), 1,
			&encoded_seq1[0], hits);

	EXPECT_EQ(1, hits.size());
}

TEST(ChainFilterTest, ChainFilterReuse) {
	DnaType type;
	AlphabetCoder coder(type);
	AlphabetCoder::Code delimiter_code = coder.GetMaxCode() + 1;
	ScoreMatrix score_matrix("test_matrix", 4, 2, -1);

	string seq0("AAAAAAC");
	vector<AlphabetCoder::Code> encoded_seq0(seq0.size() + 2);
	encoded_seq0[0] = delimiter_code;
	coder.Encode(&seq0[0], seq0.size(), &encoded_seq0[1]);
	encoded_seq0[encoded_seq0.size() - 1] = delimiter_code;

	string seq1("AACAAAG");
	vector<AlphabetCoder::Code> encoded_seq1(seq1.size() + 2);
	encoded_seq1[0] = delimiter_code;
	coder.Encode(&seq1[0], seq1.size(), &encoded_seq1[1]);
	encoded_seq1[encoded_seq1.size() - 1] = delimiter_code;

	vector<SeedSearcherCommon::Hit> hits;
	SeedSearcherCommon::Hit h;
	h.query_sequence_position = 1;
	h.database_sequence_position = 1;
	hits.push_back(h);
	h.query_sequence_position = 2;
	h.database_sequence_position = 2;
	hits.push_back(h);

	ChainFilter chain_filter(delimiter_code, score_matrix);
	chain_filter.Filter(&encoded_seq0[0], encoded_seq0.size(), 1,
			&encoded_seq1[0], hits);

	hits.clear();
	h.query_sequence_position = 1;
	h.database_sequence_position = 1;
	hits.push_back(h);
	h.query_sequence_position = 2;
	h.database_sequence_position = 2;
	hits.push_back(h);
	chain_filter.Filter(&encoded_seq0[0], encoded_seq0.size(), 1,
			&encoded_seq1[0], hits);
	EXPECT_EQ(1, hits.size());
}


TEST(ChainFilterTest, ChainFilterUnconnected) {
	DnaType type;
	AlphabetCoder coder(type);
	AlphabetCoder::Code delimiter_code = coder.GetMaxCode() + 1;
	ScoreMatrix score_matrix("test_matrix", 4, 2, -1);

	string seq0("AAAAAAC");
	vector<AlphabetCoder::Code> encoded_seq0(seq0.size() + 2);
	encoded_seq0[0] = delimiter_code;
	coder.Encode(&seq0[0], seq0.size(), &encoded_seq0[1]);
	encoded_seq0[encoded_seq0.size() - 1] = delimiter_code;

	string seq1("AACAAAG");
	vector<AlphabetCoder::Code> encoded_seq1(seq1.size() + 2);
	encoded_seq1[0] = delimiter_code;
	coder.Encode(&seq1[0], seq1.size(), &encoded_seq1[1]);
	encoded_seq1[encoded_seq1.size() - 1] = delimiter_code;

	vector<SeedSearcherCommon::Hit> hits;
	SeedSearcherCommon::Hit h;
	h.query_sequence_position = 1;
	h.database_sequence_position = 1;
	hits.push_back(h);
	h.query_sequence_position = 4;
	h.database_sequence_position = 4;
	hits.push_back(h);

	ChainFilter chain_filter(delimiter_code, score_matrix);
	chain_filter.Filter(&encoded_seq0[0], encoded_seq0.size(), 1,
			&encoded_seq1[0], hits);

	EXPECT_EQ(2, hits.size());
}
