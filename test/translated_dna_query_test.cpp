#include <gtest/gtest.h>
#include <stdint.h>
#include <iostream>
#include "../src/fasta_sequence_reader.h"
#include "../src/sequence.h"
#include "../src/translator.h"
#include "../src/sequence_no_filter.h"
#include "../src/score_matrix_reader.h"
#include "../src/statistics.h"
#include "../src/protein_type.h"
#include "../src/translated_dna_query.h"

using namespace std;

class TranslatedDnaQueryTest: public ::testing::Test {
protected:
	typedef std::tr1::shared_ptr<Query> QueryPtr;

	virtual void SetUp() {
		Sequence sequence("test0",
				"AGCGAGAGCGAGTGGTGGGCAAAAAAACCTATACTGCAAAATTTTATGAAAGGTGCTTATTGTCCTCTGAATGAT");
		ProteinType protein_type;
		Translator translator;
		SequenceNoFilter no_filter;
		const string default_protein_matrix_name = "BLOSUM62";
		const string default_protein_matrix =
				"   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *\nA  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4\n R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4\n N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4\n D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4\n C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4\n Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4\n E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4\n G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4\n H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4\n I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4\n L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4\n K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4\n M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4\n F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4\n P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4\n S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4\n T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4\n W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4\n Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4\n V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4\n B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4\n Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4\n X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4\n * -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1\n";
		ScoreMatrixReader score_matrix_reader;
		vector<int> matrix;
		unsigned int number_letters;
		istringstream default_protein_score_matrix_is(default_protein_matrix);
		score_matrix_reader.Read(default_protein_score_matrix_is, protein_type,
				matrix, number_letters);
		ScoreMatrix score_matrix(default_protein_matrix_name, &matrix[0],
				number_letters);
		Statistics statistics(protein_type);
		Statistics::KarlinParameters ungapped_karlin_parameters;
		statistics.CalculateUngappedIdealKarlinParameters(score_matrix,
				&ungapped_karlin_parameters);

		AlphabetCoder coder(protein_type);
		query_ = QueryPtr(
				new TranslatedDnaQuery(sequence, coder.GetMaxCode() + 1,
						translator, &no_filter, score_matrix,
						ungapped_karlin_parameters));

	}

	virtual void TearDown() {
	}
	QueryPtr query_;
};

TEST_F(TranslatedDnaQueryTest, GetDistanceFromStartDelimiterFrame0) {
	EXPECT_EQ(10, query_->GetDistanceFromStartDelimiter(10));
}

TEST_F(TranslatedDnaQueryTest, GetDistanceFromStartDelimiterFrame1) {
	EXPECT_EQ(1, query_->GetDistanceFromStartDelimiter(27));
}

TEST_F(TranslatedDnaQueryTest, GetDistanceFromStartDelimiterFrame3) {
	EXPECT_EQ(1, query_->GetDistanceFromStartDelimiter(79));
}

TEST_F(TranslatedDnaQueryTest, GetDistanceFromEndDelimiterFrame0) {
	EXPECT_EQ(16, query_->GetDistanceFromEndDelimiter(10));
}

TEST_F(TranslatedDnaQueryTest, GetDistanceFromEndDelimiterFrame1) {
	EXPECT_EQ(25, query_->GetDistanceFromEndDelimiter(27));
	cout << (int) query_->GetSequence()[26] << endl;
}

TEST_F(TranslatedDnaQueryTest, GetDistanceFromEndDelimiterFrame3) {
	EXPECT_EQ(25, query_->GetDistanceFromEndDelimiter(79));
}

