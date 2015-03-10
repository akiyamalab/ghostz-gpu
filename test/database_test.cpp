/*
 * database_test.cpp
 *
 *  Created on: 2012/11/08
 *      Author: shu
 */

#include <gtest/gtest.h>
#include <string>
#include <stdint.h>
#include <fstream>
#include <tr1/memory>
#include "../src/alphabet_coder.h"
#include "../src/score_matrix.h"
#include "../src/score_matrix_reader.h"
#include "../src/sequence_type.h"
#include "../src/protein_type.h"
#include "../src/database.h"
#include "../src/reduced_alphabet_file_reader.h"
#include "../src/reduced_alphabet_coder.h"

using namespace std;

class DatabaseTest: public ::testing::Test {
protected:
	typedef Database<SeedSearcher> DatabaseType;
	virtual void SetUp() {
		in_.open("./test/test_protein.fa");
		DatabaseType::Parameters parameters;
		ProteinType protein_type;
		AlphabetCoder coder(protein_type);
		parameters.chunk_size = 60;
		parameters.sequence_type_ptr = std::tr1::shared_ptr<SequenceType>(
				new ProteinType());
		parameters.chunk_build_option = DatabaseType::ChunkBuildOptionAllBuild;
		const string default_protein_matrix_name = "BLOSUM62";
		const string default_protein_matrix =
				"   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *\nA  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4\n R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4\n N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4\n D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4\n C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4\n Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4\n E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4\n G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4\n H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4\n I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4\n L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4\n K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4\n M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4\n F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4\n P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4\n S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4\n T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4\n W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4\n Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4\n V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4\n B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4\n Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4\n X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4\n * -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1\n";

		istringstream default_protein_score_matrix_is(default_protein_matrix);
		ScoreMatrixReader score_matrix_reader;
		vector<int> matrix;
		unsigned int number_letters;
		score_matrix_reader.Read(default_protein_score_matrix_is, protein_type,
				matrix, number_letters);
		parameters.seed_search_parameters_build_parameters.score_matrix =
				ScoreMatrix(default_protein_matrix_name, &matrix[0],
						number_letters);
		parameters.sequence_delimiter = coder.GetMaxCode() + 1;
		parameters.chunk_build_option = DatabaseType::ChunkBuildOptionAllBuild;
		const std::string reduced_alphabet = "A KR EDNQ C G H ILVM FYW P ST";
		std::istringstream in(reduced_alphabet);
		std::vector<std::string> alphabet_sets;
		ReducedAlphabetFileReader reduced_alphabet_reader;
		reduced_alphabet_reader.Read(in, alphabet_sets);
		ReducedAlphabetCoder reduced_alphabet_coder(protein_type,
				alphabet_sets);
		AlphabetCoder::Code max_code = coder.GetMaxRegularLetterCode();
		parameters.seed_search_parameters_build_parameters.clustering = false;
		parameters.seed_search_parameters_build_parameters.subsequence_length =
				0;
		parameters.seed_search_parameters_build_parameters.max_indexing_code =
				max_code;
		parameters.seed_search_parameters_build_parameters.sequence_delimiter =
				parameters.sequence_delimiter;
		parameters.seed_search_parameters_build_parameters.seed_threshold = 2;
		parameters.seed_search_parameters_build_parameters.threshold = 0.0;
		parameters.seed_search_parameters_build_parameters.number_threads = 1;
		parameters.seed_search_parameters_build_parameters.hash_code_map.resize(
				coder.GetMaxCode() + 1,
				reduced_alphabet_coder.GetUnknownCode());
		parameters.seed_search_parameters_build_parameters.similairty_code_map.resize(
				coder.GetMaxCode() + 1,
				reduced_alphabet_coder.GetUnknownCode());
		for (AlphabetCoder::Code code = coder.GetMinCode(); code <= max_code;
				++code) {
			char c = coder.Decode(code);
			AlphabetCoder::Code reduced_code = reduced_alphabet_coder.Encode(c);
			parameters.seed_search_parameters_build_parameters.hash_code_map[code] =
					reduced_code;
			parameters.seed_search_parameters_build_parameters.similairty_code_map[code] =
					reduced_code;
		}
		parameters.seed_search_parameters_build_parameters.max_indexing_code =
				reduced_alphabet_coder.GetMaxRegularLetterCode();
		database_ = std::tr1::shared_ptr<DatabaseType>(
				new DatabaseType(in_, parameters));
	}

	virtual void TearDown() {
		in_.close();
	}
	ifstream in_;
	std::tr1::shared_ptr<DatabaseType> database_;
};

TEST_F(DatabaseTest, GetNumberChunks) {
	EXPECT_EQ(4, database_->GetNumberChunks());
}

TEST_F(DatabaseTest, GetNumberSequencesInChunk) {
	EXPECT_EQ(2, database_->GetNumberSequencesInChunk());
}

TEST_F(DatabaseTest, GetConcatenatedSequenceLength) {
	EXPECT_EQ(53, database_->GetConcatenatedSequenceLength());
}

TEST_F(DatabaseTest, GetDatabaseTotalLenght) {
	EXPECT_EQ(198, database_->GetDatabaseTotalLenght());
}

TEST_F(DatabaseTest, GetNumberTotalSequences) {
	EXPECT_EQ(8, database_->GetNumberTotalSequences());
}

TEST_F(DatabaseTest, GetMaxSequenceLength) {
	EXPECT_EQ(25, database_->GetMaxSequenceLength());
}

TEST_F(DatabaseTest, SaveAndLoad) {
	string database_filename = "database_test";
	database_->Save(database_filename);
	DatabaseType loaded_database(database_filename);
	EXPECT_EQ(database_->GetNumberChunks(), loaded_database.GetNumberChunks());
	EXPECT_EQ(database_->GetNumberSequencesInChunk(),
			loaded_database.GetNumberSequencesInChunk());
	EXPECT_EQ(database_->GetConcatenatedSequenceLength(),
			loaded_database.GetConcatenatedSequenceLength());
	EXPECT_EQ(database_->GetDatabaseTotalLenght(),
			loaded_database.GetDatabaseTotalLenght());
	EXPECT_EQ(database_->GetNumberTotalSequences(),
			loaded_database.GetNumberTotalSequences());
	EXPECT_EQ(database_->GetMaxSequenceLength(),
			loaded_database.GetMaxSequenceLength());
}

