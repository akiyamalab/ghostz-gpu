/*
 * database_build_main.cpp
 *
 *  Created on: 2013/10/24
 *      Author: shu
 */

#include "database_build_main.h"
#include <iostream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <unistd.h>
#include <tr1/memory>
#include <sstream>
#include <cstdio>
#include "protein_type.h"
#include "dna_type.h"
#include "logger.h"
#include "score_matrix_reader.h"
#include "reduced_alphabet_file_reader.h"
#include "aligner.h"

using namespace std;

DatabaseBuildMain::DatabaseBuildMain() {

}

DatabaseBuildMain::~DatabaseBuildMain() {

}

int DatabaseBuildMain::Run(int argc, char* argv[]) {
	Logger *logger = Logger::GetInstance();
	logger->Log("building database...");
	Aligner aligner;
	try {
		string input_filename;
		string database_filename;
		Aligner::DatabaseParameters parameters;
		bool ret = BuildParameters(argc, argv, input_filename,
				database_filename, parameters);
		if (ret) {
			aligner.BuildDatabase(input_filename, database_filename,
					parameters);
			logger->Log("finished");
			return 0;
		}
	} catch (exception &e) {
		Logger *logger = Logger::GetInstance();
		logger->ErrorLog(e.what());
	}
	return 1;
}

bool DatabaseBuildMain::BuildParameters(int argc, char* argv[],
		string &input_filename, string &database_filename,
		Aligner::DatabaseParameters &parameters) {
	const string default_hash_alphabet_sets = "A KR EDNQ C G H ILVM FYW P ST";
	const string default_similarity_alphabet_sets =
			"A KR EDNQ C G H ILVM FYW P ST";
//const string default_similarity_alphabet_sets = "LVIMC ASGTP FYW EDNQ KRH";
	const string default_protein_matrix_name = "BLOSUM62";
	const string default_protein_matrix =
			"   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *\nA  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4\n R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4\n N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4\n D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4\n C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4\n Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4\n E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4\n G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4\n H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4\n I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4\n L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4\n K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4\n M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4\n F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4\n P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4\n S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4\n T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4\n W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4\n Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4\n V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4\n B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4\n Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4\n X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4\n * -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1\n";

	int c;
	extern char *optarg;
	extern int optind;
	optind = 1;
	parameters.clustering = true;
	parameters.clustering_subsequence_length = 10;
	parameters.number_threads = 1;
	parameters.chunk_size = 1 << 30;
	parameters.seed_threshold = 39;
	parameters.sequence_type_ptr = std::tr1::shared_ptr<SequenceType>(
			new ProteinType());
	parameters.hash_alphabet_sets.resize(0);
	Logger *logger = Logger::GetInstance();
	parameters.chunk_build_option = DatabaseType::ChunkBuildOptionAllBuild;
	istringstream default_protein_score_matrix_is(default_protein_matrix);
	ScoreMatrixReader score_matrix_reader;
	vector<int> matrix;
	unsigned int number_letters;
	score_matrix_reader.Read(default_protein_score_matrix_is,
			*(parameters.sequence_type_ptr), matrix, number_letters);
	parameters.score_matrix = ScoreMatrix(default_protein_matrix_name,
			&matrix[0], number_letters);
	while ((c = getopt(argc, argv, "a:c:C:i:l:L:o:s:t:")) != -1) {
		switch (c) {
		case 'a':
			parameters.number_threads = atoi(optarg);
			break;
		case 'c':
			if (atoi(optarg) == -1) {
				parameters.chunk_build_option =
						DatabaseType::ChunkBuildOptionNotBuild;
			} else if (atoi(optarg) >= 0) {
				parameters.chunk_build_option = atoi(optarg);
			} else {
				logger->ErrorLog("invalid option, -c is >= -1.");
				return false;
			}
			break;
		case 'C':
			if (optarg[0] == 'T') {
				parameters.clustering = true;
			} else if (optarg[0] == 'F') {
				parameters.clustering = false;
			} else {
				logger->ErrorLog("invalid option, -C is T or F.");
				return false;
			}
			break;

		case 'i':
			input_filename = optarg;
			break;

		case 'l':
			parameters.chunk_size = atoi(optarg);
			break;

		case 'L':
			parameters.clustering_subsequence_length = atoi(optarg);
			break;

		case 'o':
			database_filename = optarg;
			break;

		case 's':
			parameters.seed_threshold = atoi(optarg);
			break;

		case 't':
			if (optarg[0] == 'p') {
				parameters.sequence_type_ptr =
						std::tr1::shared_ptr<SequenceType>(new ProteinType());
			} else {
				logger->ErrorLog("invalid option, -q is p.");
				return false;
			}
			break;
		default:
			logger->ErrorLog("\nTry `ghostz --help' for more information.");
			return false;
		}
	}

	istringstream hash_in(default_hash_alphabet_sets);
	ReducedAlphabetFileReader reduced_alphabet_reader;
	reduced_alphabet_reader.Read(hash_in, parameters.hash_alphabet_sets);
	istringstream similarity_in(default_similarity_alphabet_sets);
	reduced_alphabet_reader.Read(similarity_in,
			parameters.similarity_alphabet_sets);
	return true;
}

