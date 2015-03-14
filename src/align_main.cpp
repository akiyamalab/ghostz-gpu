/*
 * align_main.cpp
 *
 *  Created on: 2013/10/24
 *      Author: shu
 */

#include "align_main.h"
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
#include "aligner_gpu.h"

using namespace std;

AlignMain::AlignMain() {

}

AlignMain::~AlignMain() {

}

int AlignMain::Run(int argc, char* argv[]) {
	Logger *logger = Logger::GetInstance();
	logger->Log("searching...");

	AlignerCommon::AligningCommonParameters parameters;
	string queries_filename;
	string database_filename;
	string output_filename;

	//try {
	bool ret = BuildParameters(argc, argv, queries_filename, database_filename,
			output_filename, parameters);
	if (ret) {
#ifndef GPU
		Aligner aligner;
		aligner.Align(queries_filename, database_filename, output_filename,
				parameters);
#else
		if (parameters.number_gpus == 0) {
			Aligner aligner;
			aligner.Align(queries_filename, database_filename, output_filename,
					parameters);
		} else {
			AlignerGpu aligner;
			aligner.Align(queries_filename, database_filename, output_filename,
					parameters);
		}
#endif
		logger->Log("finished");
		return 0;
	}

	/*
	 } catch (exception &e) {
	 logger->ErrorLog(e.what());
	 }*/

	return 1;
}

bool AlignMain::BuildParameters(int argc, char* argv[], string &input_filename,
		string &database_filename, string &output_filename,
		Aligner::AligningParameters &parameters) {
	int c;
	extern char *optarg;
	extern int optind;
	optind = 1;
	string score_matrix_filename;
	Logger *logger = Logger::GetInstance();
	const string default_protein_matrix_name = "BLOSUM62";
	const string default_protein_matrix =
			"   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *\nA  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4\n R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4\n N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4\n D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4\n C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4\n Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4\n E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4\n G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4\n H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4\n I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4\n L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4\n K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4\n M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4\n F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4\n P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4\n S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4\n T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4\n W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4\n Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4\n V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4\n B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4\n Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4\n X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4\n * -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1\n";

	parameters.number_threads = 1;
	parameters.filter = true;
	parameters.queries_file_sequence_type_ptr = std::tr1::shared_ptr<
			SequenceType>(new ProteinType());
	parameters.aligning_sequence_type_ptr = std::tr1::shared_ptr<SequenceType>(
			new ProteinType());
	parameters.queries_chunk_size = 1 << 27;
	ScoreMatrixReader score_matrix_reader;
	vector<int> matrix;
	unsigned int number_letters;
	istringstream default_protein_score_matrix_is(default_protein_matrix);
	score_matrix_reader.Read(default_protein_score_matrix_is,
			*(parameters.aligning_sequence_type_ptr), matrix, number_letters);
	parameters.score_matrix = ScoreMatrix(default_protein_matrix_name,
			&matrix[0], number_letters);
	parameters.gap_open = 11;
	parameters.gap_extension = 1;
	parameters.number_gpus = -1;
	parameters.normalized_presearched_ungapped_extension_cutoff = 7.0;
	parameters.normalized_presearched_gapped_extension_trigger = 22.0;
	parameters.normalized_presearched_gapped_extension_cutoff = 15.0;
	parameters.normalized_result_gapped_extension_cutoff = 25.0;
	parameters.max_number_results = 10;
	parameters.max_number_one_subject_results = 1;
	while ((c = getopt(argc, argv, "a:b:d:F:g:h:l:i:o:q:t:v:")) != -1) {
		switch (c) {
		case 'a':
			parameters.number_threads = atoi(optarg);
			break;
		case 'b':
			parameters.max_number_results = atoi(optarg);
			break;
		case 'd':
			database_filename = optarg;
			break;
		case 'F':
			if (optarg[0] == 'T') {
				parameters.filter = true;
			} else if (optarg[0] == 'F') {
				parameters.filter = false;
			} else {
				logger->ErrorLog("invalid option, -F is T or F.");
				return false;
			}
			break;
		case 'g':
			parameters.number_gpus = atoi(optarg);
			break;
		case 'l':
			parameters.queries_chunk_size = atoi(optarg);
			break;
		case 'i':
			input_filename = optarg;
			break;
		case 'o':
			output_filename = optarg;
			break;
		case 'q':
			if (optarg[0] == 'p') {
				parameters.queries_file_sequence_type_ptr =
						std::tr1::shared_ptr<SequenceType>(new ProteinType());
			} else if (optarg[0] == 'd') {
				parameters.queries_file_sequence_type_ptr =
						std::tr1::shared_ptr<SequenceType>(new DnaType());
			} else {
				logger->ErrorLog("invalid option, -q is p or d.");
				return false;
			}
			break;
		case 't':
			if (optarg[0] == 'p') {
				parameters.aligning_sequence_type_ptr = std::tr1::shared_ptr<
						SequenceType>(new ProteinType());
			} else if (optarg[0] == 'd') {
				parameters.aligning_sequence_type_ptr = std::tr1::shared_ptr<
						SequenceType>(new DnaType());
			} else {
				logger->ErrorLog("invalid option, -q is p or d.");
				return false;
			}
			break;
		case 'v':
			parameters.max_number_one_subject_results = atoi(optarg);
			break;
		default:
			logger->ErrorLog("\nTry `ghostz --help' for more information.");
			return false;
		}
	}
	return true;
}

