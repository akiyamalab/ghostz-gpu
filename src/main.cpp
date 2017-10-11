
#include <cstring>
#include "common.h"
#include "database_build_main.h"
#include "align_main.h"
#include "logger.h"

using namespace std;

void PrintUsage() {
	Logger *logger = Logger::GetInstance();
	logger->Log(
			string(
					"GHOSTZ-GPU - homology search tool. version "
							+ common::kVersion)
							);
	logger->Log(
			"\n"
					"Command and Options\n"
					"db: convert a FASTA file to GHOSTX format database files\n"
					"\n"
					"  ghostz_gpu db [-i dbFastaFile] [-o dbName] [-C clustering] [-l chunkSize]\n"
					"            [-L clusteringSubsequenceLength] [-s seedThreshold] [-a numThreads]\n"
					"\n"
					"  Options:\n"
					"  (Required)\n"
					"    -i STR    Protein sequences in FASTA format for a database\n"
					"    -o STR    The name of the database\n"
					"\n"
					"  (Optional)\n"
					"    -C STR    Clustering, T (enable) or F (disable) [T]\n"
					"    -l INT    Chunk size of the database (bytes) [1073741824 (=1GB)]\n"
					"    -L INT    Length of a subsequence for clustering [10]\n"
					"    -s INT    The seed threshold [39]\n"
			                "    -a INT    The number of threads [1]\n"
					"\n"
					"\n"
					"aln:  Search homologues of queries from database\n"
					"\n"
					"  ghostz_gpu aln [-i queries] [-o output] [-d database] [-v maxNumAliSub]\n"
					"             [-b maxNumAliQue] [-h hitsSize] [-l queriesChunkSize] [-q queryType]\n"
					"             [-t databaseType] [-F filter] [-a numThreads] [-g numGPUs]\n"
					"\n"
					"  Options:\n"
					"  (Required)\n"
					"    -i STR    Sequences in FASTA format\n"
					"    -o STR    Output file\n"
					"    -d STR    database name (must be formatted)\n"
					"\n"
					"  (Optional)\n"
					"    -v INT    Maximum number of alignments for each subject [1]\n"
					"    -b INT    Maximum number of the output for a query [10]\n"
					"\n"
					"    -l INT    Chunk size of the queries (bytes) [134217728 (=128MB)]\n"
					"    -q STR    Query sequence type, p (protein) or d (dna) [p]\n"
					"    -t STR    Database sequence type, p (protein) or d (dna) [p]\n"
					"    -F STR    Filter query sequence, T (enable) or F (disable) [T] \n"
					"    -a INT    The number of threads [1]\n"
					"    -g INT    The number of GPUs [the number of available GPUs]");
}

int main(int argc, char* argv[]) {
	if (argc < 2 || strcmp(argv[1], "-h") == 0
			|| strcmp(argv[1], "--help") == 0) {
		PrintUsage();
		exit(1);
	}

	if (strcmp(argv[1], "db") == 0) {
	        DatabaseBuildMain main;
		return main.Run(argc - 1, argv + 1);
	} else if (strcmp(argv[1], "aln") == 0) {
		AlignMain main;
		return main.Run(argc - 1, argv + 1);
	} else {
		cerr << "unrecognized command " << argv[1] << endl;
		return 1;
	}
	return 0;
}
