/*
 * mpi_common.cpp
 *  Created on 2017/05/16
 *  Author:goto
 */

#include"mpi_common.h"
#include"aligner_mpi.h"
#include"fasta_sequence_reader.h"
#include"protein_type.h"
#include"dna_type.h"
#include"sequence_type.h"
#include"queries.h"
#include"logger.h"
#include"score_matrix_reader.h"
#include"reduced_alphabet_file_reader.h"


#include"mpi.h"
#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
 using namespace std;




void MPICommon::debug(int argc,char *argv[]){
	
if(argc>=2){
		//cout<<argv[1]<<endl;
	}
	string input_file;
	string database_file;
	string output_file;
	AligningParameters parameter;
	BuildParameters(argc,argv,input_file,database_file,output_file,
					parameter);
	
	
	
	ifstream isf(input_file.c_str());
	
	ofstream osf("1chunk.fasta");
	Queries::Parameters queries_parameters;
	
	AlignerCommon::BuildQueriesParameters(parameter,queries_parameters);
	Queries queries(queries_parameters);
	int i=0;
	bool eof_flag=false;
	while(!eof_flag){
		ifstream::pos_type begin = isf.tellg();
		ifstream::pos_type isf_ptr=queries.GetNextChunkPosition(isf);
		if(isf_ptr==-1){
			isf.clear();
			isf.seekg(begin);
			isf_ptr=isf.seekg(0,std::ios_base::end).tellg();
			eof_flag=true;
		   
		}
		
		cout<<"begin,isf_ptr:"<<begin<<","<<isf_ptr<<endl;
		stringstream ss;
		ss<<output_file.c_str()<<"_"<<i;
		ofstream ofs(ss.str().c_str(),ios::binary | ios::out);

		char *array = new char[isf_ptr-begin];
		isf.clear();

		isf.seekg(begin);
		isf.read(array,isf_ptr-begin);
		//cout<<array<<endl;
		ofs.write(array,isf_ptr-begin);
		ofs.flush();
		ofs.close();
		isf.seekg(isf_ptr);
		i++;
		delete [] array;
	}
	
	//cout<<name<<"\n"<<sequence<<endl;
	   uint32_t length;
	
	
	std::tr1::shared_ptr<SequenceType> queries_sequence_type;
	std::tr1::shared_ptr<SequenceType> database_sequence_type;
	queries_sequence_type  = std::tr1::shared_ptr<SequenceType>(new DnaType());
	database_sequence_type = std::tr1::shared_ptr<SequenceType>(new ProteinType());

	//length =Queries::GetSequenceLength(queries_sequence_type,database_sequence_type,
	//							   sequence.size());
	//cout<<length<<endl;
	
	
							   
								   

}
void MPICommon::Run(string &queries_filename,string &database_filename,
					string &output_filename,
					AligningParameters &parameter,MPIParameter &mpi_parameter){
	
	if(mpi_parameter.rank==0){
		RunMaster(queries_filename,database_filename,output_filename,
				  parameter,mpi_parameter);
	}else{
		RunWorker(parameter,mpi_parameter);
		}
}

void MPICommon::RunGPU(string &queries_filename,string &database_filename,
					string &output_filename,
					AligningParameters &parameter,MPIParameter &mpi_parameter){
	
	if(mpi_parameter.rank==0){
		RunMaster(queries_filename,database_filename,output_filename,
				  parameter,mpi_parameter);
	}else{
		RunWorkerGPU(parameter,mpi_parameter);
		}
}

void MPICommon::RunMaster(string &queries_filename,string &database_filename,
						  string &output_filename,
						  AligningParameters &parameter,MPIParameter &mpi_parameter){

	BuildQueryChunkPointers(queries_filename,_chunk_pointer_list,_chunk_size_list,parameter);
	
	int i=0;
	for(i=0;i<_chunk_pointer_list.size();i++){
		cout<<i<<" ptr:"<<_chunk_pointer_list[i]<<" size:"<<_chunk_size_list[i]<<endl;
	}
	
}
void MPICommon::RunWorker(AligningParameters &parameter,MPIParameter &mpi_parameter){

}


void MPICommon::RunWorkerGPU(AligningParameters &parameter,MPIParameter &mpi_parameter){

}

					  
void MPICommon::BuildQueryChunkPointers(string &queries_filename,vector<int> &chunk_pointer_list,
										vector<int> &chunk_size_list,AligningParameters &parameter){
	
	ifstream isf(queries_filename.c_str());
	Queries::Parameters queries_parameters;
	AlignerCommon::BuildQueriesParameters(parameter,queries_parameters);
	Queries queries(queries_parameters);

	int i=0;
	bool eof_flag=false;
	while(!eof_flag){
		ifstream::pos_type begin = isf.tellg();
		ifstream::pos_type isf_ptr=queries.GetNextChunkPosition(isf);
		if(isf_ptr==-1){
			isf.clear();
			isf.seekg(begin);
			isf_ptr=isf.seekg(0,std::ios_base::end).tellg();
			eof_flag=true;
		   
		}
		
		int size=isf_ptr-begin;
		char *array = new char[size];
		isf.clear();

		isf.seekg(begin);
		isf.read(array,size);
		isf.seekg(isf_ptr);
		i++;
		chunk_pointer_list.push_back(begin);
		chunk_size_list.push_back(size);
		delete [] array;
	}
}



bool MPICommon::BuildParameters(int argc, char* argv[], string &input_filename,
								 string &database_filename, string &output_filename,
								 AligningParameters &parameters) {
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
