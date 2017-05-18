/*
 * mpi_common.h
 * Created on 2017/05/16
 *    Auther:goto
 */

#ifndef MPI_COMMON_H_
#define MPI_COMMON_H_

#include"aligner_common.h"
#include<string>

class MPICommon{
 public:
	struct MPIParameter{
		int rank;
		int loadedDatabaseChunk;
		int loadedQueryBlock;
   	
	};
	
	typedef AlignerCommon::AligningCommonParameters AligningParameters;
	typedef AlignerCommon::Result Result;
	void Run(std::string &queries_filename, std::string &database_filename,
			 std::string &output_filename,
			 AligningParameters &parameter,MPIParameter &mpi_parameter);
	void RunGPU(std::string &queries_filename,std::string &database_filename,
				std::string &output_file,
				AligningParameters &parameter,MPIParameter &mpi_parameter);
 private:
	void RunMaster(std::string &queries_filename,std::string &database_filename,
				   std::string &output_filename,
				   AligningParameters &parameter,MPIParameter &mpi_parameter);
	void RunWorker(AligningParameters &parameter,MPIParameter &mpi_parameter);

	void RunWorkerGPU(AligningParameters &parameter,MPIParameter &mpi_parameter);

	//	void RequestTask();
	//void LoadResource();
	
	
	
};



#endif
