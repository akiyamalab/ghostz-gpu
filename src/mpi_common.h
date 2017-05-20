/*
 * mpi_common.h
 * Created on 2017/05/16
 *    Auther:goto
 */

#ifndef MPI_COMMON_H_
#define MPI_COMMON_H_

#include "aligner_common.h"
#include "resource_deliverer.h"
#include "mpi_resource.h"
#include <string>
#include <fstream>

class MPICommon{
 public:
	struct MPIParameter{
		int rank;
		int size;
		int loadedDatabaseChunk;
		int loadedQueryBlock;
		
	};
	/*
	struct QueryResource{
		char *data;
		int ptr;
		int size;
		int chunk_id;
		bool available;
	};
	
	struct MasterResources{
		std::string query_filename;
		std::vector<QueryResource> query_list;
		
	};

	struct WorkerResources{
		std::vector<QueryResource> query_list;
	};
	
	*/
	typedef MPIResource::QueryResource QueryResource;
	typedef MPIResource::MasterResources MasterResources;
	typedef MPIResource::WorkerResources WorkerResources;
	
	typedef AlignerCommon::AligningCommonParameters AligningParameters;
	typedef AlignerCommon::Result Result;
	
	
	void debug(int argc,char *argv[]);

	void Run(std::string &queries_filename, std::string &database_filename,
			 std::string &output_filename,
			 AligningParameters &parameter,MPIParameter &mpi_parameter);
	void RunGPU(std::string &queries_filename,std::string &database_filename,
				std::string &output_file,
				AligningParameters &parameter,MPIParameter &mpi_parameter);
	static void LoadQueryResource(MasterResources &resources, int chunk_id);
	static void UnloadQueryResource(MasterResources &resources,int chunk_id);
private:
	
	
	void RunMaster(std::string &queries_filename,std::string &database_filename,
				   std::string &output_filename,
				   AligningParameters &parameter,MPIParameter &mpi_parameter);
	void RunWorker(AligningParameters &parameter,MPIParameter &mpi_parameter);

	void RunWorkerGPU(AligningParameters &parameter,MPIParameter &mpi_parameter);
	
	void SetupMasterResources(std::string &queries_filename, std::string &database_filename,
							  MasterResources &resources, AligningParameters &parameter,
							  MPIParameter &mpi_parameter);
	
	void SetupWorkerResources(WorkerResources &resources);

	bool BuildParameters(int argc,char *argv[],std::string &input_filename,
						 std::string &database_filename,std::string &output_filename,
						 AligningParameters &par1ameters);
	
	void BuildQueryChunkPointers(std::string &queries_filename,
								 std::vector<int> &chunk_pointer_list,
								 std::vector<int> &chunk_size_list,AligningParameters &parameter);
	//	void RequestTask();
	//void LoadResource();
	
	
	
};



#endif
