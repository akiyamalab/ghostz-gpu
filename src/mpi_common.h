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
#include "database.h"
#include "load_balancer.h"

#include "mpi.h"
#include <string>
#include <fstream>
#include <sys/time.h>

#define F_GUIDED

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
	typedef MPIResource::DatabaseResource DatabaseResource;
	typedef MPIResource::DatabaseInfo DatabaseInfo;
	typedef MPIResource::MasterResources MasterResources;
	typedef MPIResource::WorkerResources WorkerResources;
	typedef MPIResource::AlignmentTask AlignmentTask;
	
	typedef AlignerCommon::AligningCommonParameters AligningParameters;
	typedef AlignerCommon::Result Result;
	typedef Database<SeedSearcher> DatabaseType;
	typedef Database<SeedSearcherGpu> DatabaseTypeGpu;
	

	void Run(std::string &queries_filename, std::string &database_filename,
			 std::string &output_filename, std::string &tmp_dirname,
			 AligningParameters &parameter, MPIParameter &mpi_parameter);
	void RunGPU(std::string &queries_filename,std::string &database_filename,
				std::string &output_file, std::string &tmp_dirname,
				AligningParameters &parameter,MPIParameter &mpi_parameter);
	static void LoadQueryResource(MasterResources &resources, int chunk_id);
	static void UnloadQueryResource(MasterResources &resources,int chunk_id);
private:
	int count_terminate;
	void RunMaster(std::string &queries_filename,std::string &database_filename,
				   std::string &output_filename,std::string &tmp_dirname,
				   AligningParameters &parameter,MPIParameter &mpi_parameter);
	void RunWorker(std::string &queries_filename,std::string &database_filename,
				   std::string &output_filename,std::string &tmp_dirname,
				   AligningParameters &parameter,MPIParameter &mpi_parameter,
				   bool useGPU);
	

	void SetupDatabaseResourcesMaster(std::string &database_filename,
							  MasterResources &resources, AligningParameters &parameter,
							  MPIParameter &mpi_parameter);
	
	void SetupDatabaseResourcesWorker(std::string &database_filename,std::string &tmp_dirname,WorkerResources &resources,
									  AligningParameters &parameter,MPIParameter &mpi_parameter);

	void SetupQueryResourcesMaster(std::string &queries_filename,MasterResources &resources,
								   AligningParameters &parameter,MPIParameter &mpi_parameter);
	
	void SetupQueryResourcesWorker(WorkerResources &resources,
								   AligningParameters &parameter,MPIParameter &mpi_parameter);
	
	void AcceptCommand(MasterResources &resources,LoadBalancer &balancer);
	
	
	void BuildQueryChunkPointers(std::string &queries_filename,
								 std::vector<uint64_t> &chunk_pointer_list,
								 std::vector<uint64_t> &chunk_size_list,AligningParameters &parameter);
	
	void BuildGuidedQueryChunkPointers(std::string &queries_filename,
									   std::vector<uint64_t> &chunk_pointer_list,
									   std::vector<uint64_t> &chunk_size_list,AligningParameters &parameter);
									   
	
	double timevalToMillisec(struct timeval t){
		return (t.tv_sec *1000.0 + t.tv_usec* 0.001);
	}
	
};



#endif
