/* 
 * mpi_resource.h
 *    Created 2017/05/20
 *     Author:goto
 *
 */
#ifndef _MPI_RESOURCE_H_
#define _MPI_RESOURCE_H_

#include <queue>
#include <stdint.h>
#include "mpi.h"

#include "alphabet_coder.h"
#include "unistd.h"

class MPIResource{
 public:
    struct AlignmentTask{
        int query_chunk;
        int database_chunk;
    } ;
	
	struct QueryResource{
        char *data;
        int ptr;
        uint64_t size;
        int chunk_id;
        bool available;
    };
	
	struct DatabaseInfo{
		int number_chunks;
		uint32_t max_sequence_length;
		uint64_t database_length;
		uint64_t number_sequences;
		AlphabetCoder::Code sequence_delimiter;
	};

	struct DatabaseResource{
		char *inf; //informationFile
		uint64_t inf_size;
		char *nam; //Namesfilename
		uint64_t nam_size;
		char *off; //OffsetFile
		uint64_t off_size;
		char *seq; //sequenceFile
		uint64_t seq_size;
		char *scp;//SeedSearcharCommonParemeter
		uint64_t scp_size;
		char *sdp;//SeedSercharParameter
		uint64_t sdp_size;
		
		
		int chunk_id;
		bool available;
		
	};
	

    struct MasterResources{ 
		std::string query_filename;
		std::vector<QueryResource> query_list;             //size= number of query chunk
		
		std::string database_filename;
		std::vector<DatabaseResource> database_list;       //size= number pf database chunk
		DatabaseInfo database_info;
		std::string output_filename;
		
		std::vector<std::queue<AlignmentTask> > task_pool; //size= number of database chunk 
		std::vector<int> node_loading_database;            //size= number of (master+)worker node
		std::vector<float> task_ratio;                      //size= number of database chunk
		
		std::vector<MPI::Intercomm> comm_list;
		
    };
    struct WorkerResources{
		std::vector<QueryResource> query_list;
		std::string database_filename;
		DatabaseInfo database_info;
		std::vector<DatabaseResource> database_list;
		MPI::Intercomm subgroup_comm;
		
		
    };
	

	enum Command{
        CMD_RequestQuery,
        CMD_RequestDatabase,
        CMD_RequestTask,
        ACK,
        NACK
    };
	
	static int BcastDatabase(DatabaseResource &database, MPI::Intercomm comm, int root);
	static int BcastDatabaseInfo(DatabaseInfo &info,MPI::Intercomm comm ,int root);
	static void BcastLargeData(char **ptr,uint64_t size,MPI::Intercomm comm,int root);
	
	static int AcceptCommand(MasterResources &resources,int *cmd,int *target);
	
	static int RequestQuery(WorkerResources &resources,int target_chunk);
	static int AcceptRequestQuery(MasterResources &resources,int target_chunk,int dst_rank);
	static void RequestTask(int target_chunk,MPIResource::AlignmentTask &task);
	static int AcceptRequestTask(AlignmentTask task,int dst_rank);

	
	static int RecvQuery(QueryResource &query_resource, MPI::Intercomm comm,int src_rank);
	static int SendQuery(QueryResource &query_resource, MPI::Intercomm comm,int dst_rank);
	static void LoadQueryResource(MasterResources &resources, int chunk_id);
	static void UnloadQueryResource(MasterResources &resources, int chunk_id);
	static void LoadDatabaseInfo(DatabaseInfo &database_info,std::string database_info_filename);
	static void LoadDatabaseResource(WorkerResources &resources, int chunk_id);
	static void UnloadDatabaseResource(WorkerResources &resources , int chunk_id);
	
	static void loadFileData(std::string filename, char **ptr, uint64_t *size);
	
};




#endif
