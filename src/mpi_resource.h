/* 
 * mpi_resource.h
 *    Created 2017/05/20
 *     Author:goto
 *
 */
#ifndef _MPI_RESOURCE_H_
#define _MPI_RESOURCE_H_

#include <queue>

class MPIResource{
 public:
	
	struct AlignmentTask{
		int query_chunk;
		int database_chunk;
	};
	
	struct QueryResource{
        char *data;
        int ptr;
        int size;
        int chunk_id;
        bool available;
    };
	

    struct MasterResources{ 
		std::string query_filename;
		std::vector<QueryResource> query_list;             //size= number of query chunk
		std::vector<std::queue<AlignmentTask> > task_pool; //size= number of database chunk 
		std::vector<int> node_loading_database;            //size= number of (master+)worker node
		std::vector<float> est_speed;                      //size= number of database chunk
		
		
    };
    struct WorkerResources{
		std::vector<QueryResource> query_list;
    };
	

	enum Command{
        RequestQuery,
        RequestDatabase,
        RequestTask,
        ACK,
        NACK
    };

};




#endif
