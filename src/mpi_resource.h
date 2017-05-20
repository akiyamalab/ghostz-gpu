/* 
 * mpi_resource.h
 *    Created 2017/05/20
 *     Author:goto
 *
 */
#ifndef _MPI_RESOURCE_H_
#define _MPI_RESOURCE_H_

class MPIResource{
 public:

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

};




#endif
