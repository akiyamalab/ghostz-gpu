/*
 * resource_deliverer.h
 *    Created on 2017/05/20
 *       Author:goto
 * 
 */

#ifndef _RESOURCE_DELIVERER_H_
#define _RESOURCE_DELIVERER_H_

#include "mpi_common.h"

#include "mpi.h"
#include <iostream>
#include <boost/thread.hpp>


class ResourceDeliverer {
 public :
	typedef MPICommon::QueryResource QueryResource;
	typedef MPICommon::MasterResources MasterResources;
	typedef MPICommon::WorkerResources WorkerResources;
	


	//worker method
    int RequestQuery(int chunk_id,WorkerResources &resources);
	int RequestDatabase(int chunk_id,WorkerResources &resorurces);
	
	
	//master method
	void CreateResponseThread(MasterResources &resources,int size);
	
 private :
	struct ThreadParameter{
		MasterResources resources;
		boost::mutex *load;
		
	};

	enum Command{
		CMD_RequestQuery,
		CMD_RequestDatabase,
		CMD_RequestTask,
		ACK,
		NACK
	};

	//worker
	int RequestResource(char type, int chunk_id,WorkerResources &resources);
	
	int RecvQuery(QueryResource &query);

	
	static void ResponseThread(int thread_id, ResourceDeliverer::ThreadParameter &parameter);
	int SendQuery(int target,QueryResource &query);
	
  

};
#endif
