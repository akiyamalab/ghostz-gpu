/*
 *
 * load_balancer.h
 *   Created on :2017/06/20
 *     Author:goto
 *
 */
#ifndef __LOAD_BALANCER_H__
#define __LOAD_BALANCER_H__
#include "mpi_resource.h"
#include <queue>
#include <vector>
#include <sstream>

class LoadBalancer{
 public:
	typedef MPIResource::AlignmentTask AlignmentTask;
	typedef MPIResource::MasterResources MasterResources;
	LoadBalancer(){
	};
	LoadBalancer(int n_query_chunk,int n_db_chunk,int n_process);
	int GetSwitchTargetChunk();
	void SetDatabaseLoadingMap(int worker_rank,int target_chunk);
	AlignmentTask GetNext(int target_chunk);
	std::string print();
	int GetRemainTask();
 private:
	int n_db_chunk_;
	int n_query_chunk_;
	int n_process_;
	std::vector<std::queue<AlignmentTask> > queues;
	MasterResources master_;
	std::vector<int> db_map;
	

} ;

#endif /* __LOAD_BALANCER_H__*/
