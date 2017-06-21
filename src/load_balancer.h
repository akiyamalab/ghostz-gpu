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

class LoadBalancer{
 public:
	typedef MPIResource::AlignmentTask AlignmentTask;
	typedef MPIResource::MasterResources MasterResources;
	
	LoadBalancer(MasterResources &resources,int num_process);
	int GetSwitchTargetChunk(int worker_rank);
	void SetDatabaseLoadingMap(int worker_rank,int target_chunk);
	AlignmentTask getNext(int target_chunk);
	
	
 private:
	int db_chunk_;
	std::vector<std::queue<AlignmentTask> > queues;
	MasterResources master_;
	std::vector<int> db_map;
	

} ;

#endif /* __LOAD_BALANCER_H__*/
