/*
 * load_balancer.cpp
 *   Created on :2017/06/20
 *     Author:goto
 *
 */

#include "load_balancer.h"
using namespace std;
LoadBalancer::LoadBalancer(MasterResources &resources,int num_process):
	master_(resources){
	db_chunk_ = master_.database_info.number_chunks;
	queues.resize(db_chunk_);
	for(int i=0;i<queues.size();i++){
		for(int j=0;j<resources.query_list.size();j++){
			AlignmentTask task;
			task.database_chunk=i;
			task.query_chunk=j;
			queues[i].push(task);
		}
	}
	db_map.resize(num_process);
	
}

void LoadBalancer::SetDatabaseLoadingMap(int worker_rank,int target_chunk){
	db_map[worker_rank]=target_chunk;
}

int LoadBalancer::GetSwitchTargetChunk(int worker_rank){
	vector<int> num_load(queues.size(),0);
	for(int i=0;i<db_map.size();i++){
		num_load[db_map[i]]++;
	}
	float task_ratio,max_task_ratio;
	max_task_ratio=0;
   
	int switch_target=0;
	
	
	for(int i=0;i<queues.size();i++){
		if(db_map[i]==0 && queues[i].size()>0){
			return i;// no loaded db is alway target.
		}
		task_ratio=queues[i].size()/(float)db_map[i];
		if(max_task_ratio < task_ratio){
			max_task_ratio=task_ratio;
			switch_target = i;
		}			
	} 
	return switch_target;
}

LoadBalancer::AlignmentTask LoadBalancer::getNext(int target_chunk){
	AlignmentTask task;
	if(queues[target_chunk].empty()){
		task.query_chunk=-1;
		task.database_chunk=-1;
		return task;
	}
		
	task=queues[target_chunk].front();  
	queues[target_chunk].pop();
	return task;
}
