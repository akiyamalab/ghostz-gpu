/*
 * load_balancer.cpp
 *   Created on :2017/06/20
 *     Author:goto
 *
 */

#include "load_balancer.h"
using namespace std;

LoadBalancer::LoadBalancer(int n_query_chunk,int n_db_chunk,int n_process):
	n_query_chunk_(n_query_chunk),n_db_chunk_(n_db_chunk),n_process_(n_process){
	queues.resize(n_db_chunk_);
	for(int i=0;i<n_db_chunk_;i++){
		for(int j=0;j<n_query_chunk_;j++){
			AlignmentTask task;
			task.database_chunk=i;
			task.query_chunk=j;
			queues[i].push(task);
		}
	}
	db_map.resize(n_process_);
	
}

void LoadBalancer::SetDatabaseLoadingMap(int worker_rank,int target_chunk){
	db_map[worker_rank]=target_chunk;
}

int LoadBalancer::GetSwitchTargetChunk(){
	vector<int> num_load(queues.size(),0);
	for(int i=0;i<db_map.size();i++){
		num_load[db_map[i]]++;
	}
	float task_ratio,max_task_ratio;
	max_task_ratio=0;
   
	int switch_target=0;
	
	
	for(int i=0;i<queues.size();i++){
		if(num_load[i]==0 && queues[i].size()>0){
			return i;// no loaded db is alway target.
		}
		task_ratio=queues[i].size()/(float)num_load[i];
		if(max_task_ratio < task_ratio){
			max_task_ratio=task_ratio;
			switch_target = i;
		}			
	} 
	return switch_target;
}

LoadBalancer::AlignmentTask LoadBalancer::GetNext(int target_chunk){
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

int LoadBalancer::GetRemainTask(){
	int count=0;
	for(int i=0;i<queues.size();i++){
		count+=queues[i].size();
	}
	return count;
}
std::string LoadBalancer::print(){
	stringstream ss;
	ss.str("");
	ss<<"db_chk\ttasks\trank\n";
	for(int i=0;i<n_db_chunk_;i++){
		ss<<i<<"\t"<<queues[i].size()<<"\t";
		for(int j=1;j<n_process_;j++){
			if(db_map[j]==i){
				ss<<j<<" ";
			}
		}
		ss<<"\n";
	}
	return ss.str();
}
