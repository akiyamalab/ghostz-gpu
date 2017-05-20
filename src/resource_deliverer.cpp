/*
 *  resource_deliverer.cpp
 *     Created on 2017/05/20
 *       Author:goto
 *
 */

#include "resource_deliverer.h"

#include "mpi.h"
#include "mpi_common.h"
#include <iostream>
#include <boost/thread.hpp>

using namespace std;
int ResourceDeliverer::RequestQuery(int chunk_id,WorkerResources &resources){
#if 1
	cout<<"send request"<<endl;
#endif

	if(resources.query_list[chunk_id].available){
		return 0;
	}
	int cmd[2];
	cmd[0]=CMD_RequestQuery;
	cmd[1]=chunk_id;
	
	MPI::COMM_WORLD.Send(cmd,2,MPI::INT,0,0);
	MPI::COMM_WORLD.Recv(cmd,2,MPI::INT,0,0);
	if(cmd[0]!=ACK){
		return 1;
	}
	resources.query_list[chunk_id].data = new char[cmd[1]];
	MPI::COMM_WORLD.Recv(resources.query_list[chunk_id].data,cmd[1],MPI::CHAR,0,0);

	
	
	
	
	return 0;	
}

int ResourceDeliverer::RequestDatabase(int chunk_id,WorkerResources &resources){
	return 0;
}

int ResourceDeliverer::RequestResource(char type, int chunk_id,WorkerResources &resources){
	return 0;
}



void ResourceDeliverer::CreateResponseThread(MasterResources &resources,int size){
	ThreadParameter param;
	param.resources=resources;
	//	boost::thread_group threads;
	
	for(int i=1;i<size;i++){
		threads.create_thread(boost::bind(&ResponseThread, i, param));
	}

	
}
void ResourceDeliverer::DestoryResponseThread(){
	threads.interrupt_all();
}

void ResourceDeliverer::ResponseThread(int thread_id,ThreadParameter &parameter){
	int cmd[2];
	int target;
	MPI::COMM_WORLD.Recv(cmd,2,MPI::INT,thread_id,0);
	
	switch(cmd[0]){
	case CMD_RequestQuery:
		target=cmd[1];
		if(!parameter.resources.query_list[target].available){// target query not loaded
			cmd[0]=NACK;
			MPI::COMM_WORLD.Send(cmd,2,MPI::INT,thread_id,0);
			break;
		}

		cmd[0]=ACK;
		cmd[1]=parameter.resources.query_list[target].size;
		
		MPI::COMM_WORLD.Send(cmd,2,MPI::INT,thread_id,0);
		MPI::COMM_WORLD.Send(parameter.resources.query_list[target].data,cmd[1],MPI::CHAR,thread_id,0);
	  
		break;

	case CMD_RequestDatabase:
	  
		break;

	case CMD_RequestTask:
	
		break;
	default:
		
		break;	
		
	}
	
	
}
int ResourceDeliverer::RecvQuery(QueryResource &query){
	int size;
	MPI::COMM_WORLD.Recv(&size,1,MPI::INT,0,0);
	query.data = new char[size];
	MPI::COMM_WORLD.Recv(query.data,size,MPI::CHAR,0,0);
	query.size=size;
	query.available=true;	
	return 0;
}

int ResourceDeliverer::SendQuery(int target,QueryResource &query){
	MPI::COMM_WORLD.Send(&(query.size),1,MPI::INT,target,0);
	MPI::COMM_WORLD.Send(query.data,query.size,MPI::CHAR,target,0);
	
	return 0;
}
