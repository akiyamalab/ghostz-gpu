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


	if(resources.query_list[chunk_id].available){
		return 0;
	}
	int cmd[2];
	cmd[0]=MPIResource::RequestQuery;
	cmd[1]=chunk_id;
	
	MPI::COMM_WORLD.Send(cmd,2,MPI::INT,0,0);
	MPI::COMM_WORLD.Recv(cmd,2,MPI::INT,0,0);
	//cout<<cmd[0]<<":"<<cmd[1]<<endl;
	if(cmd[0]!=MPIResource::ACK){
		return 1;
	}
	//cout<<"recv query size:"<<cmd[1]<<endl;
	resources.query_list[chunk_id].size = cmd[1]; 
	resources.query_list[chunk_id].data = new char[cmd[1]];
	
	MPI::COMM_WORLD.Recv(resources.query_list[chunk_id].data,cmd[1],MPI::CHAR,0,0);

	
	
	
	
	return 0;	
}

int ResourceDeliverer::RequestDatabase(int chunk_id,WorkerResources &resources){
	if(resources.database_list[chunk_id].available){
		return 0;
	}
	int cmd[2];
	cmd[0]=MPIResource::RequestDatabase;
	cmd[1]=chunk_id;
	
	MPI::COMM_WORLD.Send(cmd,2,MPI::INT,0,0);
	MPI::COMM_WORLD.Recv(cmd,2,MPI::INT,0,0);
	if(cmd[0]!=MPIResource::ACK){
		return 1;
	}

	MPI::Request req[6];
	
	//req[0]=
		MPI::COMM_WORLD.Recv((char *)&resources.database_list[chunk_id].inf_size,
						  sizeof(uint64_t),MPI::CHAR,0,0);
		//req[1]=
		MPI::COMM_WORLD.Recv((char *)&resources.database_list[chunk_id].nam_size,
						  sizeof(uint64_t),MPI::CHAR,0,0);
		//req[2]=
		MPI::COMM_WORLD.Recv((char *)&resources.database_list[chunk_id].off_size,
						  sizeof(uint64_t),MPI::CHAR,0,0);
		//req[3]=
		MPI::COMM_WORLD.Recv((char *)&resources.database_list[chunk_id].seq_size,
						  sizeof(uint64_t),MPI::CHAR,0,0);
		//req[4]=
		MPI::COMM_WORLD.Recv((char *)&resources.database_list[chunk_id].scp_size,
						  sizeof(uint64_t),MPI::CHAR,0,0);
		//req[5]=
		MPI::COMM_WORLD.Recv((char *)&resources.database_list[chunk_id].sdp_size,
						  sizeof(uint64_t),MPI::CHAR,0,0);
		/*
	for(int i=0;i<6;i++){
		req[i].Wait();
	}*/
		//	resources.master_comm.Barrier();
	//cout<<"worker.requestdatabase ACK"<<endl;
	resources.database_list[chunk_id].inf = 
		new char[resources.database_list[chunk_id].inf_size];

	resources.database_list[chunk_id].nam = 
		new char[resources.database_list[chunk_id].nam_size];

	resources.database_list[chunk_id].off = 
		new char[resources.database_list[chunk_id].off_size];

	resources.database_list[chunk_id].seq = 
		new char[resources.database_list[chunk_id].seq_size];

	resources.database_list[chunk_id].scp = 
		new char[resources.database_list[chunk_id].scp_size];

	resources.database_list[chunk_id].sdp = 
		new char[resources.database_list[chunk_id].sdp_size];

	//req[0]=
		MPI::COMM_WORLD.Recv(resources.database_list[chunk_id].inf,
						  resources.database_list[chunk_id].inf_size,MPI::CHAR,0,0);
	//req[1]=
		MPI::COMM_WORLD.Recv(resources.database_list[chunk_id].nam,
						  resources.database_list[chunk_id].nam_size,MPI::CHAR,0,0);
	//req[2]=
		MPI::COMM_WORLD.Recv(resources.database_list[chunk_id].off,
						  resources.database_list[chunk_id].off_size,MPI::CHAR,0,0);
	//req[3]=
		MPI::COMM_WORLD.Recv(resources.database_list[chunk_id].seq,
						  resources.database_list[chunk_id].seq_size,MPI::CHAR,0,0);
	//req[4]=
		MPI::COMM_WORLD.Recv(resources.database_list[chunk_id].scp,
						  resources.database_list[chunk_id].scp_size,MPI::CHAR,0,0);
	//req[5]=
		MPI::COMM_WORLD.Recv(resources.database_list[chunk_id].sdp,
						  resources.database_list[chunk_id].sdp_size,MPI::CHAR,0,0);
		/*
	for(int i=0;i<6;i++){
		req[i].Wait();
		}*/
	
	return 0;
}

int ResourceDeliverer::RequestResource(char type, int chunk_id,WorkerResources &resources){
	return 0;
}



void ResourceDeliverer::CreateResponseThread(MasterResources &resources,int size){
	ThreadParameter param;
	param.resources=resources;
	//	boost::thread_group threads
	
	for(int i=1;i<size;i++){
		threads.create_thread(boost::bind(&ResponseThread, i, param,resources));
	}

	
}
void ResourceDeliverer::DestroyResponseThread(){
	threads.interrupt_all();
}

void ResourceDeliverer::ResponseThread(int thread_id,ThreadParameter &parameter,MasterResources &resources){
	int cmd[2];
	int target;
	MPI::COMM_WORLD.Recv(cmd,2,MPI::INT,thread_id,0);
	
	switch(cmd[0]){
	case MPIResource::RequestQuery:
		target=cmd[1];
		cout<<"target"<<target<<":"<<resources.query_list[target].data<<endl;
		if(resources.query_list[target].available==false){// target query not loaded
			cmd[0]=MPIResource::NACK;
			cout<<resources.query_list[target].available<<"NACK"<<endl;
			MPI::COMM_WORLD.Send(cmd,2,MPI::INT,thread_id,0);
			break;
		}

		cmd[0]=MPIResource::ACK;
		cmd[1]=resources.query_list[target].size;
		
		MPI::COMM_WORLD.Send(cmd,2,MPI::INT,thread_id,0);
		MPI::COMM_WORLD.Send(resources.query_list[target].data,cmd[1],MPI::CHAR,thread_id,0);
	  
		break;

	case MPIResource::RequestDatabase:
	  
		break;

	case MPIResource::RequestTask:
	
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
