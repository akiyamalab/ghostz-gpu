/* 
 *  mpi_resource.cpp
 *     Created on 2017/06/17
 *      Author:goto
 */



#include "mpi_resource.h"


using namespace std;


int MPIResource::BcastDatabase(DatabaseResource &database,MPI::Intercomm comm,int root){

	comm.Bcast((char *)&database.inf_size,sizeof(uint64_t),MPI::CHAR,0);
	comm.Bcast((char *)&database.nam_size,sizeof(uint64_t),MPI::CHAR,0);
	comm.Bcast((char *)&database.off_size,sizeof(uint64_t),MPI::CHAR,0);
	comm.Bcast((char *)&database.seq_size,sizeof(uint64_t),MPI::CHAR,0);
	comm.Bcast((char *)&database.scp_size,sizeof(uint64_t),MPI::CHAR,0);
	comm.Bcast((char *)&database.sdp_size,sizeof(uint64_t),MPI::CHAR,0);
	if(comm.Get_rank()!=0){
		database.inf = new char[database.inf_size];
		database.nam = new char[database.nam_size];
		database.off = new char[database.off_size];
		database.seq = new char[database.seq_size];
		database.scp = new char[database.scp_size];
		database.sdp = new char[database.sdp_size];
	}
	comm.Bcast(database.inf,database.inf_size,MPI::CHAR,0);
	comm.Bcast(database.nam,database.nam_size,MPI::CHAR,0); 
	comm.Bcast(database.off,database.off_size,MPI::CHAR,0);
	comm.Bcast(database.seq,database.seq_size,MPI::CHAR,0);
	comm.Bcast(database.scp,database.scp_size,MPI::CHAR,0);
	comm.Bcast(database.sdp,database.sdp_size,MPI::CHAR,0);
	
	return 0;
}
int MPIResource::AcceptCommand(MasterResources &resources){
	int cmd[2];
	MPI::Status status;
	MPI::COMM_WORLD.Recv(&cmd,sizeof(cmd),MPI::INT,MPI::ANY_SOURCE,MPI::ANY_TAG,status);
	
	int dst_rank = status.Get_source();
	switch(cmd[0]){
    case MPIResource::CMD_RequestQuery :
		
        AcceptRequestQuery(resources,cmd[1],dst_rank);

        break;
    case MPIResource::CMD_RequestDatabase :
        //AcceptRequestDatabase(cmd,resources,status);
        break;
    case MPIResource::CMD_RequestTask :
        //AcceptRequestTask(cmd,resources,status);
        break;
    default :
        break;
    }

	
	
}

int MPIResource::RequestQuery(WorkerResources &resources, int target_chunk){
	if(resources.query_list[target_chunk].available){
		return 0;
	}
	int cmd[2];
	cmd[0]=MPIResource::CMD_RequestQuery;
	cmd[1]=target_chunk;
	
	MPI::COMM_WORLD.Send(cmd,sizeof(cmd),MPI::INT,0,0);
	return RecvQuery(resources.query_list[target_chunk],MPI::COMM_WORLD,0);
	
}

int MPIResource::AcceptRequestQuery(MasterResources &resources, int target_chunk, int dst_rank){
	return SendQuery(resources.query_list[target_chunk],MPI::COMM_WORLD,dst_rank);
}

int MPIResource::RecvQuery(QueryResource &query_resource, MPI::Intercomm  comm, int src_rank){
	comm.Recv((char *) &(query_resource.size),sizeof(query_resource.size),MPI::CHAR,src_rank,0);
	query_resource.data = new char[query_resource.size];
	comm.Recv(query_resource.data,query_resource.size,MPI::CHAR,src_rank,0);
	return 0;
}

int MPIResource::SendQuery(QueryResource &query_resource, MPI::Intercomm comm, int dst_rank){
	comm.Send((char *) &(query_resource.size),sizeof(query_resource.size),MPI::CHAR,dst_rank,0);
	comm.Send(query_resource.data,query_resource.size,MPI::CHAR,dst_rank,0);
	return 0;
}
