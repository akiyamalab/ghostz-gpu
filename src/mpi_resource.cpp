/* 
 *  mpi_resource.cpp
 *     Created on 2017/06/17
 *      Author:goto
 */



#include "mpi_resource.h"
#include <fstream>
#include <sstream>
using namespace std;


int MPIResource::BcastDatabase(DatabaseResource &database,MPI::Intercomm comm,int root){
	comm.Bcast((char *)&database.inf_size,sizeof(uint64_t),MPI::CHAR,root);
	comm.Bcast((char *)&database.nam_size,sizeof(uint64_t),MPI::CHAR,root);
	comm.Bcast((char *)&database.off_size,sizeof(uint64_t),MPI::CHAR,root);
	comm.Bcast((char *)&database.seq_size,sizeof(uint64_t),MPI::CHAR,root);
	comm.Bcast((char *)&database.scp_size,sizeof(uint64_t),MPI::CHAR,root);
	comm.Bcast((char *)&database.sdp_size,sizeof(uint64_t),MPI::CHAR,root);
	if(comm.Get_rank()!=0){
		database.inf = new char[database.inf_size];
		database.nam = new char[database.nam_size];
		database.off = new char[database.off_size];
		database.seq = new char[database.seq_size];
		database.scp = new char[database.scp_size];
		database.sdp = new char[database.sdp_size];
	}

	BcastLargeData(&database.inf,database.inf_size,comm,root);
	BcastLargeData(&database.nam,database.nam_size,comm,root);
	BcastLargeData(&database.off,database.off_size,comm,root);
	BcastLargeData(&database.seq,database.seq_size,comm,root);
	BcastLargeData(&database.scp,database.scp_size,comm,root);
	BcastLargeData(&database.sdp,database.sdp_size,comm,root);

	/*
	comm.Bcast(database.inf,database.inf_size,MPI::CHAR,0);
	comm.Bcast(database.nam,database.nam_size,MPI::CHAR,0); 
	comm.Bcast(database.off,database.off_size,MPI::CHAR,0);
	comm.Bcast(database.seq,database.seq_size,MPI::CHAR,0);
	comm.Bcast(database.scp,database.scp_size,MPI::CHAR,0);
	comm.Bcast(database.sdp,database.sdp_size,MPI::CHAR,0);
	*/
	return 0;
}

int MPIResource::BcastDatabaseInfo(DatabaseInfo &info,MPI::Intercomm comm,int root){
	comm.Bcast((char *)&info.number_chunks,
			   sizeof(info.number_chunks),MPI::CHAR,root);
	comm.Bcast((char *)&info.max_sequence_length,
			   sizeof(info.max_sequence_length),MPI::CHAR,root);
	comm.Bcast((char *)&info.database_length,
			   sizeof(info.database_length),MPI::CHAR,root);
	comm.Bcast((char *)&info.number_sequences,
			   sizeof(info.number_sequences),MPI::CHAR,root);
	comm.Bcast((char *)&info.sequence_delimiter,
			   sizeof(info.sequence_delimiter),MPI::CHAR,root);
	
	
}

void MPIResource::BcastLargeData(char **ptr,uint64_t size,MPI::Intercomm comm,int root){
    char *ptr_ = *ptr;
    int iter=size/(1024*1024*1024);
    for(int i=0;i<iter;i++){
        comm.Bcast(ptr_,1024*1024*1024,MPI::CHAR,root);
        ptr_+= 1024*1024*1024;
        size-= 1024*1024*1024;
    }
    if(size>0){
        comm.Bcast(ptr_,size,MPI::CHAR,root);
    }
}

int MPIResource::AcceptCommand(MasterResources &resources, int *cmd,int *target){
	int buf[2];
	MPI::Status status;
	MPI::COMM_WORLD.Recv(buf,sizeof(buf),MPI::INT,MPI::ANY_SOURCE,MPI::ANY_TAG,status);
	
	int dst_rank = status.Get_source();
	*cmd=buf[0];
	*target=buf[1];
	return dst_rank;
	
	
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

void MPIResource::RequestTask(int target_chunk,MPIResource::AlignmentTask &task){
	int cmd[2];
	int cmd_[2];
	cmd[0]=MPIResource::CMD_RequestTask;
	cmd[1]=target_chunk;
	
	MPI::COMM_WORLD.Send(cmd,sizeof(cmd),MPI::INT,0,0);
	
	MPI::COMM_WORLD.Recv(cmd_,sizeof(cmd_),MPI::INT,0,0);
	//cout<<"recv "<<cmd_[0]<<endl;	
	task.query_chunk=cmd_[0];
	task.database_chunk=cmd_[1];
	//*task_query=0;
	//*task_db=0; 
	//cout<<"recv "<<task.query_chunk<<":"<<task.database_chunk<<endl;
	//return task;
	return ;
}

int MPIResource::AcceptRequestTask(AlignmentTask task, int dst_rank){
	int cmd[2];
	cmd[0]=task.query_chunk;
	cmd[1]=task.database_chunk;
	MPI::COMM_WORLD.Send(cmd,sizeof(cmd),MPI::INT,dst_rank,0);
	return 0;
	
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


void MPIResource::LoadQueryResource(MasterResources &resources,int chunk_id){
	
    if(chunk_id > resources.query_list.size()){
        cerr<<"query chunk id error:"<<chunk_id<<endl;
        MPI_Abort(MPI_COMM_WORLD,1);
    }
    if(resources.query_list[chunk_id].available){
        return ;
    }
    resources.query_list[chunk_id].data
        = new char[resources.query_list[chunk_id].size];

    ifstream in(resources.query_filename.c_str());
    in.seekg(resources.query_list[chunk_id].ptr);
    in.read(resources.query_list[chunk_id].data,resources.query_list[chunk_id].size);
    in.close();
    resources.query_list[chunk_id].available=true;
    //cout<<"query_chunk:"<<chunk_id<<"loaded."<<endl;
}

void MPIResource::UnloadQueryResource(MasterResources &resources,int chunk_id){
	if(chunk_id > resources.query_list.size()){
        cerr<<"query chunk id error:"<<chunk_id<<endl;
        MPI_Abort(MPI_COMM_WORLD,1);
    }
    if(!resources.query_list[chunk_id].available){
        return ;
    }

    delete [] resources.query_list[chunk_id].data;
    resources.query_list[chunk_id].available=false;

}
void MPIResource::UnloadQueryResource(WorkerResources &resources,int chunk_id){
	if(chunk_id > resources.query_list.size()){
        cerr<<"query chunk id error:"<<chunk_id<<endl;
        MPI_Abort(MPI_COMM_WORLD,1);
    }
    if(!resources.query_list[chunk_id].available){
        return ;
    }

    delete [] resources.query_list[chunk_id].data;
    resources.query_list[chunk_id].available=false;
	
}

void MPIResource::LoadDatabaseInfo(DatabaseInfo &database_info,std::string database_info_filename){
	ifstream in;
	in.open(database_info_filename.c_str(),ios::binary);	
    in.read((char *)&database_info.number_chunks,
			sizeof(database_info.number_chunks));
    in.read((char *)&database_info.max_sequence_length,
			sizeof(database_info.max_sequence_length));
    in.read((char *)&database_info.database_length,
			sizeof(database_info.database_length));
    in.read((char *)&database_info.number_sequences,
			sizeof(database_info.number_sequences));
    in.read((char *)&database_info.sequence_delimiter,
			sizeof(database_info.sequence_delimiter));
    in.close();
}
void MPIResource::LoadDatabaseResource(WorkerResources &resources,int chunk_id){
    if(chunk_id > resources.database_list.size()){
        cerr<<"database chunk id error:"<<chunk_id<<endl;
        MPI_Abort(MPI_COMM_WORLD,1);
    }
    if(resources.database_list[chunk_id].available){
        cout<<"already loaded.:"<<chunk_id<<endl;
		return ;
    }
    stringstream ss;
    //.inf
    ss<<resources.database_filename<<"_"<<chunk_id<<".inf";
    loadFileData(ss.str().c_str(),
                 &(resources.database_list[chunk_id].inf),
                 &(resources.database_list[chunk_id].inf_size));

    //.nam
    ss.str("");
    ss<<resources.database_filename<<"_"<<chunk_id<<".nam";
    loadFileData(ss.str().c_str(),
                 &(resources.database_list[chunk_id].nam),
                 &(resources.database_list[chunk_id].nam_size));
    //.off
    ss.str("");
    ss<<resources.database_filename<<"_"<<chunk_id<<".off";
    loadFileData(ss.str().c_str(),
                 &(resources.database_list[chunk_id].off),
                 &(resources.database_list[chunk_id].off_size));
    //.seq
    ss.str("");
    ss<<resources.database_filename<<"_"<<chunk_id<<".seq";
    loadFileData(ss.str().c_str(),
                 &(resources.database_list[chunk_id].seq),
                 &(resources.database_list[chunk_id].seq_size));
    //.scp
    ss.str("");
    ss<<resources.database_filename<<"_"<<chunk_id<<".scp";
    loadFileData(ss.str().c_str(),
                 &(resources.database_list[chunk_id].scp),
                 &(resources.database_list[chunk_id].scp_size));
    //.sdp
    ss.str("");
    ss<<resources.database_filename<<"_"<<chunk_id<<".sdp";
    loadFileData(ss.str().c_str(),
                 &(resources.database_list[chunk_id].sdp),
                 &(resources.database_list[chunk_id].sdp_size));
    resources.database_list[chunk_id].available=true;
}
void MPIResource::UnloadDatabaseResource(WorkerResources &resources,int chunk_id){
    if(chunk_id > resources.database_list.size()){
        cerr<<"database chunk id error:"<<chunk_id<<endl;
        MPI_Abort(MPI_COMM_WORLD,1);
    }
    if(!resources.database_list[chunk_id].available){
        return ;
	}
	//	cout<<"unloading.:"<<chunk_id<<endl;
    delete [] resources.database_list[chunk_id].inf;
    delete [] resources.database_list[chunk_id].nam;
    delete [] resources.database_list[chunk_id].off;
    delete [] resources.database_list[chunk_id].seq;
    delete [] resources.database_list[chunk_id].scp;
    delete [] resources.database_list[chunk_id].sdp;

    resources.database_list[chunk_id].available=false;

}

void MPIResource::loadFileData(std::string filename,char **ptr,uint64_t *size){
    ifstream in(filename.c_str());
    uint64_t begin,end;
    uint64_t size_;
    in.seekg(0,ios::end);
    end = in.tellg();
    in.clear();
    in.seekg(0,ios::beg);
    begin=in.tellg();
    size_=end-begin;

    *ptr = new char[size_];
	// cout<<filename.c_str()<<" size:"<<size_<<endl;
    in.read(*ptr,size_);
    in.close();
    *size=size_;

}
