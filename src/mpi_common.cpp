/*
 * mpi_common.cpp
 *  Created on 2017/05/16
 *  Author:goto
 */

#include "mpi_common.h"
#include "aligner_mpi.h"
#include "fasta_sequence_reader.h"
#include "protein_type.h"
#include "dna_type.h"
#include "sequence_type.h"
#include "queries.h"
#include "logger.h"
#include "score_matrix_reader.h"
#include "reduced_alphabet_file_reader.h"
#include "mpi_resource.h"
#include  "aligner_mpi.h"
#include "aligner.h"

#include "mpi.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;


void MPICommon::debug(int argc,char *argv[]){
	
	if(argc>=2){
		//cout<<argv[1]<<endl;
	}
	string input_file;
	string database_file;
	string output_file;
	AligningParameters parameter;
	BuildParameters(argc,argv,input_file,database_file,output_file,
					parameter);
	
	
	
	ifstream isf(input_file.c_str());
	
	ofstream osf("1chunk.fasta");
	Queries::Parameters queries_parameters;
	
	AlignerCommon::BuildQueriesParameters(parameter,queries_parameters);
	Queries queries(queries_parameters);
	int i=0;
	bool eof_flag=false;
	while(!eof_flag){
		ifstream::pos_type begin = isf.tellg();
		ifstream::pos_type isf_ptr=queries.GetNextChunkPosition(isf);
		if(isf_ptr==-1){
			isf.clear();
			isf.seekg(begin);
			isf_ptr=isf.seekg(0,std::ios_base::end).tellg();
			eof_flag=true;
		   
		}
		
		cout<<"begin,isf_ptr:"<<begin<<","<<isf_ptr<<endl;
		stringstream ss;
		ss<<output_file.c_str()<<"_"<<i;
		ofstream ofs(ss.str().c_str(),ios::binary | ios::out);

		char *array = new char[isf_ptr-begin];
		isf.clear();

		isf.seekg(begin);
		isf.read(array,isf_ptr-begin);
		//cout<<array<<endl;
		ofs.write(array,isf_ptr-begin);
		ofs.flush();
		ofs.close();
		isf.seekg(isf_ptr);
		i++;
		delete [] array;
	}
	
	//cout<<name<<"\n"<<sequence<<endl;
	   uint32_t length;
	
	
	std::tr1::shared_ptr<SequenceType> queries_sequence_type;
	std::tr1::shared_ptr<SequenceType> database_sequence_type;
	queries_sequence_type  = std::tr1::shared_ptr<SequenceType>(new DnaType());
	database_sequence_type = std::tr1::shared_ptr<SequenceType>(new ProteinType());

	//length =Queries::GetSequenceLength(queries_sequence_type,database_sequence_type,
	//							   sequence.size());
	//cout<<length<<endl;
	
	
							   
								   

}
void MPICommon::Run(string &queries_filename,string &database_filename,
					string &output_filename,
					AligningParameters &parameter,MPIParameter &mpi_parameter){
	
	if(mpi_parameter.rank==0){
		RunMaster(queries_filename,database_filename,output_filename,
				  parameter,mpi_parameter);
	}else{
		RunWorker(parameter,mpi_parameter);
		}
}

void MPICommon::RunGPU(string &queries_filename,string &database_filename,
					string &output_filename,
					AligningParameters &parameter,MPIParameter &mpi_parameter){
	
	if(mpi_parameter.rank==0){
		RunMaster(queries_filename,database_filename,output_filename,
				  parameter,mpi_parameter);
	}else{
		RunWorkerGPU(parameter,mpi_parameter);
	}
}



void MPICommon::RunMaster(string &queries_filename,string &database_filename,
						  string &output_filename,
						  AligningParameters &parameter,MPIParameter &mpi_parameter){
	
	//Master Process Init
	//cout<<"master start:"<<mpi_parameter.rank<<endl;
	MasterResources resources;
	resources.query_filename=queries_filename;
	resources.database_filename=database_filename;
	resources.output_filename;

	
	
	
	//***********************//
	//Start Init Phase

	
	SetupMasterResources(queries_filename,database_filename,resources,parameter,mpi_parameter);
#if 0
	
	for(int i=0;i<resources.query_list.size();i++){
		cout<<resources.query_list[i].chunk_id;
		cout<<" ptr:"<<resources.query_list[i].ptr;
		cout<<" size:"<<resources.query_list[i].size<<endl;
	}
#endif
	
	//Create Subgroup comm
	for(int i=0;i<resources.database_list.size();i++){
		MPI::COMM_WORLD.Split(1,0);
	}
	
	for(int i=0;i<resources.query_list.size();i++){
		LoadQueryResource(resources,i);
		
	}   //loading all query on master memory
	
	//	cout<<resources.query_list[0].data;
	
	

	//End Init Phase
	//***********************//
	//Start Search Phase 
		for(int i=0;i<mpi_parameter.size-1;i++){
			MPIResource::AcceptCommand(resources);
		}	
		//End Search Phase
	//***********************//
	//Start Report Phase
	
	//finalize
	MPI::COMM_WORLD.Barrier();	
	
}
void MPICommon::RunWorker(AligningParameters &parameter,MPIParameter &mpi_parameter){
	
	//Worker Process Init
	//cout<<"worker start:"<<mpi_parameter.rank<<endl;
	WorkerResources resources;
	//***********************//
	//Start Init Phase
	
	SetupWorkerResources(resources,mpi_parameter);
	//Create subgroup comm
	int rank=mpi_parameter.rank;
	int db_chunk_size = resources.database_list.size();
	bool submaster = false;
	int target_chunk=(rank-1)%db_chunk_size;
	if(target_chunk==rank-1){
		submaster =true;
	}
	
	for(int i=0;i<db_chunk_size;i++){
		if(i==target_chunk){
			if(submaster){
				resources.subgroup_comm = MPI::COMM_WORLD.Split(0,0);
			}else{
				resources.subgroup_comm = MPI::COMM_WORLD.Split(0,rank);
			}
		}else{
			MPI::COMM_WORLD.Split(1,rank);
		}
	}
	
	if(submaster){
		//load db from filesystem
		LoadDatabaseResource(resources,target_chunk);

	}
	//Broadcast db to subgroup
	uint64_t sum =0;
	MPIResource::BcastDatabase(resources.database_list[target_chunk],resources.subgroup_comm,0);


	DatabaseType database(resources.database_list[target_chunk],resources.database_info);
	DatabaseType database_(resources.database_filename);
	//cout<<"database.GetChunkId = "<<database.GetChunkId()<<endl;
	
	//End Init Phase
	//***********************//
	//Start Search Phase
	AlignmentTask task;
	task.query_chunk=0;
	task.database_chunk=0;
	
	
	if(!resources.query_list[task.query_chunk].available){
		MPIResource::RequestQuery(resources,task.query_chunk);
	}

	
	
	AlignerMPI aligner;
	Aligner aligner_;
	vector<vector<Result> > results_list;
	string q_name=string("/work1/t2ggenome/goto/work/ghostz/mpi/test/query/ERR315856.fasta_randN1k");
	string out_name = string("out");
	aligner.Search(resources.query_list[task.query_chunk],database,results_list,parameter,mpi_parameter);
	aligner_.Align(q_name,  resources.database_filename,out_name,parameter);
	for(int i=0;i<results_list.size();i++){
		for(int j=0;j<results_list[i].size();j++){
			//cout<<results_list[i][j].subject_name<<endl;
		}}
	//End Search Phase
	//***********************//
	//Start Report Phase

	//End Report Phase

	
	//finalize

	MPI::COMM_WORLD.Barrier();
}


void MPICommon::RunWorkerGPU(AligningParameters &parameter,MPIParameter &mpi_parameter){

}

void MPICommon::SetupMasterResources(string &queries_filename,string &database_filename,
									 MasterResources &resources,AligningParameters &parameter,
									 MPIParameter &mpi_parameter){
	
	//init communicator;
	
	/*
	resources.comm_list.push_back(MPI::COMM_SELF);
	for(int i=1;i<mpi_parameter.size;i++){
		resources.comm_list.push_back(MPI::COMM_WORLD.Split(0,0));
		cout<<"comm size"<<resources.comm_list[i].Get_size()<<endl;
		}*/

	

	//init query resource
	vector<int> pointer_list;
	vector<int> size_list;
	
	BuildQueryChunkPointers(queries_filename,pointer_list,size_list,parameter);
	int query_chunk_number = pointer_list.size();
	for(int i=0;i<query_chunk_number;i++){
		QueryResource query_resource;
		query_resource.size=size_list[i];
		query_resource.ptr=pointer_list[i];
		query_resource.chunk_id=i;
		query_resource.available=false;
		resources.query_list.push_back(query_resource);
	}


	//init database resource
	string inf_file = resources.database_filename+".inf";
	DatabaseInfo databaseinfo;
	/*
	int number_chunks;
	uint32_t max_sequence_length;
	uint64_t database_length;
	uint64_t number_sequences;
	AlphabetCoder::Code sequence_delimiter;
	*/
	
	ifstream in;
	int dbinfo_size[5]; 
	dbinfo_size[0]=sizeof(databaseinfo.number_chunks);
	dbinfo_size[1]=sizeof(databaseinfo.max_sequence_length);	
	dbinfo_size[2]=sizeof(databaseinfo.database_length);
	dbinfo_size[3]=sizeof(databaseinfo.number_sequences);
	dbinfo_size[4]=sizeof(databaseinfo.sequence_delimiter);

	in.open(inf_file.c_str(),ios::binary);
	in.read((char *)&databaseinfo.number_chunks,dbinfo_size[0]);
	in.read((char *)&databaseinfo.max_sequence_length,dbinfo_size[1]);
	in.read((char *)&databaseinfo.database_length,dbinfo_size[2]);
	in.read((char *)&databaseinfo.number_sequences,dbinfo_size[3]);
	in.read((char *)&databaseinfo.sequence_delimiter,dbinfo_size[4]);
	in.close();
	
	MPI::COMM_WORLD.Bcast((char *)&databaseinfo.number_chunks,dbinfo_size[0],MPI::CHAR,0);
	MPI::COMM_WORLD.Bcast((char *)&databaseinfo.max_sequence_length,dbinfo_size[1],MPI::CHAR,0);
	MPI::COMM_WORLD.Bcast((char *)&databaseinfo.database_length,dbinfo_size[2],MPI::CHAR,0);
	MPI::COMM_WORLD.Bcast((char *)&databaseinfo.number_sequences,dbinfo_size[3],MPI::CHAR,0);
	MPI::COMM_WORLD.Bcast((char *)&databaseinfo.sequence_delimiter,dbinfo_size[4],MPI::CHAR,0);

	

	resources.database_info=databaseinfo;
	
	
#if 0

	cout<<"number database chunk:"<<databaseinfo.number_chunks<<endl;
	//cout<<"number sequences:"<<number_sequences<<endl;
	

#endif
	for(int i=0;i<databaseinfo.number_chunks;i++){
		DatabaseResource database;
		database.chunk_id=i;
		database.available=false;
		resources.database_list.push_back(database);
	}
	

	
	
	//init task pool
	queue<AlignmentTask> task_queue;
	
	resources.task_pool.push_back(task_queue);
	

	

	//init worker process
	
	int size[3];
	size[0]=resources.query_list.size();
	size[1]=resources.database_list.size();
	size[2]=sizeof(database_filename.c_str());
	MPI::COMM_WORLD.Bcast(size,3,MPI::INT,0);
	MPI::COMM_WORLD.Bcast((char *)database_filename.c_str(),size[2],MPI::CHAR,0);
	//deliverer.CreateResponseThread(resources,mpi_parameter.size);
	//int dbinfo_size[5];
	//MPI::Bcast(dbinfo_size,5,MPI::INT,0);
	
	
}
void MPICommon::SetupWorkerResources(WorkerResources &resources,MPIParameter &mpi_parameter){
	int query_chunk_size;
	int database_chunk_size;
	char* database_filename;
	int buf[3];
	DatabaseInfo databaseinfo;
	//init communicator;
	
	/*
	for(int i=1;i<mpi_parameter.size;i++){
		if(i==mpi_parameter.rank){
			resources.master_comm = MPI::COMM_WORLD.Split(0,i);
			cout<<"worker comm:"<<resources.master_comm.Get_size()<<endl;
			cout<<resources.master_comm<<endl;
    		}
		MPI::COMM_WORLD.Split(1,i);
		
	}
	*/
	
	MPI::COMM_WORLD.Bcast((char *)&databaseinfo.number_chunks,sizeof(databaseinfo.number_chunks),MPI::CHAR,0);
	MPI::COMM_WORLD.Bcast((char *)&databaseinfo.max_sequence_length,sizeof(databaseinfo.max_sequence_length),MPI::CHAR,0);
	MPI::COMM_WORLD.Bcast((char *)&databaseinfo.database_length,sizeof(databaseinfo.database_length),MPI::CHAR,0);
	MPI::COMM_WORLD.Bcast((char *)&databaseinfo.number_sequences,sizeof(databaseinfo.number_sequences),MPI::CHAR,0);
	MPI::COMM_WORLD.Bcast((char *)&databaseinfo.sequence_delimiter,sizeof(databaseinfo.sequence_delimiter),MPI::CHAR,0);

	
	MPI::COMM_WORLD.Bcast(buf,3,MPI::INT,0);
	query_chunk_size=buf[0];
	database_chunk_size=buf[1];
	database_filename= new char[buf[2]];
	MPI::COMM_WORLD.Bcast(database_filename,buf[2],MPI::CHAR,0);
	resources.database_filename=string(database_filename);
	resources.database_info=databaseinfo;
	//cout<<"query,database:"<<query_chunk_size<<","<<database_chunk_size<<endl;
	//cout<<"databasename:"<<resources.database_filename<<endl;
	//setup query
	for(int i=0;i<query_chunk_size;i++){
		QueryResource query;
		query.chunk_id=i;
		query.available=false;
		resources.query_list.push_back(query);
	}
	//setup database
	
	for(int i=0;i<database_chunk_size;i++){
		DatabaseResource database;
		database.chunk_id=i;
		database.available=false;
		resources.database_list.push_back(database);
	}
	
	
	
	
}
void MPICommon::BuildQueryChunkPointers(string &queries_filename,vector<int> &chunk_pointer_list,
										vector<int> &chunk_size_list,AligningParameters &parameter){
	
	ifstream isf(queries_filename.c_str());
	Queries::Parameters queries_parameters;
	AlignerCommon::BuildQueriesParameters(parameter,queries_parameters);
	Queries queries(queries_parameters);

	int i=0;
	bool eof_flag=false;
	while(!eof_flag){
		ifstream::pos_type begin = isf.tellg();
		ifstream::pos_type isf_ptr=queries.GetNextChunkPosition(isf);
		if(isf_ptr==-1){
			isf.clear();
			isf.seekg(begin);
			isf_ptr=isf.seekg(0,std::ios_base::end).tellg();
			eof_flag=true;
		   
		}
	
		int size=isf_ptr-begin;
		char *array = new char[size];
		isf.clear();

		isf.seekg(begin);
		isf.read(array,size);
		isf.seekg(isf_ptr);
		i++;
		chunk_pointer_list.push_back(begin);
		chunk_size_list.push_back(size);
		delete [] array;
	}
}
void MPICommon::LoadQueryResource(MasterResources &resources,int chunk_id){
   
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
void MPICommon::UnloadQueryResource(MasterResources &resources,int chunk_id){
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

void MPICommon::LoadDatabaseResource(WorkerResources &resources,int chunk_id){
	if(chunk_id > resources.database_list.size()){
		cerr<<"database chunk id error:"<<chunk_id<<endl;
		MPI_Abort(MPI_COMM_WORLD,1);
	}
	if(resources.database_list[chunk_id].available){
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
void MPICommon::UnloadDatabaseResource(WorkerResources &resources,int chunk_id){
	if(chunk_id > resources.database_list.size()){
		cerr<<"database chunk id error:"<<chunk_id<<endl;
		MPI_Abort(MPI_COMM_WORLD,1);
	}
	if(!resources.database_list[chunk_id].available){
		return ;
	}
	delete [] resources.database_list[chunk_id].inf;
	delete [] resources.database_list[chunk_id].nam;
	delete [] resources.database_list[chunk_id].off;
	delete [] resources.database_list[chunk_id].seq;
	delete [] resources.database_list[chunk_id].scp;
	delete [] resources.database_list[chunk_id].sdp;
	
	resources.database_list[chunk_id].available=false;

}

void MPICommon::loadFileData(std::string filename,char **ptr,uint64_t *size){
	ifstream in(filename.c_str());
	uint64_t begin,end;
	int size_;
	in.seekg(0,ios::end);
	end = in.tellg();
	in.clear();
	in.seekg(0,ios::beg);
	begin=in.tellg();
	size_=end-begin;
	
	*ptr = new char[size_];
	//cout<<filename<<" size:"<<size_<<endl;	
	in.read(*ptr,size_);
	in.close();
	*size=size_;

}


void MPICommon::AcceptCommand(MasterResources &resources){
	int cmd[2];
	
	MPI::Status status;
	MPI::COMM_WORLD.Recv(&cmd,2,MPI::INT,MPI::ANY_SOURCE,MPI::ANY_TAG,status);
	
	switch(cmd[0]){
	case MPIResource::CMD_RequestQuery :
		// cmd[0] = RequestQuery cmd[1] = target_chunk
		AcceptRequestQuery(cmd,resources,status);
	  
		break;
	case MPIResource::CMD_RequestDatabase :
		AcceptRequestDatabase(cmd,resources,status);
		break;
	case MPIResource::CMD_RequestTask :
		AcceptRequestTask(cmd,resources,status);
		break;
	default :
		break;
	}
	
}

void MPICommon::AcceptRequestQuery(int cmd[2],MasterResources &resources,MPI::Status &status){
	int target_rank=status.Get_source();
	int target_chunk=cmd[1];

	if(resources.query_list[target_chunk].available==false){// target query not loaded
		cmd[0]=MPIResource::NACK;
		cout<<"NACK"<<endl;
		MPI::COMM_WORLD.Send(cmd,2,MPI::INT,target_rank,0);
		return ;
	}
	cmd[0]=MPIResource::ACK;
	cmd[1]=resources.query_list[target_chunk].size;

	MPI::COMM_WORLD.Send(cmd,2,MPI::INT,target_rank,0);
	MPI::COMM_WORLD.Send(resources.query_list[target_chunk].data,cmd[1],MPI::CHAR,target_rank,0);
	
}
void MPICommon::AcceptRequestDatabase(int cmd[2],MasterResources &resources,MPI::Status &status){
	int target_rank=status.Get_source();
	int target_chunk=cmd[1];

	if(resources.database_list[target_chunk].available==false){
		cmd[0]=MPIResource::NACK;
		MPI::COMM_WORLD.Send(cmd,2,MPI::INT,target_rank,0);
		return ;
	}
   
	cmd[0]=MPIResource::ACK;

		cout<<"acceptDatabase Send."<<endl;
   
	MPI::COMM_WORLD.Send(cmd,2,MPI::INT,target_rank,0);

	
	MPI::COMM_WORLD.Send((char *) &resources.database_list[target_chunk].inf_size,
						  sizeof(&resources.database_list[target_chunk].inf_size),
						  MPI::CHAR,target_rank,0);
	MPI::COMM_WORLD.Send((char *) &resources.database_list[target_chunk].nam_size,
						  sizeof(&resources.database_list[target_chunk].nam_size),
						  MPI::CHAR,target_rank,0);
	MPI::COMM_WORLD.Send((char *) &resources.database_list[target_chunk].off_size,
						  sizeof(&resources.database_list[target_chunk].off_size),
						  MPI::CHAR,target_rank,0);
	MPI::COMM_WORLD.Send((char *) &resources.database_list[target_chunk].seq_size,
						  sizeof(&resources.database_list[target_chunk].seq_size),
						  MPI::CHAR,target_rank,0);
	MPI::COMM_WORLD.Send((char *) &resources.database_list[target_chunk].scp_size,
						  sizeof(&resources.database_list[target_chunk].scp_size),
						  MPI::CHAR,target_rank,0);
	MPI::COMM_WORLD.Send((char *) &resources.database_list[target_chunk].sdp_size,
						  sizeof(&resources.database_list[target_chunk].sdp_size),
						  MPI::CHAR,target_rank,0);
	//sources.comm_list[target_rank].Barrier();

	MPI::COMM_WORLD.Send(resources.database_list[target_chunk].inf,
						  resources.database_list[target_chunk].inf_size,
						  MPI::CHAR,target_rank,0);
	MPI::COMM_WORLD.Send(resources.database_list[target_chunk].nam,
						  resources.database_list[target_chunk].nam_size,
						  MPI::CHAR,target_rank,0);
	MPI::COMM_WORLD.Send(resources.database_list[target_chunk].off,
						  resources.database_list[target_chunk].off_size,
						  MPI::CHAR,target_rank,0);
	MPI::COMM_WORLD.Send(resources.database_list[target_chunk].seq,
						  resources.database_list[target_chunk].seq_size,
						  MPI::CHAR,target_rank,0);
	MPI::COMM_WORLD.Send(resources.database_list[target_chunk].scp,
						  resources.database_list[target_chunk].scp_size,
						  MPI::CHAR,target_rank,0);
	MPI::COMM_WORLD.Send(resources.database_list[target_chunk].sdp,
						  resources.database_list[target_chunk].sdp_size,
						  MPI::CHAR,target_rank,0);

						  
	
}
void MPICommon::AcceptRequestTask(int cmd[2],MasterResources &resource,MPI::Status &status){
	int target_rank=status.Get_source();

}

void MPICommon::UpdateTaskBalance(MasterResources &resources,MPIParameter &mpi_resource){
	int mpi_size = mpi_resource.size;
	int database_chunk = resources.database_list.size();
	int *chunk_count = new int[database_chunk];

	for(int i=0;i<database_chunk;i++){
		chunk_count[i]=0;		
	}
	for(int i=0;i<resources.node_loading_database.size();i++){
		
	}
    
	


}

void MPICommon::GetNextTask(MasterResources &resources,int target,AlignmentTask &task){

	
	
}







bool MPICommon::BuildParameters(int argc, char* argv[], string &input_filename,
								 string &database_filename, string &output_filename,
								 AligningParameters &parameters) {
	 int c;
	 extern char *optarg;
	 extern int optind;
	 optind = 1;
	 string score_matrix_filename;
	 Logger *logger = Logger::GetInstance();
	 const string default_protein_matrix_name = "BLOSUM62";
	 const string default_protein_matrix =
		 "   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *\nA  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4\n R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4\n N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4\n D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4\n C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4\n Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4\n E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4\n G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4\n H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4\n I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4\n L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4\n K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4\n M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4\n F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4\n P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4\n S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4\n T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4\n W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4\n Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4\n V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4\n B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4\n Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4\n X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4\n * -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1\n";

	 parameters.number_threads = 1;
	 parameters.filter = true;
	 parameters.queries_file_sequence_type_ptr = std::tr1::shared_ptr<
	 SequenceType>(new ProteinType());
	 parameters.aligning_sequence_type_ptr = std::tr1::shared_ptr<SequenceType>(
																				new ProteinType());
	 parameters.queries_chunk_size = 1 << 27;
	 ScoreMatrixReader score_matrix_reader;
	 vector<int> matrix;
	 unsigned int number_letters;
	 istringstream default_protein_score_matrix_is(default_protein_matrix);
	 score_matrix_reader.Read(default_protein_score_matrix_is,
							  *(parameters.aligning_sequence_type_ptr), matrix, number_letters);
	 parameters.score_matrix = ScoreMatrix(default_protein_matrix_name,
										   &matrix[0], number_letters);
	 parameters.gap_open = 11;
	 parameters.gap_extension = 1;
	 parameters.number_gpus = -1;
	 parameters.normalized_presearched_ungapped_extension_cutoff = 7.0;
	 parameters.normalized_presearched_gapped_extension_trigger = 22.0;
	 parameters.normalized_presearched_gapped_extension_cutoff = 15.0;
	 parameters.normalized_result_gapped_extension_cutoff = 25.0;
	 parameters.max_number_results = 10;
	 parameters.max_number_one_subject_results = 1;
     while ((c = getopt(argc, argv, "a:b:d:F:g:h:l:i:o:q:t:v:")) != -1) {
		 switch (c) {
		 case 'a':
			 parameters.number_threads = atoi(optarg);
			 break;
		 case 'b':
			 parameters.max_number_results = atoi(optarg);
			 break;
		 case 'd':
			 database_filename = optarg;
			 break;
		 case 'F':
			 if (optarg[0] == 'T') {
				 parameters.filter = true;
			 } else if (optarg[0] == 'F') {
				 parameters.filter = false;
			 } else {
				 logger->ErrorLog("invalid option, -F is T or F.");
				 return false;
			 }
			 break;
		 case 'g':
			 parameters.number_gpus = atoi(optarg);
			 break;
		 case 'l':
			 parameters.queries_chunk_size = atoi(optarg);
			 break;
		 case 'i':
			 input_filename = optarg;
			 break;
		 case 'o':
			 output_filename = optarg;
			 break;
		 case 'q':
			 if (optarg[0] == 'p') {
				 parameters.queries_file_sequence_type_ptr =
					 std::tr1::shared_ptr<SequenceType>(new ProteinType());
			 } else if (optarg[0] == 'd') {
				 parameters.queries_file_sequence_type_ptr =
					 std::tr1::shared_ptr<SequenceType>(new DnaType());
			 } else {
				 logger->ErrorLog("invalid option, -q is p or d.");
				 return false;
			 }
			 break;
		 case 't':
			 if (optarg[0] == 'p') {
				 parameters.aligning_sequence_type_ptr = std::tr1::shared_ptr<
				 SequenceType>(new ProteinType());
			 } else if (optarg[0] == 'd') {
				 parameters.aligning_sequence_type_ptr = std::tr1::shared_ptr<
				 SequenceType>(new DnaType());
			 } else {
				 logger->ErrorLog("invalid option, -q is p or d.");
				 return false;
			 }
			 break;
		 case 'v':
			 parameters.max_number_one_subject_results = atoi(optarg);
			 break;
		 default:
			 logger->ErrorLog("\nTry `ghostz --help' for more information.");
			 return false;
		 }
	 }
	 return true;
}
