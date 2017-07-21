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
#include "aligner_mpi.h"
#include "aligner.h"
#include "result_summarizer.h"

#include "mpi.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <sys/stat.h>

#define F_TIMER
#ifdef F_TIMER
#include <sys/time.h>

#endif


using namespace std;



void MPICommon::Run(string &queries_filename,string &database_filename,
					string &output_filename,string &tmp_dirname,
					AligningParameters &parameter,MPIParameter &mpi_parameter){
	
	if(mpi_parameter.rank==0){
		RunMaster(queries_filename,database_filename,output_filename,tmp_dirname,
				  parameter,mpi_parameter);
	}else{
		RunWorker(parameter,mpi_parameter,output_filename,tmp_dirname);
		}
}

void MPICommon::RunGPU(string &queries_filename,string &database_filename,
					   string &output_filename,string &tmp_dirname,
					   AligningParameters &parameter,MPIParameter &mpi_parameter){
	
	if(mpi_parameter.rank==0){
		RunMaster(queries_filename,database_filename,output_filename,tmp_dirname,
				  parameter,mpi_parameter);
	}else{
		RunWorkerGPU(parameter,mpi_parameter,output_filename,tmp_dirname);
	}
}



void MPICommon::RunMaster(string &queries_filename,string &database_filename,
						  string &output_filename,string &tmp_dirname,
						  AligningParameters &parameter,MPIParameter &mpi_parameter){
#ifdef F_TIMER
	struct timeval init_tv;
	struct timeval tv;
	gettimeofday(&init_tv,NULL);
	float timer_second;
#endif
	//Master Process Init
	cout<<"master start:"<<mpi_parameter.rank<<endl;
	MasterResources resources;
	resources.query_filename=queries_filename;
	resources.database_filename=database_filename;
	resources.output_filename=output_filename;
	if(tmp_dirname.size()==0){
		cout<<"tmp directory error."<<endl;
		MPI::COMM_WORLD.Abort(1);
	}
	struct stat st;
	mkdir(tmp_dirname.c_str(),0755);

	if(stat(tmp_dirname.c_str(),&st)==-1){
		cout<<"tmp directory error"<<endl;
		MPI::COMM_WORLD.Abort(1);
	}
	
	//***********************//
	//Start Init Phase

	
	SetupMasterResources(queries_filename,database_filename,resources,parameter,mpi_parameter);
#ifdef F_TIMER
	gettimeofday(&tv,NULL);
	timer_second = (tv.tv_sec-init_tv.tv_sec)*1000 + (tv.tv_usec - init_tv.tv_usec)*0.001;	
	printf("%.0f\t[ms]\tmaster: SetupResource\n",timer_second);
#endif
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
		MPIResource::LoadQueryResource(resources,i);
		
	}   //loading all query on master memory
	
	//	cout<<resources.query_list[0].data;
	
	
#ifdef F_TIMER
	gettimeofday(&tv,NULL);
	timer_second = (tv.tv_sec-init_tv.tv_sec)*1000 + (tv.tv_usec - init_tv.tv_usec)*0.001;	
	printf("%.0f\t[ms]\tmaster: Init Phase\n",timer_second);
#endif
   //End Init Phase
	//***********************//
	//Start Search Phase 
	int task_remain=0;
	LoadBalancer balancer(resources.query_list.size(),resources.database_list.size(),mpi_parameter.size);
	for(int i=1;i<mpi_parameter.size;i++){
		
		balancer.SetDatabaseLoadingMap(i,(i-1)%resources.database_list.size());
 		//cout<<"rank:"<<i<<"  db:"<<(i-1)%resources.database_list.size()<<endl;								   
	}
	/*
	int num_q=6,num_d=3,num_p=4;
	LoadBalancer balancer(num_q,num_d,num_p);
	for(int i=1;i<num_p;i++){
		balancer.SetDatabaseLoadingMap(i,(i-1)%(num_d));
		
		
		}
	*/
	
	cout<<"total task:"<<balancer.GetRemainTask()<<endl;
	cout<<balancer.print()<<endl;
	count_terminate=0;

	for(;;){
		AcceptCommand(resources,balancer);
		if(count_terminate == mpi_parameter.size-1){
			break;
		}
		//		cout<<balancer.print()<<endl;
		
	}	

#ifdef F_TIMER
	gettimeofday(&tv,NULL);
	timer_second = (tv.tv_sec-init_tv.tv_sec)*1000 + (tv.tv_usec - init_tv.tv_usec)*0.001;	
	printf("%.0f\t[ms]\tmaster: Search Phase\n",timer_second);
#endif
	MPI::COMM_WORLD.Barrier();
	//End Search Phase
	//***********************//
	//Start Report Phase
	string tmp("tmp");
	ResultSummarizer summary(tmp_dirname,output_filename);
	vector<vector<Result> > results_list;
 	AlignmentTask task;
	summary.GatherResultMaster(resources.query_list.size(),resources.database_list.size(),
							   parameter,resources.database_info);
	
	
	//finalize
	MPI::COMM_WORLD.Barrier();	
	
}
void MPICommon::RunWorker(AligningParameters &parameter,MPIParameter &mpi_parameter,
						  string &output_filename,string &tmp_dirname){
	
	//Worker Process Init
#ifdef F_TIMER
	struct timeval init_tv;
	struct timeval tv;
	gettimeofday(&init_tv,NULL);
	float timer_second;
#endif

	cout<<"worker start:"<<mpi_parameter.rank<<endl;
	WorkerResources resources;
	int rank=mpi_parameter.rank;
	struct stat st;
	mkdir(tmp_dirname.c_str(),0755);

	if(stat(tmp_dirname.c_str(),&st)==-1){
		cout<<"tmp directory error"<<endl;
		MPI::COMM_WORLD.Abort(1);
	}
	//***********************//
	//Start Init Phase
	
	SetupWorkerResources(resources,mpi_parameter);

	//Create subgroup comm
	
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
	//	cout<<"rank:"<<rank<<"comm created."<<endl;
	if(submaster){
		//load db from filesystem
		MPIResource::LoadDatabaseResource(resources,target_chunk);

	}
	//Broadcast db to subgroup
	uint64_t sum =0;
	MPIResource::BcastDatabase(resources.database_list[target_chunk],resources.subgroup_comm,0);
	//cout<<"rank:"<<rank<<"db recved."<<endl;

	DatabaseType database(resources.database_list[target_chunk],resources.database_info);
	//DatabaseType database_(resources.database_filename);
	//cout<<"database.GetChunkId = "<<database.GetChunkId()<<endl;
#ifdef F_TIMER
	gettimeofday(&tv,NULL);
	timer_second = (tv.tv_sec-init_tv.tv_sec)*1000 + (tv.tv_usec - init_tv.tv_usec)*0.001;
	printf("%.0f\t[ms]\trank:%d Init Phase\n",timer_second,rank);
#endif
	//End Init Phase
	//***********************//
	//Start Search Phase
	AlignerMPI aligner;
	AlignmentTask task;
	string tmpDir="tmp";
	ResultSummarizer summary(tmp_dirname,output_filename);
	
	vector<vector<Result> > results_list;
	for(;;){
		int task_query=0,task_db=0;
		MPIResource::RequestTask(target_chunk,task);
		task_query=task.query_chunk;
		task_db=task.database_chunk;
	
		//cout<<"rank:"<<rank<<" recv: "<<task_query<<","<<task_db<<endl;

		if(task.query_chunk==-1){
			cout<<"break rank:"<<rank<<endl;
			break;
		}
		if(task.database_chunk!=target_chunk){
			//Loadatabase
			cout<<"rank:"<<rank<<" switching db "<<target_chunk<<" to "<<task.database_chunk<<endl;
			//cout<<"rank:"<<rank<<" unload db:"<<target_chunk<<endl;
			MPIResource::UnloadDatabaseResource(resources,target_chunk);
			target_chunk=task.database_chunk;
			//cout<<"rank:"<<rank<<" load db:"<<target_chunk<<endl;
			MPIResource::LoadDatabaseResource(resources,target_chunk);
			database.SetChunk(resources.database_list[target_chunk]);
			//cout<<"rank:"<<rank<<" end Loading."<<endl;
		}
		
		if(!resources.query_list[task.query_chunk].available){
			MPIResource::RequestQuery(resources,task.query_chunk);
		}
		aligner.Search(resources.query_list[task.query_chunk],
					   database,results_list,parameter,mpi_parameter);
		
		MPIResource::UnloadQueryResource(resources,task.query_chunk);
		summary.SaveResultFile(results_list,task);
			
		results_list.clear();
	
#ifdef F_TIMER
		gettimeofday(&tv,NULL);
		timer_second = (tv.tv_sec-init_tv.tv_sec)*1000 + (tv.tv_usec - init_tv.tv_usec)*0.001;
		printf("%.0f\t[ms]\trank:%d Search(%d,%d)\n",timer_second,rank,task_query,task_db);
#endif
	}
	
	
	
	Aligner aligner_;
	string q_name=string("/work1/t2ggenome/goto/work/ghostz/mpi/test/query/ERR315856.fasta_randN1k");
	string out_name = string("out");
	//
	//	aligner_.Align(q_name,  resources.database_filename,out_name,parameter);
	for(int i=0;i<results_list.size();i++){
		for(int j=0;j<results_list[i].size();j++){
			//cout<<results_list[i][j].subject_name<<endl;
		}}
	//unload resources;
	MPIResource::UnloadDatabaseResource(resources,target_chunk);
   /*
	for(int i=0;i<resources.query_list.size();i++){
		if(resources.query_list[i].size>0){
			cout<<"rank:"<<rank<<" query:"<<i<<"  True "<<resources.query_list[i].size<<endl;
		}else{
			cout<<"rank:"<<rank<<" query:"<<i<<"  False"<<endl;
			
		}
	}
	//*/
#ifdef F_TIMER
	gettimeofday(&tv,NULL);
	timer_second = (tv.tv_sec-init_tv.tv_sec)*1000 + (tv.tv_usec - init_tv.tv_usec)*0.001;
	printf("%.0f\t[ms]\trank:%d Search Phase\n",timer_second,rank);
#endif
	MPI::COMM_WORLD.Barrier();
	//End Search Phase
	//***********************//
	//Start Report Phase
	summary.GatherResultWorker(rank,parameter,resources.database_info);
	//End Report Phase

	
	//finalize

	MPI::COMM_WORLD.Barrier();
}


void MPICommon::RunWorkerGPU(AligningParameters &parameter,MPIParameter &mpi_parameter,
							 string &output_filename, string &tmp_dirname){

}

void MPICommon::SetupMasterResources(string &queries_filename,string &database_filename,
									 MasterResources &resources,AligningParameters &parameter,
									 MPIParameter &mpi_parameter){
	//init query resource
	vector<int> pointer_list;
	vector<int> size_list;
	struct stat st;
	if(stat(queries_filename.c_str(),&st)==-1){
		cerr<<"Query file does not exist. Abort."<<endl;
		MPI::COMM_WORLD.Abort(1);
	}
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
	if(stat(inf_file.c_str(),&st)==-1){
		cerr<<"Database file does not exist. Abort."<<endl;
		MPI::COMM_WORLD.Abort(1);
	
	}
	
	MPIResource::LoadDatabaseInfo(databaseinfo,inf_file);
	
	MPIResource::BcastDatabaseInfo(databaseinfo,MPI::COMM_WORLD,0);
	

	resources.database_info=databaseinfo;
	
	
#if 1

	cout<<"database chunk:"<<databaseinfo.number_chunks<<endl;
	//cout<<"number sequences:"<<number_sequences<<endl;
	

#endif
	for(int i=0;i<databaseinfo.number_chunks;i++){
		DatabaseResource database;
		database.chunk_id=i;
		database.available=false;
		resources.database_list.push_back(database);
	}
	
	

	//init worker process
	const char *filename = database_filename.c_str();
	cout<<"master:"<<filename<<":size="<<sizeof(filename)<<endl;
	int size[3];
	size[0]=resources.query_list.size();
	size[1]=resources.database_list.size();
	size[2]=database_filename.length();
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
	
	
	MPIResource::BcastDatabaseInfo(databaseinfo,MPI::COMM_WORLD,0);
		
	MPI::COMM_WORLD.Bcast(buf,3,MPI::INT,0);
	query_chunk_size=buf[0];
	database_chunk_size=buf[1];
	database_filename= new char[buf[2]];
	MPI::COMM_WORLD.Bcast(database_filename,buf[2],MPI::CHAR,0);
	resources.database_filename=string(database_filename,database_filename+buf[2]);
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
	bool flag=true;
	std::ifstream::pos_type end_ptr;

	while(flag){
		ifstream::pos_type begin = isf.tellg();
		ifstream::pos_type isf_ptr;
		flag=queries.GetNextChunkPosition(isf,&isf_ptr);
		
		/*if(isf_ptr==-1){
			isf.clear();
			isf.seekg(begin);
			isf_ptr=isf.seekg(0,std::ios_base::end).tellg();
			eof_flag=true;
		   
			}*/
	
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


void MPICommon::AcceptCommand(MasterResources &resources,LoadBalancer &balancer){
	int cmd;//(command,target_chunk;
	int target;
	int dst_rank;
	dst_rank=MPIResource::AcceptCommand(resources,&cmd,&target);
	//cout<<"cmd "<<cmd<<" "<<target<<" "<<MPIResource::CMD_RequestTask<<endl;
    switch(cmd){
    case MPIResource::CMD_RequestQuery :

        MPIResource::AcceptRequestQuery(resources,target,dst_rank);

        break;
    case MPIResource::CMD_RequestDatabase :
        //AcceptRequestDatabase(cmd,resources,status);
        break;
    case MPIResource::CMD_RequestTask :

        AlignmentTask task;
		if(balancer.GetRemainTask()==0){
			task.query_chunk=-1;
			task.database_chunk=-1;
			count_terminate++;
		}else{
			task=balancer.GetNext(target);
			
			if(task.query_chunk==-1){//no more task in target_chunk;
				int switch_target = balancer.GetSwitchTargetChunk();
				balancer.SetDatabaseLoadingMap(dst_rank,switch_target);
				task=balancer.GetNext(switch_target);
			}
		}
		//cout<<task.query_chunk<<","<<task.database_chunk<<endl;
		MPIResource::AcceptRequestTask(task,dst_rank);
        break;
    default :
        break;
	}

	
}





