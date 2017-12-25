/*
 * mpi_common.cpp
 *  Created on 2017/05/16
 *  Author:goto
 */

#include "mpi_common.h"
#include "aligner_mpi.h"
#include "aligner_gpu_mpi.h"
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
#include <limits>

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
		bool useGPU=false;
		RunWorker(queries_filename,database_filename,output_filename,tmp_dirname,
				  parameter,mpi_parameter,useGPU);
	}
}

void MPICommon::RunGPU(string &queries_filename,string &database_filename,
					   string &output_filename,string &tmp_dirname,
					   AligningParameters &parameter,MPIParameter &mpi_parameter){
	
	if(mpi_parameter.rank==0){
		RunMaster(queries_filename,database_filename,output_filename,tmp_dirname,
				  parameter,mpi_parameter);
	}else{
		bool useGPU=true;
	    RunWorker(queries_filename,database_filename,output_filename,tmp_dirname,
				  parameter,mpi_parameter,useGPU);
	}
}



void MPICommon::RunMaster(string &queries_filename,string &database_filename,
						  string &output_filename,string &tmp_dirname,
						  AligningParameters &parameter,MPIParameter &mpi_parameter){
#ifdef F_TIMER
	struct timeval init_tv;
	struct timeval tv;
	gettimeofday(&init_tv,NULL);
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
	mkdir(tmp_dirname.c_str(),0777);

	if(stat(tmp_dirname.c_str(),&st)==-1){
		cout<<"tmp directory error"<<endl;
		MPI::COMM_WORLD.Abort(1);
	}
	
	//***********************//
	//Start Init Phase

	
	SetupDatabaseResourcesMaster(database_filename,resources,parameter,mpi_parameter);
#ifdef F_TIMER 
	gettimeofday(&tv,NULL);
	printf("%.0lf\t[ms]\tmaster: SetupResource\n",timevalToMillisec(tv)-timevalToMillisec(init_tv));
#endif
	
	
	SetupQueryResourcesMaster(queries_filename,resources,parameter,mpi_parameter);
#if 0
	for(int i=0;i<resources.query_list.size();i++){
		cout<<resources.query_list[i].chunk_id;
		cout<<" ptr:"<<resources.query_list[i].ptr;
		cout<<" size:"<<resources.query_list[i].size<<endl;
	}
#endif
	
	
#ifdef F_TIMER
	gettimeofday(&tv,NULL);
	printf("%.0lf\t[ms]\tmaster: Init Phase\n",timevalToMillisec(tv)-timevalToMillisec(init_tv));
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
	printf("%.0lf\t[ms]\tmaster: Search Phase\n",timevalToMillisec(tv)-timevalToMillisec(init_tv));
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
	
	remove(tmp_dirname.c_str());
	
#ifdef F_TIMER
	gettimeofday(&tv,NULL);
	printf("%.0lf\t[ms]\tmaster: Report Phase\n",timevalToMillisec(tv)-timevalToMillisec(init_tv));
#endif
	//finalize
	MPI::COMM_WORLD.Barrier();	
	
}
void MPICommon::RunWorker(string &queries_filename,string &database_filename,
						  string &output_filename,string &tmp_dirname,
						  AligningParameters &parameter,MPIParameter &mpi_parameter,
						  bool useGPU){
	
	//Worker Process Init

#ifdef F_TIMER
	struct timeval init_tv;
	struct timeval tv;
	gettimeofday(&init_tv,NULL);
#endif
 
	cout<<"worker start:"<<mpi_parameter.rank<<endl;
	WorkerResources resources;
	resources.query_filename=queries_filename;

	int rank=mpi_parameter.rank;
	struct stat st;
	mkdir(tmp_dirname.c_str(),0755);

	if(stat(tmp_dirname.c_str(),&st)==-1){
		cout<<"tmp directory error"<<endl;
		MPI::COMM_WORLD.Abort(1);
	}
	//***********************//
	//Start Init Phase
	
	SetupDatabaseResourcesWorker(database_filename,resources,parameter,mpi_parameter);

	//Create subgroup comm
	
	int db_chunk_size = resources.database_list.size();
	int target_chunk=(rank-1)%db_chunk_size;
	
	
	DatabaseType database(resources.database_info);
	DatabaseTypeGpu database_gpu(resources.database_info);
	
	if(useGPU){
		database_gpu.SetChunk(resources.database_list[target_chunk]);
	}else{
		database.SetChunk(resources.database_list[target_chunk]);
	}
	
	SetupQueryResourcesWorker(resources,parameter,mpi_parameter);

#ifdef F_TIMER
	gettimeofday(&tv,NULL);
	printf("%.0lf\t[ms]\trank:%d Init Phase\n",timevalToMillisec(tv)-timevalToMillisec(init_tv),rank);
#endif
	//End Init Phase
	//***********************//
	//Start Search Phase

	AlignerMPI aligner;
	AlignerGpuMPI aligner_g;
	
	AlignmentTask task;
	string tmpDir="tmp";
	ResultSummarizer summary(tmp_dirname,output_filename);
	uint64_t timer_align=0 ,timer_align_start,timer_align_end;
	uint64_t timer_task=0, timer_task_start, timer_task_end;
	vector<vector<Result> > results_list;
	bool query_direct_load=true;
	for(;;){
		int task_query=0,task_db=0;
#ifdef F_TIMER
		gettimeofday(&tv,NULL);
		timer_task_start = (tv.tv_sec-init_tv.tv_sec)*1000 + (tv.tv_usec - init_tv.tv_usec)*0.001;
		
#endif
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
			if(useGPU){
				database_gpu.SetChunk(resources.database_list[target_chunk]);
			}else{
				database.SetChunk(resources.database_list[target_chunk]);
			}//cout<<"rank:"<<rank<<" end Loading."<<endl;
		}
		
		if(!resources.query_list[task.query_chunk].available){
			if(query_direct_load){
				MPIResource::LoadQueryResource(resources,task.query_chunk);
			}else{
				MPIResource::RequestQuery(resources,task.query_chunk);
			}
		}

#ifdef F_TIMER
	gettimeofday(&tv,NULL);
   
	timer_align_start = (tv.tv_sec-init_tv.tv_sec)*1000 + (tv.tv_usec - init_tv.tv_usec)*0.001;
	timer_task_end=timer_align_start;
	timer_task+=timer_task_end - timer_task_start;
	
	printf("%.0lf\t[ms]\trank:%d Start Task\n",timevalToMillisec(tv)-timevalToMillisec(init_tv),rank);
#endif

		if(useGPU){
			
			aligner_g.Search(resources.query_list[task.query_chunk],
						   database_gpu,results_list,parameter,mpi_parameter);
		}else{
			aligner.Search(resources.query_list[task.query_chunk],
					   database,results_list,parameter,mpi_parameter);
		}
#ifdef F_TIMER
	gettimeofday(&tv,NULL);
	timer_align_end = (tv.tv_sec-init_tv.tv_sec)*1000 + (tv.tv_usec - init_tv.tv_usec)*0.001;
	timer_align+= timer_align_end-timer_align_start;
#endif

		MPIResource::UnloadQueryResource(resources,task.query_chunk);
		summary.SaveResultFile(results_list,task);
			
		results_list.clear();
	
#ifdef F_TIMER
		gettimeofday(&tv,NULL);
		printf("%.0lf\t[ms]\trank:%d Search(%d,%d)\n",timevalToMillisec(tv)-timevalToMillisec(init_tv),rank,task_query,task_db);
#endif
	}
	
#ifdef F_TIMER
	cout<<"rank:"<<rank<<"\tCumulative Effective Search Time : "<<timer_align<< " [ms]"<<endl;
	cout<<"rank:"<<rank<<"\tCumulative Task Recv Time : "<<timer_task<<" [ms]"<<endl;
#endif
	

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
	printf("%.0lf\t[ms]\trank:%d Search Phase\n",timevalToMillisec(tv)-timevalToMillisec(init_tv),rank);
#endif
	MPI::COMM_WORLD.Barrier();
	//End Search Phase
	//***********************//
	//Start Report Phase
	summary.GatherResultWorker(rank,parameter,resources.database_info);
	//End Report Phase

	//remove(tmp_dirname.c_str());
#ifdef F_TIMER
	gettimeofday(&tv,NULL);
	printf("%.0lf\t[ms]\trank:%d Report Phase\n",timevalToMillisec(tv)-timevalToMillisec(init_tv),rank);
#endif
	//finalize

	MPI::COMM_WORLD.Barrier();
}




void MPICommon::SetupDatabaseResourcesMaster(string &database_filename,
											 MasterResources &resources,AligningParameters &parameter,
											 MPIParameter &mpi_parameter){
	struct stat st;
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
	

	for(int i=0;i<databaseinfo.number_chunks;i++){
		DatabaseResource database;
		database.chunk_id=i;
		database.available=false;
		resources.database_list.push_back(database);
	}
	
	//Create Subgroup comm
	for(int i=0;i<resources.database_list.size();i++){
		MPI::COMM_WORLD.Split(1,0);
	}

}
void MPICommon::SetupDatabaseResourcesWorker(string &database_filename,WorkerResources &resources,
											 AligningParameters &parameter,MPIParameter &mpi_parameter){
	
	DatabaseInfo databaseinfo;
	MPIResource::BcastDatabaseInfo(databaseinfo,MPI::COMM_WORLD,0);
	resources.database_info=databaseinfo;
	resources.database_filename=database_filename;;
	
	//setup database
	for(int i=0;i<databaseinfo.number_chunks;i++){
		DatabaseResource database;
		database.chunk_id=i;
		database.available=false;
		resources.database_list.push_back(database);
	}
	
	int db_chunk_size = resources.database_list.size();
	bool submaster = false;
	int rank = mpi_parameter.rank;
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
	

	
}
void MPICommon::SetupQueryResourcesMaster(string &queries_filename,MasterResources &resources,
										 AligningParameters &parameter,MPIParameter &mpi_parameter){
	vector<uint64_t> pointer_list;
	vector<uint64_t> size_list;
	struct stat st;
	if(stat(queries_filename.c_str(),&st)==-1){
		cerr<<"Query file does not exist. Abort."<<endl;
		MPI::COMM_WORLD.Abort(1);
	}

#ifdef F_GUIDED
	BuildGuidedQueryChunkPointers(queries_filename,pointer_list,size_list,parameter);
	parameter.queries_chunk_size = numeric_limits<unsigned int>::max();
#else
	BuildQueryChunkPointers(queries_filename,pointer_list,size_list,parameter);
#endif
	int query_chunk_number = pointer_list.size();
	for(int i=0;i<query_chunk_number;i++){
		QueryResource query_resource;
		query_resource.size=size_list[i];
		query_resource.ptr=pointer_list[i];
		query_resource.chunk_id=i;
		query_resource.available=false;
		resources.query_list.push_back(query_resource);
	}
	int query_chunk_size=resources.query_list.size();
	MPI::COMM_WORLD.Bcast(&query_chunk_size,1,MPI::INT,0);
	for(int i=0;i<query_chunk_number;i++){
		MPI::COMM_WORLD.Bcast((char*)&pointer_list[i],sizeof(pointer_list[i]),MPI::CHAR,0);
		MPI::COMM_WORLD.Bcast((char*)&size_list[i],sizeof(size_list[i]),MPI::CHAR,0);
		
	}
		
	for(int i=0;i<resources.query_list.size();i++){
		MPIResource::LoadQueryResource(resources,i);
		
	}   //loading all query on master memory

	
}

void MPICommon::SetupQueryResourcesWorker(WorkerResources &resources,AligningParameters &parameter,MPIParameter &mpi_parameter){
	int query_chunk_size;
	MPI::COMM_WORLD.Bcast(&query_chunk_size,1,MPI::INT,0);
	vector<uint64_t> pointer_list;
	vector<uint64_t> size_list;
	pointer_list.resize(query_chunk_size);
	size_list.resize(query_chunk_size);
#ifdef F_GUIDED
	parameter.queries_chunk_size = numeric_limits<unsigned int>::max();
#endif	
	
	for(int i=0;i<query_chunk_size;i++){
		MPI::COMM_WORLD.Bcast((char*)&pointer_list[i],sizeof(pointer_list[i]),MPI::CHAR,0);
		MPI::COMM_WORLD.Bcast((char*)&size_list[i],sizeof(size_list[i]),MPI::CHAR,0);
	}
	
	//setup query
	for(int i=0;i<query_chunk_size;i++){
		QueryResource query;
		query.size=size_list[i];
		query.ptr=pointer_list[i];
		query.chunk_id=i;
		query.available=false;
		resources.query_list.push_back(query);
	}

}
void MPICommon::BuildQueryChunkPointers(string &queries_filename,vector<uint64_t> &chunk_pointer_list,
										vector<uint64_t> &chunk_size_list,AligningParameters &parameter){
	
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
	
		uint64_t size=isf_ptr-begin;
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

void MPICommon::BuildGuidedQueryChunkPointers(string &queries_filename, vector<uint64_t> &chunk_pointer_list, 
											  vector<uint64_t> &chunk_size_list,AligningParameters &parameter){
	
	float switch_ptr_ratio  = 0.8;
	float switch_size_ratio = 0.5;
	ifstream in(queries_filename.c_str());
	uint64_t begin_ptr,end_ptr,file_size,current_ptr = 0;
	int chunk_id = 0;
	uint64_t chunk_size = parameter.queries_chunk_size;
	uint64_t switch_chunk_size = chunk_size * switch_size_ratio;
	uint64_t switch_ptr;
	
	in.seekg(0,ios::end);
	end_ptr= in.tellg();
	in.clear();
	in.seekg(0,ios::beg);
	begin_ptr=in.tellg();
	file_size=end_ptr - begin_ptr;
	switch_ptr = file_size * switch_ptr_ratio;

	chunk_pointer_list.push_back(0);

    while(file_size> current_ptr){
        if(current_ptr>switch_ptr){
            chunk_size=switch_chunk_size;
        }
        in.seekg(current_ptr+chunk_size);
        while(in){
            char c;
            in.get(c);
            if(c=='>'){break;}
        }
        uint64_t p=in.tellg();

        if(in){
            chunk_pointer_list.push_back(p-1);
            chunk_id++;
            chunk_size_list.push_back(chunk_pointer_list[chunk_id]-chunk_pointer_list[chunk_id-1]);
            current_ptr+=chunk_size;
        }else{
            chunk_size_list.push_back(file_size-chunk_pointer_list[chunk_id]);
            current_ptr+=chunk_size;
        }
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





