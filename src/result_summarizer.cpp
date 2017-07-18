/*
 * result_summarizer.cpp
 *    Created on 2017/07/06
 *        Author: goto
 */

#include "result_summarizer.h"
#include <iostream>



#include <sys/stat.h>
using namespace std;


ResultSummarizer::ResultSummarizer(string &tmpDirName):tmpDirName_(tmpDirName){
	
}


void ResultSummarizer::SaveResultFile(vector<vector<Result> > &results_list,AlignmentTask task){
	ofstream ofs(GetTmpFilename(task).c_str());
	char *data;
	int size;
	SerializeResult(results_list,&data,&size);
	ofs.write(data,size);
	ofs.close();
	delete [] data;
	AddList(task,size);
}

int ResultSummarizer::LoadResultFile(char **ptr, int *size,AlignmentTask task){
	struct stat st;
	if(stat(GetTmpFilename(task).c_str(),&st)==-1){
		return 1;
	}
	ifstream ifs(GetTmpFilename(task).c_str());
	int begin,end;
	ifs.seekg(0,ios::end);
	end = ifs.tellg();
	ifs.clear();
	ifs.seekg(0,ios::beg);
	begin=ifs.tellg();
	*size=end-begin;
	*ptr = new char[*size];
	ifs.read(*ptr,*size);
	ifs.close();

	return 0;
}
int ResultSummarizer::LoadResultFile(vector<vector<Result> > &result_list,AlignmentTask task){	
	char *data;
	int size;
	
	if(LoadResultFile(&data,&size,task)!=0){
		return 1;
	}
	
	DeserializeResult(result_list,data,size);
	delete [] data;
	return 0;   
}

void ResultSummarizer::SerializeResult(vector<vector<Result> > &results_list,char **ptr,int *size){
    stringstream ss;

    int n_results = results_list.size();
    ss.write((char *)&n_results,sizeof(int));

    for(int i=0;i<n_results;i++){
        int n_result = results_list[i].size();
        ss.write((char *)&n_result,sizeof(int));

        for(int j=0;j<n_result;j++){
            Result r =results_list[i][j];
            const char *subject_name = r.subject_name.c_str();
            int subject_name_size = r.subject_name.size();
            const char *query_name = r.query_name.c_str();
			int query_name_size = r.query_name.size();
			

			ss.write((char *)&r.database_chunk_id,sizeof(r.database_chunk_id));
            ss.write((char *)&r.subject_id_in_chunk,sizeof(r.subject_id_in_chunk));
            ss.write((char *)&subject_name_size,sizeof(subject_name_size));
            ss.write((char *)subject_name,subject_name_size);
            ss.write((char *)&r.score,sizeof(r.score));
            ss.write((char *)&r.alignment_length,sizeof(r.alignment_length));
            ss.write((char *)&r.mismatches,sizeof(r.mismatches));
            ss.write((char *)&r.gap_openings,sizeof(r.gap_openings));
            ss.write((char *)&r.gap_extensions,sizeof(r.gap_extensions));
            ss.write((char *)&r.hit.query_position,sizeof(r.hit.query_position));
            ss.write((char *)&r.hit.database_position,sizeof(r.hit.database_position));
            ss.write((char *)&r.start.query_position,sizeof(r.start.query_position));
            ss.write((char *)&r.start.database_position,sizeof(r.start.database_position));
            ss.write((char *)&r.end.query_position,sizeof(r.end.query_position));
            ss.write((char *)&r.end.database_position,sizeof(r.end.database_position));
            ss.write((char *)&r.hit_count,sizeof(r.hit_count));
			ss.write((char *)&query_name_size,sizeof(query_name_size));
            ss.write((char *)query_name,query_name_size);
			ss.write((char *)&r.query_length,sizeof(r.query_length));
        }

    }
    ss.seekp(0,ios::end);
    *size=ss.tellp();

    *ptr=new char[*size];
    ss.read(*ptr,*size);

}



void ResultSummarizer::DeserializeResult(vector<vector<Result> > &results_list,char *ptr,int size){
    stringstream ss;
    ss.write(ptr,size);
    int n_results;
    ss.read((char *)&n_results,sizeof(int));
	results_list.resize(n_results);
    for(int i=0;i<n_results;i++){
        int n_result ;
        ss.read((char *)&n_result,sizeof(int));


        for(int j=0;j<n_result;j++){
            Result r;
            int subject_name_size;
            char *subject_name;
			int query_name_size;
			char *query_name;
            ss.read((char *)&r.database_chunk_id,sizeof(r.database_chunk_id));
            ss.read((char *)&r.subject_id_in_chunk,sizeof(r.subject_id_in_chunk));
            ss.read((char *)&subject_name_size,sizeof(subject_name_size));
            subject_name= new char[subject_name_size];
            ss.read(subject_name,subject_name_size);
            r.subject_name = string(subject_name,subject_name+subject_name_size);
            ss.read((char *)&r.score,sizeof(r.score));
            ss.read((char *)&r.alignment_length,sizeof(r.alignment_length));
            ss.read((char *)&r.mismatches,sizeof(r.mismatches));
            ss.read((char *)&r.gap_openings,sizeof(r.gap_openings));
            ss.read((char *)&r.gap_extensions,sizeof(r.gap_extensions));
            ss.read((char *)&r.hit.query_position,sizeof(r.hit.query_position));
            ss.read((char *)&r.hit.database_position,sizeof(r.hit.database_position));
            ss.read((char *)&r.start.query_position,sizeof(r.start.query_position));
            ss.read((char *)&r.start.database_position,sizeof(r.start.database_position));
            ss.read((char *)&r.end.query_position,sizeof(r.end.query_position));
            ss.read((char *)&r.end.database_position,sizeof(r.end.database_position));
            ss.read((char *)&r.hit_count,sizeof(r.hit_count));
			
			ss.read((char *)&query_name_size,sizeof(query_name_size));
            query_name= new char[query_name_size];
            ss.read(query_name,query_name_size);
            r.query_name = string(query_name,query_name+query_name_size);
			ss.read((char *)&r.query_length,sizeof(r.query_length));
            results_list[i].push_back(r);

        }
    }
}


void ResultSummarizer::SendResult(char *data,int size,AlignmentTask task,int target_rank){
	int buf[3];
	buf[0]=size;
	buf[1]=task.database_chunk;
	buf[2]=task.query_chunk;
	
	MPI::COMM_WORLD.Send(buf,sizeof(buf),MPI::INT,target_rank,0);
	MPI::COMM_WORLD.Send(data,size,MPI::CHAR,target_rank,0);
	
}

void ResultSummarizer::RecvResult(vector<vector<Result> > &results_list,AlignmentTask &task){
	int size;
	int buf[3];
	char *data;
	MPI::Status status;
	MPI::COMM_WORLD.Recv(buf,sizeof(buf),MPI::INT,MPI::ANY_SOURCE,MPI::ANY_TAG,status);
	size=buf[0];
	task.database_chunk=buf[1];
	task.query_chunk=buf[2];
	data = new char[size]; 
	int target_rank=status.Get_source();
	MPI::COMM_WORLD.Recv(data,size,MPI::CHAR,target_rank,0);
	DeserializeResult(results_list,data,size);

}

void ResultSummarizer::AddList(AlignmentTask task, int size){
	task_list.push_back(task);
	size_list.push_back(size);

}
void ResultSummarizer::GatherResultMaster(int query_chunk_size, int database_chunk_size){
	//cout<<"query_size:"<<query_chunk_size<<endl;
	BuildResultTargetMap(query_chunk_size,MPI::COMM_WORLD.Get_size());
	int map_size = result_target_map.size();
	
		result_size_map.resize(query_chunk_size);
		
	for(int i=0;i<result_size_map.size();i++){
		result_size_map[i].resize(database_chunk_size);
	}
	
	for(int i=0;i<query_chunk_size * database_chunk_size;i++){
		//MPI::Status status;
		int buf[3];
		MPI::COMM_WORLD.Recv(buf,3,MPI::INT,MPI::ANY_SOURCE,MPI::ANY_TAG);
		//int src_rank = status.Get_source();
		int query_chunk = buf[0];
		int db_chunk=buf[1];
		int result_size = buf[2];
		result_size_map[query_chunk][db_chunk]=result_size;
		
	}

	MPI::COMM_WORLD.Bcast(&map_size,1,MPI::INT,0);
	MPI::COMM_WORLD.Bcast(&(result_target_map[0]),map_size,MPI::INT,0);
	
	vector<int> result_size_map_vec;
	for(int q=0;q<query_chunk_size;q++){
		for(int d=0;d<database_chunk_size;d++){
			result_size_map_vec.push_back(result_size_map[q][d]);
		}
	}
	int result_size_map_vec_length = result_size_map_vec.size();
	MPI::COMM_WORLD.Bcast(&result_size_map_vec_length,1,MPI::INT,0);
	MPI::COMM_WORLD.Bcast(&result_size_map_vec[0],result_size_map_vec.size(),MPI::INT,0);
 
	
	ReduceResult(0,query_chunk_size,database_chunk_size);

	
}
void ResultSummarizer::BuildResultTargetMap(int query_chunk_size,int mpi_size){
	result_target_map.resize(query_chunk_size);
	
	for(int i=0;i<result_target_map.size();i++){
		result_target_map[i]=i%mpi_size;
	}

}
void ResultSummarizer::GatherResultWorker(int rank){
	stringstream ss;
   	ss<<"rank"<<rank;
	for(int i=0;i<task_list.size();i++){
		int buf[3];
		buf[0]=task_list[i].query_chunk;
		buf[1]=task_list[i].database_chunk;
		buf[2]=size_list[i];
		MPI::COMM_WORLD.Send(buf,3,MPI::INT,0,0);
		ss<<"("<<buf[0]<<","<<buf[1]<<","<<buf[2]<<") ";		
	}
	//cout<<ss.str()<<endl;

	int query_chunk_size;
	MPI::COMM_WORLD.Bcast(&query_chunk_size,1,MPI::INT,0);
	result_target_map.resize(query_chunk_size);
	//	cout<<"map_size:"<<query_chunk_size<<endl;
	int *target_map = new int[query_chunk_size];
 	MPI::COMM_WORLD.Bcast(target_map,query_chunk_size,MPI::INT,0);
	for(int i=0;i<query_chunk_size;i++){
		result_target_map[i]=target_map[i];
	}
	
	int result_size_map_vec_length;
	MPI::COMM_WORLD.Bcast(&result_size_map_vec_length,1,MPI::INT,0);
	vector<int> result_size_map_vec;
	result_size_map_vec.resize(result_size_map_vec_length);
	
	MPI::COMM_WORLD.Bcast(&(result_size_map_vec[0]),result_size_map_vec_length,MPI::INT,0);
	int database_chunk_size = result_size_map_vec_length / query_chunk_size;
	//	cout<<"db_size:"<<database_chunk_size<<endl;

	result_size_map.resize(query_chunk_size);

	for(int q=0;q<query_chunk_size;q++){
		result_size_map[q].resize(database_chunk_size);
		for(int d=0;d<database_chunk_size;d++){
			result_size_map[q][d]=result_size_map_vec[q * database_chunk_size + d];
 
		}
	}
	
	// finish sync of size list
	
	ReduceResult(rank,query_chunk_size,database_chunk_size);
	
	
 	 
		
   	
	
}

void ResultSummarizer::ReduceResult(int rank,int query_chunk_size,int database_chunk_size){
	int count_recv=0;
	vector<int> my_target;
	vector<vector<char *> > result_data_list;
	vector<vector<MPI::Request> > result_data_req_list;
	result_data_list.resize(query_chunk_size);
	result_data_req_list.resize(query_chunk_size);
	for(int i=0;i<query_chunk_size;i++){
		result_data_list[i].resize(database_chunk_size);
		result_data_req_list[i].resize(database_chunk_size);
		
		for(int j=0;j<database_chunk_size;j++){
			result_data_list[i][j]=NULL;		

		}
	}
	
	
	stringstream ss;
	ss.str("");
	
   for(int i=0;i<result_target_map.size();i++){
		if(result_target_map[i]==rank){
			my_target.push_back(i);
			count_recv+=database_chunk_size;
			}
	}
	
	//Send
	vector<char*> ptr_list;
	vector<MPI::Request> req_list;
	for(int i=0;i<task_list.size();i++){
		int q_chunk = task_list[i].query_chunk;
		int db_chunk = task_list[i].database_chunk;
		int size ;
		char *data ;
		if(LoadResultFile(&data,&size,task_list[i])!=0){
			cout<<"Error. Cannnot open tempfile."<<endl;
		};
		int target_rank = result_target_map[q_chunk];

			   
		
		if(target_rank==rank){
			result_data_list[q_chunk][db_chunk]=data;
			count_recv--;
		}else{
			ptr_list.push_back(data);
			MPI::Request req;
			int tag = q_chunk * database_chunk_size + db_chunk;
			//Send 
			cout<<"rank "<<rank<<"("<<q_chunk<<","<<db_chunk<<")"<<" to"<<target_rank<<"  "<<size<<"byte"<<endl;
			req = MPI::COMM_WORLD.Isend(data,size,MPI::CHAR,target_rank,tag);
			
			req_list.push_back(req);
		}
		
	}
	//Recv and Create thread
	boost::thread_group threads;
	
	for(int i=0;i<my_target.size();i++){
		for(int db_chunk=0;db_chunk<database_chunk_size;db_chunk++){
			int q_chunk = my_target[i];
			if(result_data_list[q_chunk][db_chunk]==NULL){
				int size  = result_size_map[q_chunk][db_chunk];
				char *data = new char[size];
				int tag = q_chunk * database_chunk_size + db_chunk;
				result_data_list[q_chunk][db_chunk]=data;
				MPI::Status status;
			   	result_data_req_list[q_chunk][db_chunk] 
					= MPI::COMM_WORLD.Irecv(data,size,MPI::CHAR,MPI::ANY_SOURCE,tag);
				
			}//else already this rank have data.			
		}
		
	}
	ss.str("");
	ss<<"rank:"<<rank<<"  ";
	for(int i=0;i<my_target.size();i++){
		ss<<my_target[i]<<" ";
	}
	ss<<endl;
	for(int i=0;i<query_chunk_size;i++){
		ss<<i<<":";
		for(int j=0;j<database_chunk_size;j++){
	 		result_data_req_list[i][j].Wait();
			if(result_data_list[i][j]==NULL){
				ss<<"F";
			}else{
				ss<<"T";
			}
			
			
		}
		ss<<endl;
	}	
	cout<<ss.str();
	
	int first_target = my_target[0];
	vector<vector<Result> > results_list;
	vector<vector<Result> > results_list_;
	AlignmentTask task;
	task.query_chunk=first_target;
	task.database_chunk=0;
	
   	DeserializeResult(results_list,
					  result_data_list[first_target][0],result_size_map[first_target][0]);
	
	//DeserializeResult(results_list_,
	//				  result_data_list[first_target][1],result_size_map[first_target][1]);
	//cout<<"rank:"<<rank<<" Q:"<<first_target<<":"<<results_list[0].size()<<endl;
   
	for(int i=0;i<results_list[0].size();i++){
		//cout<<i<<endl;
		cout<<results_list[0][i].subject_name<<":"<<results_list[0][0].query_name<<":len "<<results_list[0][0].query_length<<endl;
		}
	
	
		
	
	
		
	cout<<"rank "<<rank<<" clean up. "<<ptr_list.size()<<endl;
	for(int i=0;i<ptr_list.size();i++){
		req_list[i].Wait();
		free(ptr_list[i]);
	}

}
	
void ResultSummarizer::ReduceResultThread(int target_query){
	
	
}
	

