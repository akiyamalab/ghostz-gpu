/*
 * result_summarizer.cpp
 *    Created on 2017/07/06
 *        Author: goto
 */

#include "result_summarizer.h"
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
}

int ResultSummarizer::LoadResultFile(vector<vector<Result> > &result_list,AlignmentTask task){
	struct stat st;
	if(stat(GetTmpFilename(task).c_str(),&st)==-1){
		return 1;
	}
	ifstream ifs(GetTmpFilename(task).c_str());
	char *data;
	int size;
	int begin,end;
	ifs.seekg(0,ios::end);
	end = ifs.tellg();
	ifs.clear();
	ifs.seekg(0,ios::beg);
	begin=ifs.tellg();
	size=end-begin;
	data = new char[size];
	ifs.read(data,size);
	ifs.close();
	DeserializeResult(result_list,data,size);
	delete [] data;
	   
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
	task.query_chunk[2];
	data = new char[size]; 
	int target_rank=status.Get_source();
	MPI::COMM_WORLD.Recv(data,size,MPI::CHAR,target_rank,0);
	DeserializeResult(results_list,data,size);

}
