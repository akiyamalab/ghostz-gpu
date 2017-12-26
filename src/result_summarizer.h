/*
 * result_summarizer.h
 *    Created on 2017/07/06
 *        Author: goto
 */

#ifndef __RESULT_SUMMARIZER_H__
#define __RESULT_SUMMARIZER_H__


#include "mpi.h"
#include <string>
#include <sstream>
#include <boost/thread.hpp>
#include <stdint.h>

#include "statistics.h"

#include "aligner_common.h"
#include "mpi_resource.h"

class ResultSummarizer{
 public:
	typedef AlignerCommon::Result Result;
	typedef AlignerCommon::AligningCommonParameters AligningParameters;
	typedef MPIResource::AlignmentTask AlignmentTask;
	typedef MPIResource::DatabaseInfo DatabaseInfo;
	ResultSummarizer(std::string &tmp_dirname,std::string &output_filename);

	struct ThreadParameters{
		int query_chunk;
		MPIResource::DatabaseInfo database_info;
		std::vector<MPI::Request> reqs;
		std::vector<char *> result_data_list;
		std::vector<uint64_t> result_size_list;
		AlignerCommon::AligningCommonParameters parameters;
		std::ostream *os;
		boost::mutex *check_request_mutex;
	};
		
		
	
	void SaveResultFile(std::vector<std::vector<AlignerCommon::Result> > &results_list,
						AlignmentTask task);
	
	int LoadResultFile(std::vector<std::vector<Result> > &results_list,
					   AlignmentTask task);
	int LoadResultFile(char **ptr,uint64_t *size,AlignmentTask task);
	void SendResult(char *data,uint64_t size,AlignmentTask task,int target_rank);
	void RecvResult(std::vector<std::vector<Result> > &results_list,
					AlignmentTask &task);


	void SerializeResult(std::vector<std::vector<Result> > &results_list,
						 char **ptr,uint64_t *size);
	void DeserializeResult(std::vector<std::vector<Result> > &results_list,
						   char *ptr,uint64_t size);

	void GatherResultMaster(int query_chunk_size,int mpi_size,
							AligningParameters &parameters, DatabaseInfo &database_info);
	void GatherResultWorker(int rank,
							AligningParameters &parameters, DatabaseInfo &database_info);
   
	void AddList(AlignmentTask task,uint64_t size);
		
	void BuildResultTargetMap(int query_chunk_size, int mpi_size);
 private:
	std::vector<AlignmentTask > task_list; // searched results file on each process
	std::vector<uint64_t> size_list;            // file size of results file on each process
	std::string tmp_dirname_;
	std::string output_filename_;
	std::vector<int> result_target_map;   // map of query block to report phase node number
 	std::vector<std::vector<uint64_t> > result_size_map; // sizeof (query,db) 's results file
	
	
	void DeleteResultFile(AlignmentTask task);
	void ReduceResult(int rank,int query_chunk_size,int database_chunk_size,
					  AligningParameters &parameters, DatabaseInfo &database_info);
	void CreateReduceResultThread(boost::thread_group threads,int query_chunk,
								  std::vector<char*> result_data_list,
								  std::vector<MPI::Request> result_data_req_list);
	void ReduceResultThread(ThreadParameters &thread_parameters);

	void ReduceResultViaFS(int rank,int query_chunk_size,int database_chunk_size,
						   AligningParameters &parameters, DatabaseInfo &database_info);
	
	void WriteOutput(std::ostream &os,std::string &query_name,
					 uint32_t query_length,
					 DatabaseInfo &database_info,AligningParameters &parameters,
					 std::vector<Result> &result_list);
	
	
	std::string GetTmpFilename(AlignmentTask task){
		std::stringstream ss;
		ss<<tmp_dirname_<<"/result"<<task.query_chunk<<"_"<<task.database_chunk;
		return ss.str();
	
	}

	
};



#endif
