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

#include "aligner_common.h"
#include "mpi_resource.h"

class ResultSummarizer{
 public:
	typedef AlignerCommon::Result Result;
	typedef MPIResource::AlignmentTask AlignmentTask;
	
	ResultSummarizer(std::string &tmpDirName);
		
		
	
	void SaveResultFile(std::vector<std::vector<AlignerCommon::Result> > &results_list,
						AlignmentTask task);
	
	int LoadResultFile(std::vector<std::vector<Result> > &results_list,
					   AlignmentTask task);
	
	void SendResult(char *data,int size,AlignmentTask task,int target_rank);
	void RecvResult(std::vector<std::vector<Result> > &results_list,
					AlignmentTask &task);


	void SerializeResult(std::vector<std::vector<Result> > &results_list,
						 char **ptr,int *size);
	void DeserializeResult(std::vector<std::vector<Result> > &results_list,
						   char *ptr,int size);
	
 private:
	std::vector<std::vector<AlignmentTask> > list_;
	std::string tmpDirName_;
	
	
	 
	
	std::string GetTmpFilename(AlignmentTask task){
		std::stringstream ss;
		ss<<tmpDirName_<<"/result"<<task.query_chunk<<"_"<<task.database_chunk;
		return ss.str();
	
	}

	
};



#endif
