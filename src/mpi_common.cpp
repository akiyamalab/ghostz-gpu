/*
 * mpi_common.cpp
 *  Created on 2017/05/16
 *  Author:goto
 */

#include"mpi_common.h"
#include"aligner_mpi.h"

#include"mpi.h"
#include<iostream>
#include<string>
using namespace std;

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

	
}
void MPICommon::RunWorker(AligningParameters &parameter,MPIParameter &mpi_parameter){

}


void MPICommon::RunWorkerGPU(AligningParameters &parameter,MPIParameter &mpi_parameter){

}
					  

