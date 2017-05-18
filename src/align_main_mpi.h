/*
 * align_main_mpi.h
 *
 *  Created on: 2017/05/16
 *      Author: goto
 */

#ifndef ALIGN_MAIN_MPI_H_
#define ALIGN_MAIN_MPI_H_

#include <string>
#include "aligner.h"
#include "aligner_gpu.h"
#include "mpi_common.h"

class AlignMainMPI {
public:
	typedef MPICommon::MPIParameter MPIParameter;
	

	AlignMainMPI();
	virtual ~AlignMainMPI();
	int Run(int argc, char* argv[],MPIParameter &mpi_paramter);

private:
	
	bool BuildParameters(int argc, char* argv[], std::string &input_filename,
			std::string &database_filename, std::string &output_filename,
			AlignerCommon::AligningCommonParameters &parameters);
};

#endif /* ALIGN_MAIN_MPI_H_ */
