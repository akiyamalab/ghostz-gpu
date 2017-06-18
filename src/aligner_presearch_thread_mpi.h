/*
 * aligner_presearch_thread_mpi.h
 *
 *  Created on: 2017/06/17
 *      Author: goto
 */

#ifndef ALIGNER_PRESEARCH_THREAD_MPI_H_
#define ALIGNER_PRESEARCH_THREAD_MPI_H_

#include "aligner_presearch_thread.h"

//class GappedExtender;
class AlignerPresearchThread;

class AlignerPresearchThreadMPI : public AlignerPresearchThread 
{
 public:
	void Run(ThreadParameters& thread_parameters);
	
 private:
};

#endif /* ALIGNER_PRESEARCH_THREAD_MPI_H_ */
