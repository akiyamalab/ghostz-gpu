/*
 * aligner_gpu_presearch_thread_mpi.h
 *
 *  Created on: 2017/10/12
 *      Author: goto
 */

#ifndef ALIGNER_GPU_PRESEARCH_THREAD_MPI_H_
#define ALIGNER_GPU_PRESEARCH_THREAD_MPI_H_

#include "aligner_gpu_presearch_thread.h"

//class GappedExtender;
class AlignerGpuPresearchThread;

class AlignerGpuPresearchThreadMPI : public AlignerGpuPresearchThread
{
 public:
    void Run(ThreadParameters& thread_parameters);

 private:
};

#endif /* ALIGNER_GPU_PRESEARCH_THREAD_MPI_H_ */

