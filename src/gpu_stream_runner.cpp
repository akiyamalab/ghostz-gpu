/*
 * gpu_stream_runner.cpp
 *
 *  Created on: Aug 6, 2014
 *      Author: shu
 */

#include "gpu_stream_runner.h"
#include "host_seeds_memory.h"
#include "device_seeds_memory.h"
#include "cuda_common.h"
#include "ungapped_extender_gpu.h"

GpuStreamRunner::GpuStreamRunner() :
		distance_calculator_(NULL), ungapped_extender_(NULL), gapped_extender_(
				NULL) {
	for (size_t i = 0; i < kNumberStream; ++i) {
		CUDA_CHECK_RETURN(cudaStreamCreate(&stream_[i]));
	}
}

GpuStreamRunner::~GpuStreamRunner() {
	for (size_t i = 0; i < kNumberStream; ++i) {
		CUDA_CHECK_RETURN(cudaStreamDestroy(stream_[i]));
	}
}

bool GpuStreamRunner::Finished() {
	return !(cudaErrorNotReady == cudaStreamQuery(stream_[0])
			|| cudaErrorNotReady == cudaStreamQuery(stream_[1]));
}

void GpuStreamRunner::Synchronize() {
	CUDA_CHECK_RETURN(cudaStreamSynchronize(stream_[0]));
	CUDA_CHECK_RETURN(cudaStreamSynchronize(stream_[1]));
}

int GpuStreamRunner::RunDistanceCalculationAsync(HostSeedsMemory &host_buffer,
		DeviceSeedsMemory &device_buffer, size_t size) {
	distance_calculator_->CalculateDistancesAsync(size,
			host_buffer.GetQueryConcatenatedPositions(),
			host_buffer.GetDatabasePositions(),
			(DistanceCalculatorGpu::Distance*) host_buffer.GetResultValues(),
			device_buffer.GetQueryConcatenatedPositions(),
			device_buffer.GetDatabasePositions(),
			(DistanceCalculatorGpu::Distance*) device_buffer.GetResultValues(),
			stream_[0]);
	return 0;
}

int GpuStreamRunner::RunUngappedExtensionWithTriggerAsync(
		HostSeedsMemory &host_buffer, DeviceSeedsMemory &device_buffer,
		size_t query_seed_size, size_t size) {
	ungapped_extender_->SetQuerySeeds(query_seed_size, size,
			host_buffer.GetQuerySeedIds(), host_buffer.GetQuerySeedStarts(),
			device_buffer.GetTempStorageSize(), device_buffer.GetTempStorage(),
			device_buffer.GetQuerySeedIds(), device_buffer.GetQuerySeedStarts(),
			device_buffer.GetQueryIds(),
			device_buffer.GetQueryConcatenatedPositions(),
			device_buffer.GetDatabasePositions(), stream_[0]);
	ungapped_extender_->ExtendWithTriggerAsync(size, NULL, NULL,
			host_buffer.GetDatabasePositions(),
			(char*) host_buffer.GetResultValues(), device_buffer.GetQueryIds(),
			device_buffer.GetQueryConcatenatedPositions(),
			device_buffer.GetDatabasePositions(),
			(char*) device_buffer.GetResultValues(),
			(int*) device_buffer.GetQuerySeedIds(), stream_[0]);
	return 0;
}

int GpuStreamRunner::RunGappedExtensionAsync(HostSeedsMemory &host_buffer,
		DeviceSeedsMemory &device_buffer, size_t foward_offset,
		size_t foward_size, size_t reverse_offset, size_t reverse_size) {

	gapped_extender_->ConvertToGappedExtensionSeedsAsync(reverse_size, true,
			host_buffer.GetQuerySeedIds() + reverse_offset,
			host_buffer.GetQueryConcatenatedPositions() + reverse_offset,
			host_buffer.GetDatabasePositions() + reverse_offset,
			device_buffer.GetQuerySeedIds() + reverse_offset,
			device_buffer.GetQuerySeedStarts() + reverse_offset,
			device_buffer.GetQueryConcatenatedPositions() + reverse_offset,
			device_buffer.GetDatabasePositions() + reverse_offset, stream_[0]);
	gapped_extender_->ExtendOneSideScoreOnlyAsync(reverse_size, true,
			host_buffer.GetQueryConcatenatedPositions() + reverse_offset,
			host_buffer.GetDatabasePositions() + reverse_offset,
			host_buffer.GetResultValues() + reverse_offset,
			device_buffer.GetQueryConcatenatedPositions() + reverse_offset,
			device_buffer.GetDatabasePositions() + reverse_offset,
			device_buffer.GetResultValues() + reverse_offset, stream_[0]);

	gapped_extender_->ConvertToGappedExtensionSeedsAsync(foward_size, false,
			host_buffer.GetQuerySeedIds() + foward_offset,
			host_buffer.GetQueryConcatenatedPositions() + foward_offset,
			host_buffer.GetDatabasePositions() + foward_offset,
			device_buffer.GetQuerySeedIds() + foward_offset,
			device_buffer.GetQuerySeedStarts() + foward_offset,
			device_buffer.GetQueryConcatenatedPositions() + foward_offset,
			device_buffer.GetDatabasePositions() + foward_offset, stream_[1]);
	gapped_extender_->ExtendOneSideScoreOnlyAsync(foward_size, false,
			host_buffer.GetQueryConcatenatedPositions() + foward_offset,
			host_buffer.GetDatabasePositions() + foward_offset,
			host_buffer.GetResultValues() + foward_offset,
			device_buffer.GetQueryConcatenatedPositions() + foward_offset,
			device_buffer.GetDatabasePositions() + foward_offset,
			device_buffer.GetResultValues() + foward_offset, stream_[1]);
	return 0;
}
