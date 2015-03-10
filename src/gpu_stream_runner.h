/*
 * gpu_stream_runner.h
 *
 *  Created on: Aug 6, 2014
 *      Author: shu
 */

#ifndef GPU_STREAM_RUNNER_H_
#define GPU_STREAM_RUNNER_H_

#include <cuda_runtime_api.h>
#include "ungapped_extender_gpu.h"
#include "distance_calculator_gpu.h"
#include "gapped_extender_gpu.h"
#include "host_seeds_memory.h"
#include "device_seeds_memory.h"

class GpuStreamRunner {
public:
	GpuStreamRunner();
	virtual ~GpuStreamRunner();

	void SetDistanceCalculator(DistanceCalculatorGpu *distance_calculator) {
		distance_calculator_ = distance_calculator;
	}

	void SetUngappedExtender(UngappedExtenderGpu *ungapped_extender) {
		ungapped_extender_ = ungapped_extender;
	}

	void SetGappedExtender(GappedExtenderGpu *gapped_extender) {
		gapped_extender_ = gapped_extender;
	}

	bool Finished();
	void Synchronize();

	int RunDistanceCalculationAsync(HostSeedsMemory &host_buffer,
			DeviceSeedsMemory &device_seed_buffer, size_t size);

	int RunUngappedExtensionWithTriggerAsync(HostSeedsMemory &host_buffer,
			DeviceSeedsMemory &device_buffer, size_t query_seed_size,
			size_t size);

	int RunGappedExtensionAsync(HostSeedsMemory &host_buffer,
			DeviceSeedsMemory &device_buffer, size_t foward_offset,
			size_t foward_size, size_t reverse_offset, size_t reverse_size);

private:
	GpuStreamRunner& operator =(const GpuStreamRunner& rhs);
	GpuStreamRunner(const GpuStreamRunner& rhs);

	static const size_t kNumberStream = 2;

	cudaStream_t stream_[kNumberStream];
	DistanceCalculatorGpu *distance_calculator_;
	UngappedExtenderGpu *ungapped_extender_;
	GappedExtenderGpu *gapped_extender_;

};

#endif /* GPU_STREAM_RUNNER_H_ */
