/*
 * gpu_stream_contoller.h
 *
 *  Created on: Aug 6, 2014
 *      Author: shu
 */

#ifndef GPU_STREAM_CONTOLLER_H_
#define GPU_STREAM_CONTOLLER_H_

#include <vector>
#include <queue>
#include <boost/thread/recursive_mutex.hpp>
#include <boost/thread.hpp>
#include "seed_searcher_gpu.h"
#include "device_seeds_memory.h"
#include "gpu_stream_runner.h"
#include "aligner_common.h"

class DistanceCalculatorGpu;
class UngappedExtenderGpu;
class GappedExtenderGpu;

class GpuStreamController {
public:
	static const size_t kMaxNumberOfConcurrentExecutions = 2;
	struct DistanceCalculationSeedsParameters {
		size_t size;
	};

	struct UngappedExtensionWithTriggerSeedsParameters {
		size_t query_seed_size;
		size_t size;
	};

	struct UngappedExtensionSeedsParameters {
		size_t reverse_offset;
		size_t reverse_length;
		size_t forward_offset;
		size_t forward_length;
	};

	struct GappedExtensionSeesParameters {
		size_t reverse_offset;
		size_t reverse_length;
		size_t forward_offset;
		size_t forward_length;
	};

	GpuStreamController();
	virtual ~GpuStreamController();

	void SetDistanceCalculator(DistanceCalculatorGpu *distance_calculator) {
		for (size_t i = 0; i < kMaxNumberOfConcurrentExecutions; ++i) {
			runners_[i].SetDistanceCalculator(distance_calculator);
		}
	}

	void SetUngappedExtender(UngappedExtenderGpu *ungapped_extender) {
		for (size_t i = 0; i < kMaxNumberOfConcurrentExecutions; ++i) {
			runners_[i].SetUngappedExtender(ungapped_extender);
		}
	}

	void SetGappedExtender(GappedExtenderGpu *gapped_extender) {
		for (size_t i = 0; i < kMaxNumberOfConcurrentExecutions; ++i) {
			runners_[i].SetGappedExtender(gapped_extender);
		}
	}

	int Init(size_t device_buffer_size);
	bool FinishedRun(int resource_id);

	void WaitRun();
	void WaitRun(int resource_id);

	int Run(DistanceCalculationSeedsParameters &parameters,
			HostSeedsMemory &host_memory);
	int Run(UngappedExtensionWithTriggerSeedsParameters &parameters,
			HostSeedsMemory &host_memory);
	int Run(UngappedExtensionSeedsParameters &parameters,
			HostSeedsMemory &host_memory);
	int Run(GappedExtensionSeesParameters &parameters,
			HostSeedsMemory &host_memory);

	void ReleaseResource(int resource_id) {
		assert(
				0 <= resource_id
						&& resource_id < kMaxNumberOfConcurrentExecutions);
		available_resource_ids_.push(resource_id);
	}

	uint32_t GetNumberOfAvailableResources() {
		return available_resource_ids_.size();
	}

private:
	std::queue<int> available_resource_ids_;
	DeviceSeedsMemory device_memory_[kMaxNumberOfConcurrentExecutions];
	GpuStreamRunner runners_[kMaxNumberOfConcurrentExecutions];
};

#endif /* GPU_STREAM_CONTOLLER_H_ */
