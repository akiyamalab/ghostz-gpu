/*
 * gpu_stream_contoller.cpp
 *
 *  Created on: Aug 6, 2014
 *      Author: shu
 */

#include <iostream>
#include "gpu_stream_controller.h"

using namespace std;

GpuStreamController::GpuStreamController() {
}

GpuStreamController::~GpuStreamController() {

}

int GpuStreamController::Init(size_t device_buffer_size) {
	for (size_t i = 0; i < kMaxNumberOfConcurrentExecutions; ++i) {
		device_memory_[i].Init(device_buffer_size);
		available_resource_ids_.push(i);
	}
	return 0;
}

bool GpuStreamController::FinishedRun(int resource_id) {
	if (resource_id >= 0) {
		return runners_[resource_id].Finished();
	} else {
		return true;
	}
}

void GpuStreamController::WaitRun() {
	for (size_t i = 0; i < kMaxNumberOfConcurrentExecutions; ++i) {
		runners_[i].Synchronize();
	}
}

void GpuStreamController::WaitRun(int resource_id) {
	if (resource_id >= 0) {
		runners_[resource_id].Synchronize();
	}
}

int GpuStreamController::Run(DistanceCalculationSeedsParameters &parameters,
		HostSeedsMemory &host_memory) {
	assert(!available_resource_ids_.empty());
	int id = available_resource_ids_.front();
	available_resource_ids_.pop();
	runners_[id].RunDistanceCalculationAsync(host_memory, device_memory_[id],
			parameters.size);
	return id;
}
int GpuStreamController::Run(
		UngappedExtensionWithTriggerSeedsParameters &parameters,
		HostSeedsMemory &host_memory) {
	assert(!available_resource_ids_.empty());
	int id = available_resource_ids_.front();
	available_resource_ids_.pop();
	runners_[id].RunUngappedExtensionWithTriggerAsync(host_memory,
			device_memory_[id], parameters.query_seed_size, parameters.size);
	return id;
}

int GpuStreamController::Run(GappedExtensionSeesParameters &parameters,
		HostSeedsMemory &host_memory) {
	assert(!available_resource_ids_.empty());
	int id = available_resource_ids_.front();
	available_resource_ids_.pop();
	runners_[id].RunGappedExtensionAsync(host_memory, device_memory_[id],
			parameters.forward_offset, parameters.forward_length,
			parameters.reverse_offset, parameters.reverse_length);
	return id;
}
