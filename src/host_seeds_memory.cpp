/*
 * host_seed_buffer.cpp
 *
 *  Created on: Aug 5, 2014
 *      Author: shu
 */

#include <cstring>
#include <iostream>
#include "host_seeds_memory.h"

using namespace std;

HostSeedsMemory::HostSeedsMemory() :
		allocated_(false), id_(0), size_(0){
}

HostSeedsMemory::HostSeedsMemory(size_t id, size_t size) :
		allocated_(false), id_(0), size_(0){
	Init(id, size);
}

HostSeedsMemory::~HostSeedsMemory() {
}

int HostSeedsMemory::Init(size_t id, size_t size) {
	allocated_ = false;
	id_ = id;
	size_ = size;
	query_seed_ids_.resize(size);
	query_seed_lengths_.resize(size + 1);
	query_ids_.resize(size);
	query_concatenated_positions_.resize(size);
	database_positions_.resize(size);
	result_values_.resize(size, 0);
	return 0;
}

