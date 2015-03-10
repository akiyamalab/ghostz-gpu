/*
 * group_loader.h
 *
 *  Created on: Aug 21, 2014
 *      Author: shu
 */

#pragma once

#include <assert.h>
#include <iterator>

template<typename TInputIterator, int TLength>
class GroupLoader {
public:
	typedef TInputIterator InputIterator;
	typedef typename std::iterator_traits<TInputIterator>::value_type Value;
	static const int kLength = TLength;

	static __forceinline__
	               __host__        __device__ size_t GetTotalSharedMemorySize(
			size_t number_threads) {
		assert((number_threads % kLength) == 0);
		return sizeof(Value)
				* (kLength * number_threads);
	}

	__forceinline__
	__device__ int GetGroupId() {
		int idx = blockIdx.x * blockDim.x + threadIdx.x;
		return idx / kLength;
	}

	__forceinline__
	__device__ int GetIdInGroup() {
		int idx = blockIdx.x * blockDim.x + threadIdx.x;
		return idx - (idx / kLength) * kLength;
	}

	__forceinline__
	__device__ int GetNumberMembers() {
		return kLength;
	}

	__forceinline__
	__device__ int GetNumberGroups() {
		return (blockDim.x * gridDim.x) / kLength;
	}

	template<typename TPositionIterator, typename TInputOutputIterator,
			typename TOutputIterator>
	__forceinline__
	__device__ void Load(InputIterator inputs, TPositionIterator positions,
			int number_positins, TInputOutputIterator shared_temp,
			TOutputIterator thread_data) {
		typedef typename std::iterator_traits<TPositionIterator>::value_type Position;
		int idx = blockIdx.x * blockDim.x + threadIdx.x;
		int group_idx = GetGroupId();
		int group_offset = idx - group_idx * kLength;
		TInputOutputIterator shared_tmp_group = &shared_temp[(threadIdx.x
				/ kLength) * (kLength * kLength)];
		for (int work_i = 0; work_i < number_positins; ++work_i) {
			Position work_p = positions[work_i];
			shared_tmp_group[work_i * kLength + group_offset] = inputs[work_p
					+ group_offset];
		}

		//__syncthreads();
		int work_start = group_offset * kLength;
		if (group_offset < number_positins) {
			for (int i = 0; i < kLength; ++i) {
				thread_data[i] = shared_tmp_group[work_start + i];
			}
		}
		return;
	}

	template<typename TPositionIterator, typename TInputOutputIterator>
	__forceinline__
	__device__ TInputOutputIterator Load(InputIterator inputs,
			TPositionIterator positions, int number_positins,
			TInputOutputIterator shared_output) {
		typedef typename std::iterator_traits<TPositionIterator>::value_type Position;
		int idx = blockIdx.x * blockDim.x + threadIdx.x;
		int group_idx = GetGroupId();
		int group_offset = idx - group_idx * kLength;
		TInputOutputIterator shared_tmp_group = &shared_output[(threadIdx.x
				/ kLength) * kLength * kLength];
		for (int work_i = 0; work_i < number_positins; ++work_i) {
			Position work_p = positions[work_i];
			shared_tmp_group[work_i * kLength + group_offset] = inputs[work_p
					+ group_offset];
		}

		//__syncthreads();
		int work_start = group_offset * kLength;
		return &shared_tmp_group[work_start];
	}
};

