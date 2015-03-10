/*
 * cuda_common.h
 *
 *  Created on: Jul 14, 2014
 *      Author: shu
 */

#pragma once

#include <cstdio>
#include <stdint.h>

#define CUDA_CHECK_RETURN(value) {											\
	cudaError_t _m_cudaStat = value;										\
	if (_m_cudaStat != cudaSuccess) {										\
		fprintf(stderr, "Error %s at line %d in file %s\n",					\
				cudaGetErrorString(_m_cudaStat), __LINE__, __FILE__);		\
		exit(1);															\
	} }

namespace cuda_common {

typedef struct {
	uint32_t query_position;
	uint32_t database_position;
} Coordinate;

enum Direction {
	kFoward, kReverse,
};
static const size_t kMaxLoadLength = 4;
}
