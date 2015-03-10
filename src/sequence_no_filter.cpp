/*
 * sequence_no_filter.cpp
 *
 *  Created on: 2013/10/23
 *      Author: shu
 */

#include "sequence_no_filter.h"


int SequenceNoFilter::Filter(const std::string &seq, std::string *masked_seq) {
	*masked_seq = seq;
	return 0;
}
