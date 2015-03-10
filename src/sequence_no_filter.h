/*
 * sequence_no_filter.h
 *
 *  Created on: 2013/10/23
 *      Author: shu
 */

#ifndef SEQUENCE_NO_FILTER_H_
#define SEQUENCE_NO_FILTER_H_

#include "sequence_filter_interface.h"

class SequenceNoFilter : public SequenceFilterInterface{
public:
	SequenceNoFilter() {
	}

	~SequenceNoFilter() {

	}
	int Filter(const std::string &seq, std::string *masked_seq);
};

#endif /* SEQUENCE_NO_FILTER_H_ */
