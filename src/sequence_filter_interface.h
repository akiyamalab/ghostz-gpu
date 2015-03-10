/*
 * sequence_filter_interface.h
 *
 *  Created on: 2013/10/23
 *      Author: shu
 */

#ifndef SEQUENCE_FILTER_INTERFACE_H_
#define SEQUENCE_FILTER_INTERFACE_H_

#include <string>

class SequenceFilterInterface {
public:

	virtual ~SequenceFilterInterface() {
	};
	virtual int Filter(const std::string &seq, std::string *masked_seq) = 0;
};

#endif /* SEQUENCE_FILTER_INTERFACE_H_ */
