/*
 * sequence_seg_filter.h
 *
 *  Created on: 2013/10/23
 *      Author: shu
 */

#ifndef SEQUENCE_SEG_FILTER_H_
#define SEQUENCE_SEG_FILTER_H_

#include "sequence_filter_interface.h"
#include "../ext/seg/src/seg.h"

class SequenceSegFilter: public SequenceFilterInterface {
public:
	SequenceSegFilter(): seg_(kDefalutWindow) {
	}
	SequenceSegFilter(int window);

	~SequenceSegFilter() {
	}
	int Filter(const std::string &seq, std::string *masked_seq);
private:
	static const int kDefalutWindow = 12;
	seg::Seg seg_;
};
#endif /* SEQUENCE_SEG_FILTER_H_ */
