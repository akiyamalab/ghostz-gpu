/*
 * sequence_seg_filter.cpp
 *
 *  Created on: 2013/10/23
 *      Author: shu
 */

#include "../ext/seg/src/seg.h"
#include "sequence_seg_filter.h"


SequenceSegFilter::SequenceSegFilter(int window) : seg_(window) {

}

int SequenceSegFilter::Filter(const std::string &seq, std::string *masked_seq) {
	int ret = seg_.MaskSequence(seq, masked_seq);
	if (ret != 0) {
		return 1;
	}
	return 0;
}

