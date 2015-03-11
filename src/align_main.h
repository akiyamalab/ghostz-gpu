/*
 * align_main.h
 *
 *  Created on: 2013/10/24
 *      Author: shu
 */

#ifndef ALIGN_MAIN_H_
#define ALIGN_MAIN_H_

#include <string>
#include "aligner.h"
#include "aligner_gpu.h"

class AlignMain {
public:
	AlignMain();
	virtual ~AlignMain();
	int Run(int argc, char* argv[]);

private:
	bool BuildParameters(int argc, char* argv[], std::string &input_filename,
			std::string &database_filename, std::string &output_filename,
			AlignerCommon::AligningCommonParameters &parameters);
};

#endif /* ALIGN_MAIN_H_ */
