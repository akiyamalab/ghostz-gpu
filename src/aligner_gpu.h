/*
 * aligner_gpu.h
 *
 *  Created on: Jul 10, 2014
 *      Author: shu
 */

#ifndef ALIGNER_GPU_H_
#define ALIGNER_GPU_H_

#include <vector>
#include <boost/thread/barrier.hpp>
#include "aligner.h"
#include "seed_searcher_gpu.h"
#include "aligner_gpu_presearch_thread.h"
#include "gpu_stream_controller.h"

class UngappedExtenderGpu;
class UngappedExtender;
class GappedExtender;
class AlignerGpuData;

class AlignerGpu {
public:
	typedef AlignerGpuPresearchThread::DatabaseType DatabaseType;
	typedef AlignerCommon::DatabaseCommonParameters<DatabaseType> DatabaseParameters;
	typedef AlignerCommon::AligningCommonParameters AligningParameters;
	typedef AlignerCommon::PresearchedResult PresearchedResult;
	typedef AlignerCommon::Result Result;
	AlignerGpu();
	virtual ~AlignerGpu();

	void Align(std::string &queries_filename, std::string &database_filename,
				std::string &output_filename, AligningParameters &parameters);

private:
	typedef AlignerCommon::AlignmentPositionLessPosition AlignmentPositionLessPosition;
	typedef AlignerCommon::PresearchedResultGreaterScore PresearchedResultGreaterScore;
	typedef AlignerCommon::ResultGreaterScore ResultGreaterScore;
	typedef AlignerCommon::AlignmentPosition AlignmentPosition;
	typedef AlignerCommon::Coordinate Coordinate;

	static const uint32_t kMaxSeedBufferSizePerGPU = 1 << 22;

	void Presearch(Queries &queries, DatabaseType &database,
			AligningParameters &parameters,
			std::vector<std::vector<PresearchedResult> > &results_list);

	void SetQueriesData(Queries &queries,
			AlphabetCoder::Code sequence_delimiter,
			AligningParameters &parameters,
			std::vector<AlphabetCoder::Code> &queries_concatenated_sequence,
			std::vector<uint32_t> &query_sequence_offsets,
			std::vector<int> &ungapped_extension_cutoffs,
			std::vector<int> &gapped_extension_triggers);


	void BuildResults(Queries &queries, DatabaseType &database,
			AligningParameters &parameters,
			std::vector<std::vector<PresearchedResult> > &presearch_results_list,
			std::vector<std::vector<Result> > &results_list);
};

#endif /* ALIGNER_GPU_H_ */
