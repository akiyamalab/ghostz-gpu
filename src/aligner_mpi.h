/*
 * aligner_mpi.h
 *
 *  Created on: 2017/05/17
 *      Author: goto
 */

#ifndef ALIGNER_MPI_H_
#define ALIGNER_MPI_H_

#include <string>
#include <vector>
#include <tr1/memory>
#include "alphabet_coder.h"
#include "edit_blocks.h"
#include "sequence_type.h"
#include "score_matrix.h"
#include "queries.h"
#include "database.h"
#include "aligner_common.h"
#include "aligner_presearch_thread_mpi.h"
#include "aligner_build_results_thread_mpi.h"
#include "mpi_common.h"
#include "mpi_resource.h"
class AlignerMPI {
public:
	typedef AlignerPresearchThread::DatabaseType DatabaseType;
	typedef AlignerCommon::DatabaseCommonParameters<DatabaseType> DatabaseParameters;
	typedef AlignerCommon::AligningCommonParameters AligningParameters;
	typedef AlignerCommon::PresearchedResult PresearchedResult;
	typedef AlignerCommon::Result Result;
	typedef MPICommon::MPIParameter MPIParameter;
	typedef MPIResource::QueryResource QueryResource;
	
	AlignerMPI() {
	}
	virtual ~AlignerMPI() {
	}

	void Search(QueryResource &query_resource, DatabaseType &database,
				std::vector<std::vector<Result> > &results_list,
				AligningParameters &parameters,MPIParameter &mpi_parameter);
 
private:
	typedef AlignerCommon::AlignmentPositionLessPosition AlignmentPositionLessPosition;
	typedef AlignerCommon::PresearchedResultGreaterScore PresearchedResultGreaterScore;
	typedef AlignerCommon::ResultGreaterScore ResultGreaterScore;
	typedef AlignerCommon::AlignmentPosition AlignmentPosition;
	typedef AlignerCommon::Coordinate Coordinate;

	

	void SetQueriesData(Queries &queries, AligningParameters &parameters,
			std::vector<int> &ungapped_extension_cutoffs,
			std::vector<int> &gapped_extension_triggers);

	void AddResults(DatabaseType &database, AligningParameters &parameters,
			std::vector<std::vector<AlignmentPosition> > &alignment_start_positions_in_database,
			std::vector<PresearchedResult> &added_results,
			std::vector<PresearchedResult> &results);

	void BuildDatabaseParameters(DatabaseParameters &parameters,
			DatabaseType::Parameters &database_parameters);
	virtual void Presearch(Queries &queries, DatabaseType &database,
			AligningParameters &parameters,
			std::vector<std::vector<PresearchedResult> > &results_list);

	void BuildResults(Queries &queries, DatabaseType &database,
			AligningParameters &parameters,
			std::vector<std::vector<PresearchedResult> > &presearch_results_list,
			std::vector<std::vector<Result> > &results_list);

	void WriteOutput(std::ostream &os, Queries &queries, DatabaseType &database,
			AligningParameters &parameters,
			std::vector<std::vector<Result> > &results);

};

#endif /* ALIGNER_H_ */
