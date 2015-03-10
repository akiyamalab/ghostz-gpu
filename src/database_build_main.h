/*
 * database_build_main.h
 *
 *  Created on: 2013/10/24
 *      Author: shu
 */

#ifndef DATABASE_BUILD_MAIN_H_
#define DATABASE_BUILD_MAIN_H_

#include <string>
#include "seed_searcher.h"
#include "aligner.h"

class DatabaseBuildMain {
public:
	typedef Database<SeedSearcher> DatabaseType;
	DatabaseBuildMain();
	virtual ~DatabaseBuildMain();
	int Run(int argc, char* argv[]);

private:
	bool BuildParameters(int argc, char* argv[], std::string &input_filename,
			std::string &database_filename, Aligner::DatabaseParameters &parameters);
};

#endif /* DATABASE_BUILD_MAIN_H_ */
