/* 
 *  mpi_resource.cpp
 *     Created on 2017/06/17
 *      Author:goto
 */



#include "mpi_resource.h"


using namespace std;


int MPIResource::BcastDatabase(DatabaseResource &database,MPI::Intercomm comm,int root){

	comm.Bcast((char *)&database.inf_size,sizeof(uint64_t),MPI::CHAR,0);
	comm.Bcast((char *)&database.nam_size,sizeof(uint64_t),MPI::CHAR,0);
	comm.Bcast((char *)&database.off_size,sizeof(uint64_t),MPI::CHAR,0);
	comm.Bcast((char *)&database.seq_size,sizeof(uint64_t),MPI::CHAR,0);
	comm.Bcast((char *)&database.scp_size,sizeof(uint64_t),MPI::CHAR,0);
	comm.Bcast((char *)&database.sdp_size,sizeof(uint64_t),MPI::CHAR,0);
	
	database.inf = new char[database.inf_size];
	database.nam = new char[database.nam_size];
	database.off = new char[database.off_size];
	database.seq = new char[database.seq_size];
	database.scp = new char[database.scp_size];
	database.sdp = new char[database.sdp_size];
	
	comm.Bcast(database.inf,database.inf_size,MPI::CHAR,0);
	comm.Bcast(database.nam,database.nam_size,MPI::CHAR,0); 
	comm.Bcast(database.off,database.off_size,MPI::CHAR,0);
	comm.Bcast(database.seq,database.seq_size,MPI::CHAR,0);
	comm.Bcast(database.scp,database.scp_size,MPI::CHAR,0);
	comm.Bcast(database.sdp,database.sdp_size,MPI::CHAR,0);
	
	return 0;
}

