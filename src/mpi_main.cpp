/*
 * mpi_main.cpp
 * 
 *   Created on :2017/05/14
 *      Auther: goto
 */


#include<iostream>
#include"mpi.h"
using namespace std;
#ifdef F_MPI
int main(int argc, char* argv[]){
	MPI_Init(&argc,&argv);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	cout<<"rank"<<rank<<endl;
	MPI_Finalize();
	return 0;
}
#endif

