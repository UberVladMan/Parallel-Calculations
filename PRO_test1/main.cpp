#define _CRT_SECURE_NO_WARNINGS

#include "claster.h"
#include "mpi.h"

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <Windows.h>
#include "Timer.h"

int ProcNum, ProcRank, RecvRank;

int main(int argc, char* argv[]) {
	CClaster obj;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	printf("\n------------------------------");

	if (ProcRank == 0) {
		obj.process_0();
		MPI_Finalize();
	}
	if (ProcRank == 1) {
		obj.process_1();
		MPI_Finalize();
	}
	if (ProcRank == 2) {
		obj.process_2();
		MPI_Finalize();
	}
	if (ProcRank == 3) {
		obj.process_3();
		MPI_Finalize();
	}
	if (ProcRank == 4) {
		obj.process_4();
		MPI_Finalize();
	}
	if (ProcRank == 5) {
		obj.process_5();
		MPI_Finalize();
	}
	if (ProcRank == 6) {
		obj.process_6();
		MPI_Finalize();
	}
	if (ProcRank == 7) {
		obj.process_7();
		MPI_Finalize();
	}
	return 0;
}
