#ifndef  CLASTER_H 
#define  CLASTER_H

#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include "targetver.h"

#include <stdio.h>
#include <stdlib.h>
#include <tchar.h>
#include "mpi.h"
#include "Timer.h"

#define ROW_A 29
#define ROW_A2 28
#define COLUMN_B 15
#define COLUMN_B2 14
#define COLUMNMATRIX 264

#define ALLOCATE_MEMORY_A(ROW) \
	subA = new int* [ROW]; \
	for (int i = 0; i < ROW; ++i) \
		subA[i] = new int[ColumnMatrix]

#define ALLOCATE_MEMORY_B() \
	subB = new int* [COLUMNMATRIX]; \
	for (int i = 0; i < COLUMNMATRIX; ++i) \
		subB[i] = new int[COLUMN_B]

#define ALLOCATE_MEMORY_temp(COLUMN) \
	temp = new int* [COLUMNMATRIX]; \
	for (int i = 0; i < COLUMNMATRIX; ++i) \
		temp[i] = new int[COLUMN_B]

#define ALLOCATE_MEMORY_R(ROW) \
	res = new int* [ROW]; \
	for (int i = 0; i < ROW; ++i) { \
		res[i] = new int[COLUMN_B]; \
		for(int j = 0; j < COLUMN_B; ++j) \
			res[i][j] = 0; \
	}

enum processors {
	PROC0, PROC1, PROC2, PROC3,
	PROC4, PROC5, PROC6, PROC7
};

class CClaster
{
private:
	int Row_A = ROW_A;
	int Row_A2 = ROW_A2;
	int Column_B = COLUMN_B;
	int Column_B2 = COLUMN_B2;
	int ColumnMatrix = COLUMNMATRIX;

	int **storage;

public:
	CClaster();
	CClaster(processors proc);

	void read_memoryA(processors proc, int **tmp);
	void read_memoryB(processors proc, int **tmp);
	void mul(int **a, int **b, int **res, int SUB_N1, int SUB_N3);
	void tmp_storage(int **tempRes, processors proc, int colPos, bool entry);
	//void printTmpMul(int** a, int** b, int** res);
	
	void process_0();
	void process_1();
	void process_2();
	void process_3();
	void process_4();
	void process_5();
	void process_6();
	void process_7();
	
};
#endif 
