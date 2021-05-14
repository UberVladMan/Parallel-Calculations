#include "claster.h"

MPI_Status Status;


CClaster::CClaster() {}

CClaster::CClaster(processors proc) {
	if (proc < 6) {
		storage = new int* [ROW_A];
		for (int i = 0; i < ROW_A; ++i) {
			storage[i] = new int[COLUMN_B * 5 + COLUMN_B2 * 3];
			for (int j = 0; j < COLUMN_B * 5 + COLUMN_B2 * 3; ++j)
				storage[i][j] = 0;
		}
	}
	else {
		storage = new int* [ROW_A2];
		for (int i = 0; i < ROW_A2; ++i) {
			storage[i] = new int[COLUMN_B * 5 + COLUMN_B2 * 3];
			for (int j = 0; j < COLUMN_B * 5 + COLUMN_B2 * 3; ++j)
				storage[i][j] = 0;
		}
	}

}


void CClaster::mul(int** a, int** b, int** res, int SUB_N1, int SUB_N3) {
	for (int row = 0; row < SUB_N1; row++) {
		for (int col = 0; col < SUB_N3; col++) {
			res[row][col] = 0;
			for (int inner = 0; inner < COLUMNMATRIX; inner++) {
				res[row][col] += a[row][inner] * b[inner][col];
			}
		}
	}
}

void CClaster::read_memoryA(processors proc, int **tmp) //читаю з А
{
	int sub_n1 = proc < 6 ? Row_A : Row_A2;
	int rowStart = proc < 6 ? proc * sub_n1 : 6 * Row_A + (proc - 6) * Row_A2;

	FILE* memory;
	memory = fopen("matrixA.txt", "r");

	while (rowStart--) {
		while (fgetc(memory) != '\n');
	}
	for (int row = 0; row < sub_n1; ++row) {
		for (int col = 0; col < ColumnMatrix; ++col)
			fscanf(memory, "%d", &tmp[row][col]);
	}

	fclose(memory);
}

void CClaster::read_memoryB(processors proc, int **tmp) {
	int sub_n3 = proc < 5 ? Column_B : Column_B2;
	int colStart = proc < 5 ? proc * sub_n3 : 5 * Column_B + (proc - 5) * Column_B2;
	
	FILE* memory;
	memory = fopen("matrixB.txt", "r");

	for (int row = 0; row < ColumnMatrix; ++row) {
		for (int i = 0; i < colStart;) {
			if (fgetc(memory) == '\t')
				++i;
		}
		for (int col = 0; col < sub_n3; ++col) {
			fscanf(memory, "%d", &tmp[row][col]);
		}
		while (fgetc(memory) != '\n');
	}

	fclose(memory);
}

void CClaster::tmp_storage(int **tempRes, processors proc, int colPos, bool end)
{
	int sub_n3 = colPos < 5 ? Column_B : Column_B2;
	int lastCol = colPos < 5 ? colPos * sub_n3 : 5 * Column_B + (colPos - 5) * Column_B2;

	for (int row = 0; row < (proc < 6 ? ROW_A : ROW_A2 ); row++) {
		int c = 0;
		for (int col = lastCol; col < lastCol + sub_n3; col++) {
			storage[row][col] = tempRes[row][c];
			++c;
		}
	}

	if (end == true) {
		FILE* memory;
		memory = fopen("Result.txt", proc == 0 ? "w" : "a");
		for (int row = 0; row < (proc < 6 ? ROW_A : ROW_A2); row++) {
			for (int col = 0; col < COLUMN_B * 5 + COLUMN_B2 * 3; col++)
				fprintf(memory, "%d\t", storage[row][col]);
			fputc('\n', memory);
		}
		fclose(memory);
	}
}

/*void CClaster::printTmpMul(int **a, int** b, int **res) {
	printf("\ntmp mul:\n");
	for (int i = 0; i < ROW_A; ++i) {
		for (int j = 0; j < COLUMNMATRIX; ++j) {
			printf("%d\t", a[i][j]);
		}
		printf("\t");
		for (int j = 0; j < COLUMN_B; ++j) {
			printf("%d\t", b[i][j]);
		}
		printf("\t");
		for (int j = 0; j < COLUMN_B; ++j) {
			printf("%d\t", res[i][j]);
		}
		putchar('\n');
	}
	for (int i = ROW_A; i < COLUMNMATRIX; ++i) {
		for (int j = 0; j < COLUMNMATRIX + 1; ++j) {
			putchar('\t');
		}
		for (int j = 0; j < COLUMN_B; ++j) {
			printf("%d\t", b[i][j]);
		}
		putchar('\n');
	}
	putchar('\n');
}*/

void CClaster::process_0() {
	Timer timer;
	double clockBegin = timer.getCurrentTime();

	CClaster obj(PROC0);

	int** subA;
	int** subB;
	int** res;

	ALLOCATE_MEMORY_A(ROW_A);
	ALLOCATE_MEMORY_B();
	ALLOCATE_MEMORY_R(ROW_A);

	int flag;
	MPI_Initialized(&flag); //дививмся чи ініціалізувалось
	if (flag)
		printf("\nMPI is initialized (proc 0) - TRUE");
	else
		printf("\nMPI is NOT initialized (proc 0) - False");

	printf("\nProcess 0 started work");

	printf("\nProcess 0 read memory");
	obj.read_memoryA(PROC0, subA);
	obj.read_memoryB(PROC0, subB);
	printf("\nProcess 0 END read memory");

	int sendRec = 1;
	MPI_Send(&sendRec, 1, MPI_INT, 1, 777, MPI_COMM_WORLD);
	printf("\nProcessing data (proc 0) ...");

	obj.mul(subA, subB, res, ROW_A, COLUMN_B);
	obj.tmp_storage(res, PROC0, 0, false);

	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Send(subB[i], COLUMN_B, MPI_INT, 5, 0, MPI_COMM_WORLD);

	int processes[6] = { 7, 6, 4, 3, 2, 1 };
	int sub_n3;
	for (int j = 0; j < 6; ++j) {
		sub_n3 = processes[j] < 5 ? COLUMN_B : COLUMN_B2;

		for (int i = 0; i < COLUMNMATRIX; ++i)
			MPI_Recv(subB[i], sub_n3, MPI_INT, 7, processes[j], MPI_COMM_WORLD, &Status);

		obj.mul(subA, subB, res, ROW_A, sub_n3);
		obj.tmp_storage(res, PROC0, processes[j], false);

		for (int i = 0; i < COLUMNMATRIX; ++i)
			MPI_Send(subB[i], sub_n3, MPI_INT, 5, processes[j], MPI_COMM_WORLD);
	}
	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Recv(subB[i], COLUMN_B2, MPI_INT, 7, 5, MPI_COMM_WORLD, &Status);

	obj.mul(subA, subB, res, ROW_A, COLUMN_B2);

	printf("\nProcess 0 started write memory!");
	obj.tmp_storage(res, PROC0, 5, true);
	printf("\nProcess 0 ended write memory!");

	MPI_Send(&sendRec, 1, MPI_INT, 1, 777, MPI_COMM_WORLD);
	printf("\nProcess 0 END work!!!");

	MPI_Send(&clockBegin, 1, MPI_DOUBLE, 7, 777, MPI_COMM_WORLD);
}

void CClaster::process_1() {
	CClaster obj(PROC1);

	int** subA;
	int** subB;
	int** res;

	ALLOCATE_MEMORY_A(ROW_A);
	ALLOCATE_MEMORY_B();
	ALLOCATE_MEMORY_R(ROW_A);

	int flag;
	MPI_Initialized(&flag); //дививмся чи ініціалізувалось
	if (flag)
		printf("\nMPI is initialized (proc 1) - TRUE");
	else
		printf("\nMPI is NOT initialized (proc 1) - False");

	printf("\nProcess 1 started work");

	int sendRec;
	MPI_Recv(&sendRec, 1, MPI_INT, 0, 777, MPI_COMM_WORLD, &Status); //очікування , на дозвіл читання з пам'яті.
	if (sendRec == 1)
		printf("\nProcess 0 give access to memory");

	printf("\nProcess 1 read memory");
	obj.read_memoryA(PROC1, subA);
	obj.read_memoryB(PROC1, subB);
	printf("\nProcess 1 END read memory");

	MPI_Send(&sendRec, 1, MPI_INT, 2, 777, MPI_COMM_WORLD);
	printf("\nProcessing data (proc 1) ...");

	obj.mul(subA, subB, res, ROW_A, COLUMN_B);
	obj.tmp_storage(res, PROC1, 1, false);

	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Send(subB[i], COLUMN_B, MPI_INT, 2, 1, MPI_COMM_WORLD);

	int processes[6] = { 5, 0, 7, 6, 4, 3 };
	int sub_n3;
	for (int j = 0; j < 6; ++j) {
		sub_n3 = processes[j] < 5 ? COLUMN_B : COLUMN_B2;
		for (int i = 0; i < COLUMNMATRIX; ++i)
			MPI_Recv(subB[i], sub_n3, MPI_INT, 5, processes[j], MPI_COMM_WORLD, &Status);

		obj.mul(subA, subB, res, ROW_A, sub_n3);
		obj.tmp_storage(res, PROC1, processes[j], false);

		for (int i = 0; i < COLUMNMATRIX; ++i)
			MPI_Send(subB[i], sub_n3, MPI_INT, 2, processes[j], MPI_COMM_WORLD);
	}
	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Recv(subB[i], COLUMN_B, MPI_INT, 5, 2, MPI_COMM_WORLD, &Status);

	obj.mul(subA, subB, res, ROW_A, COLUMN_B);

	MPI_Recv(&sendRec, 1, MPI_INT, 0, 777, MPI_COMM_WORLD, &Status); //очікування , на дозвіл читання з пам'яті.
	if (sendRec == 1)
		printf("\nProcess 0 give access to memory");

	printf("\nProcess 1 started write memory!");
	obj.tmp_storage(res, PROC1, 2, true);
	printf("\nProcess 1 ended write memory!");

	MPI_Send(&sendRec, 1, MPI_INT, 2, 777, MPI_COMM_WORLD);
	printf("\nProcess 1 END work!!!");
}

void CClaster::process_2() {
	CClaster obj(PROC2);

	int** subA;
	int** subB, **temp;
	int** res;

	ALLOCATE_MEMORY_A(ROW_A);
	ALLOCATE_MEMORY_B();
	ALLOCATE_MEMORY_temp();
	ALLOCATE_MEMORY_R(ROW_A);

	int flag;
	MPI_Initialized(&flag); //дививмся чи ініціалізувалось
	if (flag)
		printf("\nMPI is initialized (proc 2) - TRUE");
	else
		printf("\nMPI is NOT initialized (proc 2) - False");

	printf("\nProcess 2 started work");

	int sendRec;
	MPI_Recv(&sendRec, 1, MPI_INT, 1, 777, MPI_COMM_WORLD, &Status); //очікування , на дозвіл читання з пам'яті.
	if (sendRec == 1)
		printf("\nProcess 1 give access to memory");

	printf("\nProcess 2 read memory");
	obj.read_memoryA(PROC2, subA);
	obj.read_memoryB(PROC2, subB);
	printf("\nProcess 2 END read memory");

	MPI_Send(&sendRec, 1, MPI_INT, 3, 777, MPI_COMM_WORLD);
	printf("\nProcessing data (proc 2) ...");

	obj.mul(subA, subB, res, ROW_A, COLUMN_B);
	obj.tmp_storage(res, PROC2, 2, false);

	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Recv(temp[i], COLUMN_B, MPI_INT, 1, 1, MPI_COMM_WORLD, &Status);
	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Send(subB[i], COLUMN_B, MPI_INT, 3, 2, MPI_COMM_WORLD);

	obj.mul(subA, temp, res, ROW_A, COLUMN_B);
	obj.tmp_storage(res, PROC2, 1, false);

	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Recv(subB[i], COLUMN_B2, MPI_INT, 1, 5, MPI_COMM_WORLD, &Status);
	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Send(temp[i], COLUMN_B, MPI_INT, 3, 1, MPI_COMM_WORLD);

	obj.mul(subA, subB, res, ROW_A, COLUMN_B2);
	obj.tmp_storage(res, PROC2, 5, false);

	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Recv(temp[i], COLUMN_B, MPI_INT, 1, 0, MPI_COMM_WORLD, &Status);
	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Send(subB[i], COLUMN_B2, MPI_INT, 3, 5, MPI_COMM_WORLD);

	obj.mul(subA, temp, res, ROW_A, COLUMN_B);
	obj.tmp_storage(res, PROC2, 0, false);

	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Recv(subB[i], COLUMN_B2, MPI_INT, 1, 7, MPI_COMM_WORLD, &Status);
	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Send(temp[i], COLUMN_B, MPI_INT, 3, 0, MPI_COMM_WORLD);

	obj.mul(subA, subB, res, ROW_A, COLUMN_B2);
	obj.tmp_storage(res, PROC2, 7, false);

	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Recv(temp[i], COLUMN_B2, MPI_INT, 1, 6, MPI_COMM_WORLD, &Status);
	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Send(subB[i], COLUMN_B2, MPI_INT, 3, 7, MPI_COMM_WORLD);

	obj.mul(subA, temp, res, ROW_A, COLUMN_B2);
	obj.tmp_storage(res, PROC2, 6, false);

	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Recv(subB[i], COLUMN_B, MPI_INT, 1, 4, MPI_COMM_WORLD, &Status);
	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Send(temp[i], COLUMN_B2, MPI_INT, 3, 6, MPI_COMM_WORLD);

	obj.mul(subA, subB, res, ROW_A, COLUMN_B);
	obj.tmp_storage(res, PROC2, 4, false);

	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Recv(temp[i], COLUMN_B, MPI_INT, 1, 3, MPI_COMM_WORLD, &Status);
	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Send(subB[i], COLUMN_B, MPI_INT, 3, 4, MPI_COMM_WORLD);

	obj.mul(subA, temp, res, ROW_A, COLUMN_B);

	MPI_Recv(&sendRec, 1, MPI_INT, 1, 777, MPI_COMM_WORLD, &Status); //очікування , на дозвіл читання з пам'яті.
	if (sendRec == 1)
		printf("\nProcess 1 give access to memory");

	printf("\nProcess 2 started write memory!");
	obj.tmp_storage(res, PROC2, 3, true);
	printf("\nProcess 2 ended write memory!");

	MPI_Send(&sendRec, 1, MPI_INT, 3, 777, MPI_COMM_WORLD);
	printf("\nProcess 2 END work!!!");
}

void CClaster::process_3() {
	CClaster obj(PROC3);

	int** subA;
	int** subB;
	int** res;

	ALLOCATE_MEMORY_A(ROW_A);
	ALLOCATE_MEMORY_B();
	ALLOCATE_MEMORY_R(ROW_A);

	int flag;
	MPI_Initialized(&flag); //дививмся чи ініціалізувалось
	if (flag)
		printf("\nMPI is initialized (proc 3) - TRUE");
	else
		printf("\nMPI is NOT initialized (proc 3) - False");

	printf("\nProcess 3 started work");

	int sendRec;
	MPI_Recv(&sendRec, 1, MPI_INT, 2, 777, MPI_COMM_WORLD, &Status); //очікування , на дозвіл читання з пам'яті.
	if (sendRec == 1)
		printf("\nProcess 2 give access to memory");

	printf("\nProcess 3 read memory");
	obj.read_memoryA(PROC3, subA);
	obj.read_memoryB(PROC3, subB);
	printf("\nProcess 3 END read memory");

	MPI_Send(&sendRec, 1, MPI_INT, 4, 777, MPI_COMM_WORLD);
	printf("\nProcessing data (proc 3) ...");

	obj.mul(subA, subB, res, ROW_A, COLUMN_B);
	obj.tmp_storage(res, PROC3, 3, false);

	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Send(subB[i], COLUMN_B, MPI_INT, 4, 3, MPI_COMM_WORLD);

	int processes[6] = { 2, 1, 5, 0, 7, 6 };
	int sub_n3;
	for (int j = 0; j < 6; ++j) {
		sub_n3 = processes[j] < 5 ? COLUMN_B : COLUMN_B2;
		for (int i = 0; i < COLUMNMATRIX; ++i)
			MPI_Recv(subB[i], sub_n3, MPI_INT, 2, processes[j], MPI_COMM_WORLD, &Status);

		obj.mul(subA, subB, res, ROW_A, sub_n3);
		obj.tmp_storage(res, PROC3, processes[j], false);

		for (int i = 0; i < COLUMNMATRIX; ++i)
			MPI_Send(subB[i], sub_n3, MPI_INT, 4, processes[j], MPI_COMM_WORLD);
	}
	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Recv(subB[i], COLUMN_B, MPI_INT, 2, 4, MPI_COMM_WORLD, &Status);
	obj.mul(subA, subB, res, ROW_A, COLUMN_B);

	MPI_Recv(&sendRec, 1, MPI_INT, 2, 777, MPI_COMM_WORLD, &Status); //очікування , на дозвіл читання з пам'яті.
	if (sendRec == 1)
		printf("\nProcess 2 give access to memory");

	printf("\nProcess 3 started write memory!");
	obj.tmp_storage(res, PROC3, 4, true);
	printf("\nProcess 3 ended write memory!");

	MPI_Send(&sendRec, 1, MPI_INT, 4, 777, MPI_COMM_WORLD);
	printf("\nProcess 3 END work!!!");
}

void CClaster::process_4() {
	CClaster obj(PROC4);

	int** subA;
	int** subB, **temp;
	int** res;

	ALLOCATE_MEMORY_A(ROW_A);
	ALLOCATE_MEMORY_B();
	ALLOCATE_MEMORY_temp();
	ALLOCATE_MEMORY_R(ROW_A);

	int flag;
	MPI_Initialized(&flag); //дививмся чи ініціалізувалось
	if (flag)
		printf("\nMPI is initialized (proc 4) - TRUE");
	else
		printf("\nMPI is NOT initialized (proc 4) - False");

	printf("\nProcess 4 started work");

	int sendRec;
	MPI_Recv(&sendRec, 1, MPI_INT, 3, 777, MPI_COMM_WORLD, &Status); //очікування , на дозвіл читання з пам'яті.
	if (sendRec == 1)
		printf("\nProcess 3 give access to memory");

	printf("\nProcess 4 read memory");
	obj.read_memoryA(PROC4, subA);
	obj.read_memoryB(PROC4, subB);
	printf("\nProcess 4 END read memory");

	MPI_Send(&sendRec, 1, MPI_INT, 5, 777, MPI_COMM_WORLD);
	printf("\nProcessing data (proc 4) ...");

	obj.mul(subA, subB, res, ROW_A, COLUMN_B);
	obj.tmp_storage(res, PROC4, 4, false);

	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Recv(temp[i], COLUMN_B, MPI_INT, 3, 3, MPI_COMM_WORLD, &Status);
	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Send(subB[i], COLUMN_B, MPI_INT, 6, 4, MPI_COMM_WORLD);

	obj.mul(subA, temp, res, ROW_A, COLUMN_B);
	obj.tmp_storage(res, PROC4, 3, false);

	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Recv(subB[i], COLUMN_B, MPI_INT, 3, 2, MPI_COMM_WORLD, &Status);
	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Send(temp[i], COLUMN_B, MPI_INT, 6, 3, MPI_COMM_WORLD);

	obj.mul(subA, subB, res, ROW_A, COLUMN_B);
	obj.tmp_storage(res, PROC4, 2, false);

	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Recv(temp[i], COLUMN_B, MPI_INT, 3, 1, MPI_COMM_WORLD, &Status);
	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Send(subB[i], COLUMN_B, MPI_INT, 6, 2, MPI_COMM_WORLD);

	obj.mul(subA, temp, res, ROW_A, COLUMN_B);
	obj.tmp_storage(res, PROC4, 1, false);

	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Recv(subB[i], COLUMN_B2, MPI_INT, 3, 5, MPI_COMM_WORLD, &Status);
	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Send(temp[i], COLUMN_B, MPI_INT, 6, 1, MPI_COMM_WORLD);

	obj.mul(subA, subB, res, ROW_A, COLUMN_B2);
	obj.tmp_storage(res, PROC4, 5, false);

	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Recv(temp[i], COLUMN_B, MPI_INT, 3, 0, MPI_COMM_WORLD, &Status);
	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Send(subB[i], COLUMN_B2, MPI_INT, 6, 5, MPI_COMM_WORLD);

	obj.mul(subA, temp, res, ROW_A, COLUMN_B);
	obj.tmp_storage(res, PROC4, 0, false);

	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Recv(subB[i], COLUMN_B2, MPI_INT, 3, 7, MPI_COMM_WORLD, &Status);
	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Send(temp[i], COLUMN_B, MPI_INT, 6, 0, MPI_COMM_WORLD);

	obj.mul(subA, subB, res, ROW_A, COLUMN_B2);
	obj.tmp_storage(res, PROC4, 7, false);

	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Recv(temp[i], COLUMN_B2, MPI_INT, 3, 6, MPI_COMM_WORLD, &Status);
	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Send(subB[i], COLUMN_B2, MPI_INT, 6, 7, MPI_COMM_WORLD);

	obj.mul(subA, temp, res, ROW_A, COLUMN_B2);

	MPI_Recv(&sendRec, 1, MPI_INT, 3, 777, MPI_COMM_WORLD, &Status); //очікування , на дозвіл читання з пам'яті.
	if (sendRec == 1)
		printf("\nProcess 3 give access to memory");

	printf("\nProcess 4 started write memory!");
	obj.tmp_storage(res, PROC4, 6, true);
	printf("\nProcess 4 ended write memory!");

	MPI_Send(&sendRec, 1, MPI_INT, 5, 777, MPI_COMM_WORLD);
	printf("\nProcess 4 END work!!!");
}

void CClaster::process_5() {
	CClaster obj(PROC5);

	int** subA;
	int** subB, **temp;
	int** res;

	ALLOCATE_MEMORY_A(ROW_A);
	ALLOCATE_MEMORY_B();
	ALLOCATE_MEMORY_temp();
	ALLOCATE_MEMORY_R(ROW_A);

	int flag;
	MPI_Initialized(&flag); //дививмся чи ініціалізувалось
	if (flag)
		printf("\nMPI is initialized (proc 5) - TRUE");
	else
		printf("\nMPI is NOT initialized (proc 5) - False");

	printf("\nProcess 5 started work");

	int sendRec;
	MPI_Recv(&sendRec, 1, MPI_INT, 4, 777, MPI_COMM_WORLD, &Status); //очікування , на дозвіл читання з пам'яті.
	if (sendRec == 1)
		printf("\nProcess 4 give access to memory");

	printf("\nProcess 5 read memory");
	obj.read_memoryA(PROC5, subA);
	obj.read_memoryB(PROC5, subB);
	printf("\nProcess 5 END read memory");

	MPI_Send(&sendRec, 1, MPI_INT, 6, 777, MPI_COMM_WORLD);
	printf("\nProcessing data (proc 5) ...");

	obj.mul(subA, subB, res, ROW_A, COLUMN_B2);
	obj.tmp_storage(res, PROC5, 5, false);

	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Recv(temp[i], COLUMN_B, MPI_INT, 0, 0, MPI_COMM_WORLD, &Status);
	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Send(subB[i], COLUMN_B2, MPI_INT, 1, 5, MPI_COMM_WORLD);

	obj.mul(subA, temp, res, ROW_A, COLUMN_B);
	obj.tmp_storage(res, PROC5, 0, false);

	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Recv(subB[i], COLUMN_B2, MPI_INT, 0, 7, MPI_COMM_WORLD, &Status);
	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Send(temp[i], COLUMN_B, MPI_INT, 1, 0, MPI_COMM_WORLD);

	obj.mul(subA, subB, res, ROW_A, COLUMN_B2);
	obj.tmp_storage(res, PROC5, 7, false);

	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Recv(temp[i], COLUMN_B2, MPI_INT, 0, 6, MPI_COMM_WORLD, &Status);
	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Send(subB[i], COLUMN_B2, MPI_INT, 1, 7, MPI_COMM_WORLD);

	obj.mul(subA, temp, res, ROW_A, COLUMN_B);
	obj.tmp_storage(res, PROC5, 6, false);

	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Recv(subB[i], COLUMN_B, MPI_INT, 0, 4, MPI_COMM_WORLD, &Status);
	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Send(temp[i], COLUMN_B2, MPI_INT, 1, 6, MPI_COMM_WORLD);

	obj.mul(subA, subB, res, ROW_A, COLUMN_B);
	obj.tmp_storage(res, PROC5, 4, false);

	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Recv(temp[i], COLUMN_B, MPI_INT, 0, 3, MPI_COMM_WORLD, &Status);
	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Send(subB[i], COLUMN_B, MPI_INT, 1, 4, MPI_COMM_WORLD);

	obj.mul(subA, temp, res, ROW_A, COLUMN_B);
	obj.tmp_storage(res, PROC5, 3, false);

	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Recv(subB[i], COLUMN_B, MPI_INT, 0, 2, MPI_COMM_WORLD, &Status);
	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Send(temp[i], COLUMN_B, MPI_INT, 1, 3, MPI_COMM_WORLD);

	obj.mul(subA, subB, res, ROW_A, COLUMN_B);
	obj.tmp_storage(res, PROC5, 2, false);

	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Recv(temp[i], COLUMN_B, MPI_INT, 0, 1, MPI_COMM_WORLD, &Status);
	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Send(subB[i], COLUMN_B, MPI_INT, 1, 2, MPI_COMM_WORLD);

	obj.mul(subA, temp, res, ROW_A, COLUMN_B);

	MPI_Recv(&sendRec, 1, MPI_INT, 4, 777, MPI_COMM_WORLD, &Status); //очікування , на дозвіл читання з пам'яті.
	if (sendRec == 1)
		printf("\nProcess 4 give access to memory");

	printf("\nProcess 5 started write memory!");
	obj.tmp_storage(res, PROC5, 1, true);
	printf("\nProcess 5 ended write memory!");

	MPI_Send(&sendRec, 1, MPI_INT, 6, 777, MPI_COMM_WORLD);
	printf("\nProcess 5 END work!!!");
}

void CClaster::process_6() {
	CClaster obj(PROC6);

	int** subA;
	int** subB;
	int** res;

	ALLOCATE_MEMORY_A(ROW_A2);
	ALLOCATE_MEMORY_B();
	ALLOCATE_MEMORY_R(ROW_A2);

	int flag;
	MPI_Initialized(&flag); //дививмся чи ініціалізувалось
	if (flag)
		printf("\nMPI is initialized (proc 6) - TRUE");
	else
		printf("\nMPI is NOT initialized (proc 6) - False");

	printf("\nProcess 6 started work");

	int sendRec;
	MPI_Recv(&sendRec, 1, MPI_INT, 5, 777, MPI_COMM_WORLD, &Status); //очікування , на дозвіл читання з пам'яті.
	if (sendRec == 1)
		printf("\nProcess 5 give access to memory");

	printf("\nProcess 6 read memory");
	obj.read_memoryA(PROC6, subA);
	obj.read_memoryB(PROC6, subB);
	printf("\nProcess 6 END read memory");

	MPI_Send(&sendRec, 1, MPI_INT, 7, 777, MPI_COMM_WORLD);
	printf("\nProcessing data (proc 6) ...");

	obj.mul(subA, subB, res, ROW_A2, COLUMN_B2);
	obj.tmp_storage(res, PROC6, 6, false);

	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Send(subB[i], COLUMN_B2, MPI_INT, 7, 6, MPI_COMM_WORLD);

	int processes[6] = { 4, 3, 2, 1, 5, 0 };
	int sub_n3;
	for (int j = 0; j < 6; ++j) {
		sub_n3 = processes[j] < 5 ? COLUMN_B : COLUMN_B2;
		for (int i = 0; i < COLUMNMATRIX; ++i)
			MPI_Recv(subB[i], sub_n3, MPI_INT, 4, processes[j], MPI_COMM_WORLD, &Status);

		obj.mul(subA, subB, res, ROW_A2, sub_n3);
		obj.tmp_storage(res, PROC6, processes[j], false);

		for (int i = 0; i < COLUMNMATRIX; ++i)
			MPI_Send(subB[i], sub_n3, MPI_INT, 7, processes[j], MPI_COMM_WORLD);
	}
	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Recv(subB[i], COLUMN_B2, MPI_INT, 4, 7, MPI_COMM_WORLD, &Status);

	obj.mul(subA, subB, res, ROW_A2, COLUMN_B2);

	MPI_Recv(&sendRec, 1, MPI_INT, 5, 777, MPI_COMM_WORLD, &Status); //очікування , на дозвіл читання з пам'яті.
	if (sendRec == 1)
		printf("\nProcess 5 give access to memory");

	printf("\nProcess 6 started write memory!");
	obj.tmp_storage(res, PROC6, 7, true);
	printf("\nProcess 6 ended write memory!");

	MPI_Send(&sendRec, 1, MPI_INT, 7, 777, MPI_COMM_WORLD);
	printf("\nProcess 6 END work!!!");
}

void CClaster::process_7() {
	Timer timer;
	double clockBegin;

	CClaster obj(PROC7);

	int** subA;
	int** subB, ** temp;
	int** res;

	ALLOCATE_MEMORY_A(ROW_A2);
	ALLOCATE_MEMORY_B();
	ALLOCATE_MEMORY_temp();
	ALLOCATE_MEMORY_R(ROW_A2);

	int flag;
	MPI_Initialized(&flag); //дививмся чи ініціалізувалось
	if (flag)
		printf("\nMPI is initialized (proc 7) - TRUE");
	else
		printf("\nMPI is NOT initialized (proc 7) - False");

	printf("\nProcess 7 started work");

	int sendRec;
	MPI_Recv(&sendRec, 1, MPI_INT, 6, 777, MPI_COMM_WORLD, &Status); //очікування , на дозвіл читання з пам'яті.
	if (sendRec == 1)
		printf("\nProcess 6 give access to memory");

	printf("\nProcess 7 read memory");
	obj.read_memoryA(PROC7, subA);
	obj.read_memoryB(PROC7, subB);
	printf("\nProcess 7 END read memory");

	printf("\nProcessing data (proc 7) ...");

	obj.mul(subA, subB, res, ROW_A2, COLUMN_B2);
	obj.tmp_storage(res, PROC7, 7, false);

	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Recv(temp[i], COLUMN_B2, MPI_INT, 6, 6, MPI_COMM_WORLD, &Status);
	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Send(subB[i], COLUMN_B2, MPI_INT, 0, 7, MPI_COMM_WORLD);

	obj.mul(subA, temp, res, ROW_A2, COLUMN_B2);
	obj.tmp_storage(res, PROC7, 6, false);

	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Recv(subB[i], COLUMN_B, MPI_INT, 6, 4, MPI_COMM_WORLD, &Status);
	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Send(temp[i], COLUMN_B2, MPI_INT, 0, 6, MPI_COMM_WORLD);

	obj.mul(subA, subB, res, ROW_A2, COLUMN_B);
	obj.tmp_storage(res, PROC7, 4, false);

	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Recv(temp[i], COLUMN_B, MPI_INT, 6, 3, MPI_COMM_WORLD, &Status);
	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Send(subB[i], COLUMN_B, MPI_INT, 0, 4, MPI_COMM_WORLD);

	obj.mul(subA, temp, res, ROW_A2, COLUMN_B);
	obj.tmp_storage(res, PROC7, 3, false);

	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Recv(subB[i], COLUMN_B, MPI_INT, 6, 2, MPI_COMM_WORLD, &Status);
	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Send(temp[i], COLUMN_B, MPI_INT, 0, 3, MPI_COMM_WORLD);

	obj.mul(subA, subB, res, ROW_A2, COLUMN_B);
	obj.tmp_storage(res, PROC7, 2, false);

	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Recv(temp[i], COLUMN_B, MPI_INT, 6, 1, MPI_COMM_WORLD, &Status);
	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Send(subB[i], COLUMN_B, MPI_INT, 0, 2, MPI_COMM_WORLD);

	obj.mul(subA, temp, res, ROW_A2, COLUMN_B);
	obj.tmp_storage(res, PROC7, 1, false);

	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Recv(subB[i], COLUMN_B2, MPI_INT, 6, 5, MPI_COMM_WORLD, &Status);
	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Send(temp[i], COLUMN_B, MPI_INT, 0, 1, MPI_COMM_WORLD);

	obj.mul(subA, subB, res, ROW_A2, COLUMN_B2);
	obj.tmp_storage(res, PROC7, 5, false);

	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Recv(temp[i], COLUMN_B, MPI_INT, 6, 0, MPI_COMM_WORLD, &Status);
	for (int i = 0; i < COLUMNMATRIX; ++i)
		MPI_Send(subB[i], COLUMN_B2, MPI_INT, 0, 5, MPI_COMM_WORLD);

	obj.mul(subA, temp, res, ROW_A2, COLUMN_B);

	MPI_Recv(&clockBegin, 1, MPI_DOUBLE, 0, 777, MPI_COMM_WORLD, &Status); //receiving start time

	MPI_Recv(&sendRec, 1, MPI_INT, 6, 777, MPI_COMM_WORLD, &Status); //очікування , на дозвіл читання з пам'яті.
	if (sendRec == 1)
		printf("\nProcess 6 give access to memory");

	printf("\nProcess 7 started write memory!");
	obj.tmp_storage(res, PROC7, 0, true);
	printf("\nProcess 7 ended write memory!");

	printf("\nProcess 7 END work!!!");



	double clockEnd = timer.getCurrentTime();

	printf("\n\nThe program with MPI run for %lf seconds...\n\n", clockEnd - clockBegin);
	FILE* f1;
	f1 = fopen("resTime.txt", "a");
	fprintf(f1, "The program with MPI run for %lf seconds...\n\n", clockEnd - clockBegin);
	fclose(f1);
}