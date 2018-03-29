/*
==================================================================
Name : parallel_task2_2.c
Author : Nadezhda Kiryushkina
Version :
Copyright : No copyright notice
Description : Parallel task 2 in C, Ansistyle
==================================================================
*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#define INIT_SIZE 128
#define MAXSTEPS 40
#define EPS 0.01
#define PI 3.14159265
#define Max(a,b) ((a)>(b)?(a):(b))
double ** init_array(int m, int n);
void bounCond(int m, int n, double **u, int processId, int commSize);
double relax(int m, int n, double **u, double **unew, double wn, int processId);
void neighbors(int processId, int commSize, int *below, int *above);
int updateBoundCount(int m, int n, double **u, int processId, int below, int above );
void write_file(int m, int n, double **u, int commSize);

int main(int argc, char *argv[]) {
	int size = INIT_SIZE, segment, processId, commSize, below, above;
	long iter = 0;
	double err = EPS, localEps, currEps = 1.0, start, finish, time;
	double **u, **unew;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &processId);
	MPI_Comm_size(MPI_COMM_WORLD, &commSize);
	if(processId == 0) {
		start = MPI_Wtime();
	}
	segment = size/commSize;
	u = init_array(segment, size);
	unew = init_array(segment, size);
	bounCond(segment, size, u, processId, commSize);
	bounCond(segment, size, unew, processId, commSize);
	neighbors(processId, commSize, &below, &above);
	double radius = 1.0 - (PI * PI) / (4 * (size + 2) * (size + 2));
	radius *= radius;
	double wn = (1.0 - 0.5 * radius);
	while (currEps > err) {
		if(iter > MAXSTEPS) {
		printf("Iteration terminated %6d", MAXSTEPS);
			return 0;
		}
		localEps = relax(segment, size, u, unew, wn, processId);
		wn = (1.0 - 0.25 * radius * wn);
		MPI_Allreduce( &localEps, &currEps, 1, MPI_DOUBLE,MPI_MAX, MPI_COMM_WORLD );
		if(iter%2 == 0) {
			if(processId == 0) {
				printf("iter: %6d\tdel: %lf\tgdel: %lf\n",iter,localEps,currEps);
			}
		}
		updateBoundCount(segment, size, u, processId, below, above);
		iter++;
	}
	if (processId == 0) {
		finish=MPI_Wtime();
		time = finish - start;
		printf("Stopped at iteration %d\n",iter);
		printf("The maximum error = %f\n",currEps);
		printf("Time = %f\n",time);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	write_file(segment, size, u, commSize);
	free(u);
	free(unew);
	return 0;
}
int updateBoundCount(int segment, int size, double **u, int processId, int below, int above ) {
	MPI_Status status[6];
	MPI_Sendrecv( u[segment]+1, size, MPI_DOUBLE, below, 0, u[0]+1, size, MPI_DOUBLE, above, 0, MPI_COMM_WORLD, status );
	MPI_Sendrecv( u[1]+1, size, MPI_DOUBLE, above, 1, u[segment+1]+1, size, MPI_DOUBLE, below, 1, MPI_COMM_WORLD, status );
	return 0;
}
double** init_array(int m, int n) {
	int i, j;
	double **u;
	u = (double **) malloc((unsigned) (m+2)*sizeof(double*));
	for (i = 0; i <=m+1; i++) {
		u[i] = (double *) malloc((unsigned) (n+2)*sizeof(double));
	}
	for(i = 0; i <=m+1; i++) {
		for(j = 0; j <=n+1; j++) {
			u[i][j] = 0.0;
		}
	}
	return u;
}
void bounCond(int m, int n, double **u, int processId, int commSize) {
	int i;
	if (processId == commSize - 1) {
		for (i = 0; i <=n+1; i++) {
			u[m+1][i] = exp(-(PI*i/(n+2))/2);
		}
	}
	for (i = 0; i <=m+1; i++) {
		double rel = ((double) i + 1.0*processId*m)/ ((double)n+2);
		u[i][0] = sin(PI*rel/2);
		u[i][n+1] = sin(PI*rel/2) * exp(-PI/2);
	}
}
double relax( int m, int n, double **u, double **unew, double wn, int processId) {
	int i, j;
	double localEps = 0.0;
	for (i = 1; i <=m; i++) {
		for (j = 1; j <=n; j++) {
			double uij = ( u[i ][j+1] + u[i+1][j ] + u[i-1][j ] + u[i ][j-1] )*0.25;
			unew[i][j] = wn*uij + (1.0 - wn)*u[i][j];
			localEps = Max(localEps, fabs(unew[i][j] - u[i][j]));
		}
	}
	for (i = 1; i <=m; i++) {
		for (j = 1; j <=n; j++) {
			u[i][j] = unew[i][j];
		}
	}
	return localEps;
}
void neighbors(int processId, int commSize, int *below, int *above) {
	if(commSize == 1) {
		*above = MPI_PROC_NULL;
		*below = MPI_PROC_NULL;
		return;
	}
	if(processId == 0) {
		*below = processId+1;
		*above = MPI_PROC_NULL;
	} else if(processId == commSize - 1) {
		*below = MPI_PROC_NULL;
		*above = processId - 1;
	} else {
		*below = processId + 1;
		*above = processId - 1;
	}
}
void write_file(int m, int n, double **u, int commSize) {
	int i, j;
	char filename[50];
	MPI_File fd;
	MPI_Status status;
	int buf_size = (m * (15 * n + 2) + 2);
	char *buf = (char *) calloc(buf_size, sizeof(char));
	char s[14];
	sprintf(filename, "m=%d_proc=%d_eps=%lf.txt", n, commSize, EPS);
	MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fd);
	for (j = 0; j <=n+1; j++) {
		for (i = 0; i <=m+1; i++) {
			sprintf(s, "%6.6g ", u[i][j]);
			strcat(buf, s);
		}
		sprintf(s, "\n");
		strcat(buf, s);
	}
	MPI_File_write_ordered(fd, buf, (int) strlen(buf), MPI_CHAR, &status);
	MPI_File_close(&fd);
	free(buf);
	return;
}