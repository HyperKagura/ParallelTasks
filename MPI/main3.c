/*
==================================================================
Name : parallel_task3.cpp
Author : Nadezhda Kiryushkina
Version :
Copyright : No copyright notice
==================================================================
*/
#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#define M_SIZE 128
double bound_x_u (double y, double z)
{
return (1 - fabs(2*y - 1)) * (1 - fabs(2*z - 1));
}
double bound_x_l (double y, double z)
{
	return 0;
}
double bound_y_u (double x, double z)
{
	return 0;
}
double bound_y_l (double x, double z)
{
	return 0;
}
double bound_z_u (double x, double y)
{
	return 0;
}
double bound_z_l (double x, double y)
{
	return 0;
}

int main(int argc, char** argv)
{
	MPI_Init(&argc, &argv);
	int rank, size, size_dim;
	int rank_x, rank_y, rank_z;
	double time;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0)
	{
		std::cout << "Allocated " << size << " processors" << std::endl;
	}
	MPI_Comm comm;
	size_dim = size = pow(size+0.5, 1.0/3);
	size*=size*size;
	if (rank == 0)
	{
		std::cout << "Using " << size << " processors" << std::endl;
		time = MPI_
		Wtime();
	}
	int color;
	if (rank < size) color = 1;
	else color = 2;
	MPI_Comm_split(MPI_COMM_WORLD, color, rank, &comm);
	if (rank>=size) {
		MPI_Finalize();
		return 0;
	}
	MPI_Comm_rank(comm, &rank);
	int N;
	N = M_SIZE;
	int n = (N-1)/size_dim + 1;
	int n_x, n_y, n_z;
	n_x = n_y = n_z = n;
	rank_x = rank%size_dim;
	rank_y = (rank/size_dim)%size_dim;
	rank_z = (rank/size_dim)/size_dim;
	if (rank_x == size_dim - 1) n_x = N - n*(size_dim - 1);
	if (rank_y == size_dim - 1) n_y = N - n*(size_dim - 1);
	if (rank_z == size_dim - 1) n_z = N - n*(size_dim - 1);
	double* u = (double*)calloc((n_x+2)*(n_y+2)*(n_z+2), sizeof(double));
	double* new_u = (double*)calloc((n_x+2)*(n_y+2)*(n_z+2), sizeof(double));
	if (rank_x == size_dim - 1) {
		for (int j = 1; j<=n_y; j++) {
			for (int k = 1; k<= n_z; k++) {
				double y, z;
				y = double(rank_y*n + j)/(N+1);
				z = double(rank_z*n + k)/(N+1);
				u[n_x+1 + j*(n_x+2) + k*(n_x+2)*(n_y+2)] = new_u[n_x+1 + j*(n_x+2) +
				k*(n_x+2)*(n_y+2)] = bound_x_u(y, z);
			}
		}
	}
	if (rank_x == 0) {
		for (int j = 1; j<=n_y; j++) {
			for (int k = 1; k<= n_z; k++) {
				double y, z;
				y = double(rank_y*n + j)/(N+1);
				z = double(rank_z*n + k)/(N+1);
				u[j*(n_x+2) + k*(n_x+2)*(n_y+2)] = new_u[j*(n_x+2) + k*(n_x+2)*(n_y+2)] =
				bound_x_l(y, z);
			}
		}
	}
	if (rank_y == size_dim - 1) {
		for (int i = 1; i<=n_x; i++) {
			for (int k = 1; k<= n_z; k++) {
				double x, z;
				x = double(rank_x*n + i)/(N+1);
				z = double(rank_z*n + k)/(N+1);
				u[i + (n_y + 1)*(n_x+2) + k*(n_x+2)*(n_y+2)] = new_u[i + (n_y + 1)*(n_x+2) +
				k*(n_x+2)*(n_y+2)] = bound_y_u(x, z);
			}
		}
	}
	if (rank_y == 0) {
		for (int i = 1; i<=n_x; i++) {
			for (int k = 1; k<= n_z; k++) {
				double x, z;
				x = double(rank_x*n + i)/(N+1);
				z = double(rank_z*n + k)/(N+1);
				u[i + k*(n_x+2)*(n_y+2)] = new_u[i + k*(n_x+2)*(n_y+2)] = bound_y_l(x, z);
			}
		}
	}
	if (rank_z == size_dim 1) {
		for (int i = 1; i<=n_x; i++) {
			for (int j = 1; j<= n_y; j++) {
				double x, y;
				x = double(rank_x*n + i)/(N+1);
				y = double(rank_y*n + j)/(N+1);
				u[i + j*(n_x+2) + (n_z + 1)*(n_x+2)*(n_y+2)] = new_u[i + j*(n_x+2) + (n_z +
				1)*(n_x+2)*(n_y+2)] = bound_z_u(x, y);
			}
		}
	}
	if (rank_z == 0) {
		for (int i = 1; i<=n_x; i++) {
			for (int j = 1; j<= n_y; j++) {
				double x, y;
				x = double(rank_x*n + i)/(N+1);
				y = double(rank_y*n + j)/(N+1);
				u[i + j * (n_x + 2)] = new_u[i + j * (n_x + 2)] = bound_z_l(x, y);
			}
		}
	}
	int iter = 1;
	double rho = 0;
	//rho = 1.0 - 2 * M_PI / (N+1);
	/*rho = cos(M_PI/(N+1)) / ( 1 + sin(M_PI/(N+1)) );
	rho *= rho;*/
	if (rank == 0) std::cout << "Rho = " << rho << std::endl;
	double * buf_x = (double *) malloc(n_y * n_z * sizeof(double));
	double * buf_y = (double *) malloc(n_x * n_z * sizeof(double));
	double * buf_z = (double *) malloc(n_x * n_y * sizeof(double));
	double* err_array = (double*)calloc(2, sizeof(double));
	double w;
	MPI_Status status;
	do {
		if (iter == 0) w = 0;
		else if (iter == 1) w = 1.0l / (1.0l - rho*rho / 2);
		else w = 1.0l / (1.0l - rho*rho*w/4);
		#pragma omp parallel for
		for (int k = 1; k <= n_z; k++)
			for (int j = 1; j <= n_y; j++)
				for (int i = 1; i <= n_x; i++) {
					new_u[i + (n_x + 2) * (j + (n_y + 2) * k)] = (1.0l - w) * u[i + (n_x + 2) * (j + (n_y + 2) * k)]
					+ w * (
					u[i - 1 + (n_x + 2) * (j + (n_y + 2) * k)] + u[i + 1 + (n_x + 2) * (j + (n_y +
					2) * k)]
					+ u[i + (n_x + 2) * (j - 1 + (n_y + 2) * k)] + u[i + (n_x + 2) * (j + 1 + (n_y
					+ 2) * k)]
					+ u[i + (n_x + 2) * (j + (n_y + 2) * (k - 1))] + u[i + (n_x + 2) * (j + (n_y +
					2) * (k + 1))]
					) / 6.0l;
				}
		int dest, source;
		dest = rank - 1;
		source = rank + 1;
		if (rank_x == 0) dest += size_dim;
		if (rank_x == size_dim - 1) source = size_dim;
		for (int k = 1; k <= n_z; k++)
			for(int j = 1; j <= n_y; j ++)
				buf_x [j - 1 + (k - 1) * n_y] = new_u[1 + (n_x + 2) * (j + (n_y + 2) * k)];
		MPI_Sendrecv_replace(buf_x, n_y*n_z, MPI_DOUBLE, dest, 1, source, 1, comm, &status);
		if (rank_x < size_dim - 1) {
			#pragma omp parallel for
			for (int k = 1; k <= n_z; k++)
				for(int j = 1; j <= n_y; j ++)
					new_u[n_x + 1 + (n_x + 2) * (j + (n_y + 2) * k)] = buf_x [j - 1 + (k - 1) * n_y];
		}
		dest = rank + 1;
		source = rank - 1;
		if (rank_x == 0) source += size_dim;
		if (rank_x == size_dim - 1) dest = size_dim;
		for (int k = 1; k <= n_z; k++)
			for(int j = 1; j <= n_y; j ++)
				buf_x [j - 1 + (k - 1) * n_y] = new_u[n_x + (n_x + 2) * (j + (n_y + 2) * k)];
		MPI_Sendrecv_replace(buf_x, n_y*n_z, MPI_DOUBLE, dest, 2, source, 2, comm, &status);
		if (rank_x > 0) {
			#pragma omp parallel for
			for (int k = 1; k <= n_z; k++)
				for(int j = 1; j <= n_y; j ++)
					new_u[(n_x + 2) * (j + (n_y + 2) * k)] = buf_x [j - 1 + (k - 1) * n_y];
		}
		dest = rank - size_dim;
		source = rank + size_dim;
		if (rank_y == 0) dest += size_dim * size_dim;
		if (rank_y == size_dim - 1) source -= size_dim * size_dim;
		for (int k = 1; k <= n_z; k++)
			for(int i = 1; i <= n_x; i ++)
				buf_y [i - 1 + (k - 1) * n_x] = new_u[i + (n_x + 2) * (1 + (n_y + 2) * k)];
		MPI_Sendrecv_replace(buf_y, n_x * n_z, MPI_DOUBLE, dest, 3, source, 3, comm, &status);
		if (rank_y < size_dim 1) {
			for (int k = 1; k <= n_z; k++)
				for(int i = 1; i <= n_x; i ++)
					new_u[i + (n_x + 2) * (n_y + 1 + (n_y + 2) * k)] = buf_y [i - 1 + (k - 1) * n_x];
		}
		dest = rank + size_dim;
		source = rank - size_dim;
		if (rank_y == 0) source += size_dim * size_dim;
		if (rank_y == size_dim - 1) dest -= size_dim * size_dim;
		for (int k = 1; k <= n_z; k++)
			for(int i = 1; i <= n_x; i++)
				buf_y [i - 1 + (k - 1) * n_x] = new_u[i + (n_x + 2) * (n_y + (n_y + 2) * k)];
		MPI_Sendrecv_replace(buf_y, n_x * n_z, MPI_DOUBLE, dest, 4, source, 4, comm, &status);
		if (rank_y > 0) {
			for (int k = 1; k <= n_z; k++)
				for(int i = 1; i <= n_x; i ++)
					new_u[i + (n_x + 2) * (n_y + 2) * k] = buf_y [i - 1 + (k - 1) * n_x];
		}
		dest = rank - size_dim * size_dim;
		source = rank + size_dim * size_dim;
		if (rank_z == 0) dest += size_dim * size_dim * size_dim;
		if (rank_z == size_dim - 1) source -= size_dim * size_dim * size_dim;
		for (int j = 1; j <= n_y; j++)
			for(int i = 1; i <= n_x; i ++)
				buf_z [i - 1 + (j - 1) * n_x] = new_u[i + (n_x + 2) * (j + (n_y + 2))];
		MPI_Sendrecv_replace(buf_z, n_x * n_y, MPI_DOUBLE, dest, 5, source, 5, comm, &status);
		if (rank_z < size_dim - 1) {
			for (int j = 1; j <= n_y; j++)
				for(int i = 1; i <= n_x; i ++)
					new_u[i + (n_x + 2) * (j + (n_y + 2) * (n_z + 1))] = buf_z [i - 1 + (j - 1) * n_x];
		}
		dest = rank + size_dim * size_dim;
		source = rank - size_dim * size_dim;
		if (rank_z == 0) source += size_dim * size_dim * size_dim;
		if (rank_z == size_dim - 1) dest -= size_dim * size_dim * size_dim;
		for (int j = 1; j <= n_y; j++)
			for(int i = 1; i <= n_x; i++)
				buf_z [i - 1 + (j - 1) * n_x] = new_u[i + (n_x + 2) * (j + (n_y + 2) * n_z)];
		MPI_Sendrecv_replace(buf_z, n_x * n_y, MPI_DOUBLE, dest, 6, source, 6, comm, &status);
		if (rank_z > 0) {
			for (int j = 1; j <= n_y; j++)
				for(int i = 1; i <= n_x; i ++)
					new_u[i + (n_x + 2) * j] = buf_z [i - 1 + (j - 1) * n_x];
		}
		double * tmp = u;
		u = new_u;
		new_u = tmp;
		err_array[0] = 0;
		err_array[1] = 0;
		for (int i=0; i < n_x + 2; i++) {
			for (int j = 0; j < n_y + 2; j++) {
				for (int k = 0; k < n_z + 2; k++) {
					if ( fabs(new_u[i + (n_x + 2) * (j + (n_y + 2) * k)] u[i + (n_x + 2) * (j + (n_y +
					2) * k)]) > err_array[0])
						err_array[0] = fabs(new_u[i + (n_x + 2) * (j + (n_y + 2) * k)] u[i + (n_x +
					2) * (j + (n_y + 2) * k)]);
					if ( fabs(u[i + (n_x + 2) * (j + (n_y + 2) * k)]) > err_array[1])
						err_array[1] = fabs(u[i + (n_x + 2) * (j + (n_y + 2) * k)]);
				}
			}
		}
		MPI_Allreduce(MPI_IN_PLACE, err_array, 2, MPI_DOUBLE, MPI_MAX, comm);
		iter++;
	}
	while ((err_array[0] / err_array[1] > 0.0001) && (iter <= 1000));
	if (rank == 0) {
		std::cout << "Iterations = " << iter << ", error norm = "<< err_array[0] / err_array[1] << std::endl;
		time += MPI_Wtime();
		std::cout << "Time is: " << time << std::endl;
	}
	MPI_Finalize();
	return 0;
}