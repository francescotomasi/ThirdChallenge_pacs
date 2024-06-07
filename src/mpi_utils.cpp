#include "mpi_utils.hpp"
#include "muparser_fun.hpp"
#include <mpi.h>
#include <vector>
#include <string>

void broadcast_parameters(int rank, int& max_it, double& tol_res, std::string& forceFunction, int& gridDim) {
    MPI_Bcast(&max_it, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tol_res, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&gridDim, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&forceFunction[0], forceFunction.size() + 1, MPI_CHAR, 0, MPI_COMM_WORLD);
}

void distribute_grid(int rank, int size, int gridDim, MuparserFun& fun, double h, 
                     std::vector<double>& f_local, int& local_nrows, 
                     std::vector<int>& counts, std::vector<int>& displs) {
    std::vector<double> f_tot(gridDim * gridDim, 0.0);

    local_nrows = (gridDim % size > rank) ? gridDim / size + 1 : gridDim / size;
    f_local.resize(gridDim * local_nrows);

    if (rank == 0) {
        for (int i = 1; i < gridDim - 1; ++i) {
            for (int j = 1; j < gridDim - 1; ++j) {
                double x = i * h;
                double y = j * h;
                f_tot[i * gridDim + j] = fun(x, y);
            }
        }
        
        counts.resize(size);
        displs.resize(size);
        for (int i = 0; i < size; i++) {
            int local_r = (gridDim % size > i) ? gridDim / size + 1 : gridDim / size;
            counts[i] = gridDim * local_r;
            if (i != 0) {
                displs[i] = displs[i - 1] + counts[i - 1];
            }
        }
    }

    MPI_Bcast(counts.data(), size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(displs.data(), size, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Scatterv(f_tot.data(), counts.data(), displs.data(), MPI_DOUBLE, 
                 f_local.data(), gridDim * local_nrows, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void gather_results(int rank, int size, int gridDim, int local_nrows, 
                    std::vector<double>& U_local, std::vector<double>& U_real_local, 
                    std::vector<double>& U_real_tot, std::vector<int>& counts, 
                    std::vector<int>& displs) {
    for (int i = 1; i < local_nrows + 1; i++) {
        for (int j = 1; j < gridDim - 1; j++) {
            U_real_local[(i - 1) * gridDim + j] = U_local[i * gridDim + j];
        }
    }

    MPI_Gatherv(U_real_local.data(), gridDim * local_nrows, MPI_DOUBLE, 
                U_real_tot.data(), counts.data(), displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
}
