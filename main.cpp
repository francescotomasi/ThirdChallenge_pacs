#include "mpi_utils.hpp"
#include "grid_utils.hpp"
#include "config.hpp"
#include "algorithms.hpp"
#include "muparser_fun.hpp"
#include <vector>
#include <string>
#include <iostream>
#include <mpi.h>

int main(int argc, char** argv) {
    int rank, size, provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int max_it;
    double tol_res;
    std::string forceFunction;
    int gridDim;

    if (rank == 0) {
        read_config("data.json", max_it, tol_res, forceFunction, gridDim);
    }

    broadcast_parameters(rank, max_it, tol_res, forceFunction, gridDim);

    double h = 1.0 / (gridDim - 1);
    MuparserFun fun(forceFunction, 2);

    int local_nrows;
    std::vector<double> f_local;
    std::vector<int> counts, displs;
    distribute_grid(rank, size, gridDim, fun, h, f_local, local_nrows, counts, displs);

    if (rank == 0) {
        std::cout << "Force function da rank 0:\n";
        print_grid(f_local, local_nrows, gridDim);
    }

    std::vector<double> U_local(gridDim * (local_nrows + 2), 0.0);
    std::vector<double> U_local_new(gridDim * (local_nrows + 2), 0.0);

    // Choose algorithm
    bool use_parallel = true;  // Set this based on your requirements

    if (use_parallel) {
        run_parallel_algorithm(rank, size, gridDim, local_nrows, max_it, tol_res, 
                               f_local, U_local, U_local_new, counts, displs);
    } else {
        run_sequential_algorithm(gridDim, max_it, tol_res, f_local, U_local, U_local_new);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 1) {
        std::cout << "\nU_local da rank 1 dopo l'algoritmo:\n";
        print_grid(U_local, local_nrows + 2, gridDim);
    }

    std::vector<double> U_real_local(gridDim * local_nrows, 0.0);
    std::vector<double> U_real_tot(gridDim * gridDim, 0.0);

    gather_results(rank, size, gridDim, local_nrows, U_local, U_real_local, U_real_tot, counts, displs);

    if (rank == 0) {
        std::cout << "\nU_real_tot dopo l'algoritmo:\n";
        print_grid(U_real_tot, gridDim, gridDim);
    }

    MPI_Finalize();
    return 0;
}
