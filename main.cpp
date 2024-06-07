#include <mpi.h>
#include <omp.h>
#include <vector>
#include <nlohmann/json.hpp>
#include "muparser_fun.hpp"
#include <fstream>
#include "algorithms.hpp"

using json = nlohmann::json;

void print_grids( const std::vector<double>& grid, int row, int col) {
    std::cout << "\nGrid:\n";
    for ( int i = 0; i < row; ++i) {
        for ( int j = 0; j < col; ++j) {
            std::cout << grid[i * col + j] << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

void write_vtks(const std::string& filename, const std::vector<double>& grid, int row, int col) {
    std::ofstream file;
    file.open(filename);

    file << "# vtk DataFile Version 3.0\n";
    file << "Grid data\n";
    file << "ASCII\n";
    file << "DATASET RECTILINEAR_GRID\n";
    file << "DIMENSIONS " << col << " " << row << " 1\n";
    
    file << "X_COORDINATES " << col << " float\n";
    for (int i = 0; i < col; ++i) {
        file << i << " ";
    }
    file << "\n";
    
    file << "Y_COORDINATES " << row << " float\n";
    for (int i = 0; i < row; ++i) {
        file << i << " ";
    }
    file << "\n";
    
    file << "Z_COORDINATES 1 float\n";
    file << "0\n";
    
    file << "POINT_DATA " << row * col << "\n";
    file << "SCALARS solution float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            file << grid[i * col + j] << " ";
        }
        file << "\n";
    }
    
    file.close();
}

int main(int argc, char** argv){

    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    
    int rank,size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int max_it;
    double tol_res;
    std::string forceFunction;
    int gridDim ;

    if (rank == 0) {
        std::ifstream z("data.json");
        json j = json::parse(z);
        
        max_it = j.value("max_it", 1000);
        tol_res = j.value("tol_res", 1e-7);
        forceFunction = j.value("forceFunction"," ");
        gridDim  = j.value("gridDimension", 10);
    }

    MPI_Bcast(&max_it, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tol_res, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&gridDim, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&forceFunction, 1, MPI_CHAR, 0, MPI_COMM_WORLD);

    std::vector<double> f_tot( gridDim * (gridDim) , 0.0);
    int local_nrows = (gridDim % size > rank) ? gridDim/size +1 : gridDim/size;
    std::vector<double> f_local( gridDim * local_nrows , 0.0);
    std::vector<int> counts(size-1, 0);
    std::vector<int> displs(size-1, 0);
    std::vector<double> U_real_local( gridDim * (local_nrows)  , 0.0);
    std::vector<double> U_real_tot( gridDim * (gridDim) , 0.0);
    
    jacobi_parallel( gridDim, forceFunction, f_tot, f_local, counts, displs, U_real_local, U_real_tot, max_it, tol_res );
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank == 0){
        //std::cout << "\nU_real_tot dopo iterazioni:\n";
        //print_grids(U_real_tot, gridDim, gridDim);
        write_vtks("solution.vtk", U_real_tot, gridDim, gridDim);
    }

    MPI_Finalize();


    return 0;
}