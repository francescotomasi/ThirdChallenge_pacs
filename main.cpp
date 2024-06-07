#include <mpi.h>
#include <omp.h>
#include <vector>
#include <nlohmann/json.hpp>
#include "muparser_fun.hpp"
#include <fstream>

using json = nlohmann::json;

void print_grid( const std::vector<double>& grid, int row, int col) {
    std::cout << "\nGrid:\n";
    for ( int i = 0; i < row; ++i) {
        for ( int j = 0; j < col; ++j) {
            std::cout << grid[i * col + j] << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
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

    double h = 1.0 / (gridDim - 1);
    MuparserFun fun(forceFunction, 2);
    std::vector<double> f_tot( gridDim * (gridDim) , 0.0);
    int local_nrows = (gridDim % size > rank) ? gridDim/size +1 : gridDim/size;
    std::vector<double> f_local( gridDim * local_nrows , 0.0);
    std::vector<int> counts(size-1, 0);
    std::vector<int> displs(size-1, 0);

    if (rank == 0){
        // Update f_tot with the function values
        for ( int i = 1; i < gridDim-1  ; ++i) {
            for ( int j = 1; j < gridDim-1 ; ++j) {
                double x = (rank * gridDim + i) * h;
                double y = j * h;
                f_tot[i * gridDim + j] = fun(x, y);
            }
        }
        // Calculate counts and displacements
        for (int i = 0; i < size; i++){
            int local_r = (gridDim % size > i) ? gridDim/size +1 : gridDim/size;
            counts[i] = gridDim * local_r;
            if (i != size-1){
                displs[i+1] = displs[i] + counts[i];
            }
        }
    }
   
    MPI_Bcast(counts.data(), size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(displs.data(), size, MPI_INT, 0, MPI_COMM_WORLD);
 
    MPI_Scatterv( f_tot.data(), counts.data(), displs.data(), MPI_DOUBLE, 
                  f_local.data(), gridDim * local_nrows, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0){
        std::cout << "Force function da rank 0:\n";
        print_grid(f_local, local_nrows, gridDim);
    }

    //per algoritmo
    std::vector<double> U_local( gridDim * (local_nrows+2)  , 0.0);
    std::vector<double> U_local_new( gridDim * (local_nrows+2)  , 0.0);
   
    double err= 0.0;
    double err_loc=0.0;
    int k = 0;
    do{
        MPI_Barrier(MPI_COMM_WORLD);

        if(rank > 0){
            MPI_Sendrecv( U_local.data()+ gridDim, gridDim, MPI_DOUBLE, rank -1, 0, //invio a quello prima la seconda riga di local_new
                          U_local.data(), gridDim, MPI_DOUBLE, rank- 1, 1,  //ricevo da quello prima la sua penultima colonna
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
        }
        if(rank < size -1 ){
            MPI_Sendrecv( U_local.data() + (local_nrows)*gridDim, gridDim, MPI_DOUBLE, rank +1, 1, //invio a quello dopo la penultima riga di local_new
                          U_local.data()+ (local_nrows+1)*gridDim, gridDim, MPI_DOUBLE, rank+ 1, 0,  //ricevo da quello dopo la sua seconda colonna nella mia ultima riga
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        err_loc = 0.0;  
        MPI_Barrier(MPI_COMM_WORLD);   

        // Update each internal entry
        if(rank == 0){
            #pragma omp parallel for 
            for(int i = 2; i< local_nrows+1; i++){
                for ( int j = 1; j < gridDim - 1; ++j) {
                    U_local_new[i * gridDim + j] = 0.25 * (  U_local[(i - 1) * gridDim + j] + U_local[(i + 1) * gridDim + j] +
                                                             U_local[i * gridDim + (j - 1)] + U_local[i * gridDim + (j + 1)] +
                                                             h * h * f_local[(i-1)* gridDim + j]);
                    err_loc += (U_local_new[i * gridDim + j] - U_local[i * gridDim + j]) * (U_local_new[i * gridDim + j] - U_local[i * gridDim + j]);
                }
            }
        }
        else if( rank == size-1){
            #pragma omp parallel for 
            for(int i = 1; i< local_nrows; i++){
                for ( int j = 1; j < gridDim - 1; ++j) {
                    U_local_new[i * gridDim + j] = 0.25 * (  U_local[(i - 1) * gridDim + j] + U_local[(i + 1) * gridDim + j] +
                                                             U_local[i * gridDim + (j - 1)] + U_local[i * gridDim + (j + 1)] +
                                                             h * h * f_local[(i-1)* gridDim + j]);
                    err_loc += (U_local_new[i * gridDim + j] - U_local[i * gridDim + j]) * (U_local_new[i * gridDim + j] - U_local[i * gridDim + j]);
                }
            }
        }
        else {
            #pragma omp parallel for 
            for(int i = 1; i< local_nrows+1; i++){
                for ( int j = 1; j < gridDim - 1; ++j) {
                    U_local_new[i * gridDim + j] = 0.25 * (  U_local[(i - 1) * gridDim + j] + U_local[(i + 1) * gridDim + j] +
                                                             U_local[i * gridDim + (j - 1)] + U_local[i * gridDim + (j + 1)] +
                                                             h * h * f_local[(i-1)* gridDim + j]);
                    err_loc += (U_local_new[i * gridDim + j] - U_local[i * gridDim + j]) * (U_local_new[i * gridDim + j] - U_local[i * gridDim + j]);                
                }
            }
        }
        
        k = k+1;
        //update error
        MPI_Allreduce(&err_loc, &err, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        err = sqrt(err*h);

        MPI_Barrier(MPI_COMM_WORLD);
        U_local = U_local_new;
        MPI_Barrier(MPI_COMM_WORLD);

    }while(k< max_it && err > tol_res);

    MPI_Barrier(MPI_COMM_WORLD);
    if(rank == 1){
        std::cout << "\nU_local da rank 1 dopo " << k << " iterazioni:\n";
        print_grid(U_local, local_nrows+2, gridDim);
    }
    std::vector<double> U_real_local( gridDim * (local_nrows)  , 0.0);
    std::vector<double> U_real_tot( gridDim * (gridDim) , 0.0);
   
    for (int i = 1; i < local_nrows+1; i++){
        for (int j = 1; j < gridDim - 1; j++){
            U_real_local[(i-1)*gridDim + j] = U_local[i*gridDim + j];
        }
    }

    if(rank == 1){
        std::cout << "\nU_real_local da rank 1 dopo " << k << " iterazioni:\n";
        print_grid(U_real_local, local_nrows, gridDim);
    }
    MPI_Gatherv(U_real_local.data() , gridDim * local_nrows, MPI_DOUBLE, 
                U_real_tot.data(), counts.data(), displs.data() , MPI_DOUBLE, 0,
                MPI_COMM_WORLD);
    
    if(rank == 0){
        std::cout << "\nU_real_tot dopo " << k << " iterazioni:\n";
        print_grid(U_real_tot, gridDim, gridDim);
    }

    MPI_Finalize();


    return 0;
}