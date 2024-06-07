#include "algorithms.hpp"
#include "muparser_fun.hpp"
#include <vector>
#include <iostream>
#include <cmath>
#include <mpi.h>

void jacobi_parallel( int gridDim, std::string& forceFunction, std::vector<double>& f_tot, 
                      std::vector<double>& f_local, std::vector<int>& counts, 
                      std::vector<int>& displs, std::vector<double>& U_real_local, 
                      std::vector<double>& U_real_tot, int max_it, double tol_res ) {
    MuparserFun fun(forceFunction, 2);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int local_nrows = (gridDim % size > rank) ? gridDim/size +1 : gridDim/size;
    double h = 1.0 / (gridDim - 1);

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

  MPI_Bcast(counts.data(), size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(displs.data(), size, MPI_INT, 0, MPI_COMM_WORLD);
 
    MPI_Scatterv( f_tot.data(), counts.data(), displs.data(), MPI_DOUBLE, 
                  f_local.data(), gridDim * local_nrows, MPI_DOUBLE, 0, MPI_COMM_WORLD);

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
   
    for (int i = 1; i < local_nrows+1; i++){
        for (int j = 1; j < gridDim - 1; j++){
            U_real_local[(i-1)*gridDim + j] = U_local[i*gridDim + j];
        }
    }

    MPI_Gatherv(U_real_local.data() , gridDim * local_nrows, MPI_DOUBLE, 
                U_real_tot.data(), counts.data(), displs.data() , MPI_DOUBLE, 0,
                MPI_COMM_WORLD);
    
    MPI_Barrier(MPI_COMM_WORLD);                     
                    
                      
}
    


void run_parallel_algorithm(int rank, int size, int gridDim, int local_nrows, int max_it, double tol_res, 
                            std::vector<double>& f_local, std::vector<double>& U_local, 
                            std::vector<double>& U_local_new, std::vector<int>& counts, 
                            std::vector<int>& displs) {
    double err = 0.0;
    double err_loc = 0.0;
    int k = 0;
    int h = 1.0 / (gridDim - 1);
    do {
        // Parallel algorithm logic here
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


    } while (k < max_it && err > tol_res);
}

void run_sequential_algorithm(int gridDim, int max_it, double tol_res, 
                              std::vector<double>& f, std::vector<double>& U, 
                              std::vector<double>& U_new) {
    double err = 0.0;
    int k = 0;
    double h= 1.0 / (gridDim - 1);
    do {
        err = 0;
        for (int i = 0; i < gridDim ; ++i) {
            for (int j = 0; j < gridDim ; ++j) {
                if(i == 0 || i == gridDim-1 || j == 0 || j == gridDim-1)
                    U_new[i*gridDim + j] = 0.0;
                else
                    U_new[i*gridDim + j] = 0.25 * (  U[(i - 1) * gridDim + j] + U[(i + 1) * gridDim + j] +
                                                    U[i * gridDim + (j - 1)] + U[i * gridDim + (j + 1)] +
                                                    h * h * f[i*gridDim + j]);
            
                err += (U_new[i*gridDim + j] - U[i*gridDim + j]) * (U_new[i*gridDim + j] - U[i*gridDim + j]);
            }
        }
        err = sqrt(err*h);

        U = U_new;
        k= k + 1;
        

        // Sequential algorithm logic here
    } while (k < max_it && err > tol_res);
}
