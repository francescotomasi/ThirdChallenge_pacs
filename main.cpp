#include <mpi.h>
#include <omp.h>
#include <vector>
#include <nlohmann/json.hpp>
#include "muparser_fun.hpp"
#include <fstream>

using json = nlohmann::json;

void print_grid( const std::vector<double>& grid, int row, int col) {
    for ( int i = 0; i < row; ++i) {
        for ( int j = 0; j < col; ++j) {
            std::cout << grid[i * col + j] << " ";
        }
        std::cout << "\n";
    }
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
    std::vector<double> U_tot( gridDim * (gridDim) , 0.0);
    int local_nrows = (gridDim % size > rank) ? gridDim/size +1 : gridDim/size;
    std::vector<double> U_local( gridDim * local_nrows , 0.0);

    if (rank == 0){
        // Update U_tot with the function values
        for ( int i = 1; i < gridDim-1  ; ++i) {
            for ( int j = 1; j < gridDim-1 ; ++j) {
                double x = (rank * gridDim + i) * h;
                double y = j * h;
                U_tot[i * gridDim + j] = fun(x, y);
            }
        }
    }

    // Calculate counts and displacements
    std::vector<int> counts(size-1, 0);
    std::vector<int> displs(size-1, 0);
    if (rank == 0){
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
    MPI_Scatterv( U_tot.data(), counts.data(), displs.data(), MPI_DOUBLE, 
                  U_local.data(), gridDim * local_nrows, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // fino a qua ho mandato bene i count, displs, U_local

    //ora faccio un nuovo vettore per tenere conto delle condizioni al contorno
    std::vector<double> U_new_local( gridDim * (local_nrows+2)  , 0.0);

    for(size_t i = gridDim ; i < U_new_local.size() - gridDim ; i++){
            U_new_local[i] = U_local[i - gridDim]; 
        }

    //aggiornati correttamente U_new_local

    
    //creo una copia di U_new_local, mi servirà dell'algoritmo di Jacobi
    std::vector<double> U_new_local_copy= U_new_local;

    //copiato correttamente U_new_local in U_new_local_copy

    //jacobi implemetation
    double e;
    int k = 0;
    do{
        if(rank > 0){
            MPI_Sendrecv( U_new_local.data()+ gridDim, gridDim, MPI_DOUBLE, rank -1, 0, //invio a quello prima la seconda riga di local_new
                          U_new_local.data(), gridDim, MPI_DOUBLE, rank- 1, 1,  //ricevo da quello prima la sua penultima colonna
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
        }
        if(rank < size -1 ){
            MPI_Sendrecv( U_new_local.data() + (local_nrows)*gridDim, gridDim, MPI_DOUBLE, rank +1, 1, //invio a quello dopo la penultima riga di local_new
                          U_new_local.data()+ (local_nrows+1)*gridDim, gridDim, MPI_DOUBLE, rank+ 1, 0,  //ricevo da quello dopo la sua seconda colonna nella mia ultima riga
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        
        // Update each internal entry
        //#pragma omp parallel for
        for ( int i = 1; i < local_nrows + 1; ++i) {
            for ( int j = 1; j < gridDim - 1; ++j) {
                U_new_local_copy[i * gridDim + j] = 0.25 * (U_new_local[(i - 1) * gridDim + j] + U_new_local[(i + 1) * gridDim + j] +
                                                 U_new_local[i * gridDim + (j - 1)] + U_new_local[i * gridDim + (j + 1)] +
                                                 h * h * U_new_local_copy[i * gridDim + j]);  //  f_start[(rank*local_nrows+i) *h][j*h]);  //((rank * local_nrows + i) * h, j * h));
            }
        }
        // setting boundary conditions
        for (int j=0; j<gridDim; j++){
            U_new_local_copy[j] = 0.0;
            U_new_local_copy[(local_nrows+1)*gridDim + j] = 0.0;
        }

        U_new_local = U_new_local_copy;

        if(rank == 1 ){
        std:: cout<< "\nmatrice U_new_local in rank 1 dopo la prima iterazione" << std::endl;
        print_grid(U_new_local, local_nrows+2, gridDim);
        }
        
        // Synchronize all processes
         MPI_Barrier(MPI_COMM_WORLD);

        k=1111;
    }while(k< max_it);
    
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank == 0 ){
        std:: cout<< "\nmatrice U_new_local in rank 0 dopo le iterazioni" << std::endl;
        print_grid(U_new_local, local_nrows+2, gridDim);
    }

    for (int i= 1; i< local_nrows+1; i++){
        for (int j = 0; j< gridDim; j++){
            U_local[i*gridDim + j] = U_new_local[i*gridDim + j];
        }
    }
   MPI_Scatterv( U_tot.data(), counts.data(), displs.data(), MPI_DOUBLE, 
                  U_local.data(), gridDim * local_nrows, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Gatherv(U_local.data() + gridDim, gridDim * local_nrows, MPI_DOUBLE, 
                  U_tot.data(), counts.data(), displs.data() , MPI_DOUBLE, 0,
                  MPI_COMM_WORLD);
/*
    if(rank == 0){ 
        std::cout << "\nFinal grid after iteration: " << k << std::endl;
        print_grid(U_tot, gridDim, gridDim);
    }
*/
    
    MPI_Finalize();

    return 0;
}