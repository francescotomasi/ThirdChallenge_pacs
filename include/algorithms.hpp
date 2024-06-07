#ifndef ALGORITHMS_HPP
#define ALGORITHMS_HPP

#include <vector>
#include <string>

void run_parallel_algorithm(int rank, int size, int gridDim, int local_nrows, int max_it, double tol_res, 
                            std::vector<double>& f_local, std::vector<double>& U_local, 
                            std::vector<double>& U_local_new, std::vector<int>& counts, 
                            std::vector<int>& displs);

void run_sequential_algorithm(int gridDim, int max_it, double tol_res, 
                              std::vector<double>& f, std::vector<double>& U, 
                              std::vector<double>& U_new);

void jacobi_parallel( int gridDim, std::string& forceFunction, std::vector<double>& f_tot, 
                      std::vector<double>& f_local, std::vector<int>& counts, 
                      std::vector<int>& displs, std::vector<double>& U_real_local, 
                      std::vector<double>& U_real_tot, int max_it, double tol_res );
    

#endif // ALGORITHMS_HPP
