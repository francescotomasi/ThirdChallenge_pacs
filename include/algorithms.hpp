#ifndef ALGORITHMS_HPP
#define ALGORITHMS_HPP

#include <vector>

void run_parallel_algorithm(int rank, int size, int gridDim, int local_nrows, int max_it, double tol_res, 
                            std::vector<double>& f_local, std::vector<double>& U_local, 
                            std::vector<double>& U_local_new, std::vector<int>& counts, 
                            std::vector<int>& displs);

void run_sequential_algorithm(int gridDim, int max_it, double tol_res, 
                              std::vector<double>& f_local, std::vector<double>& U_local, 
                              std::vector<double>& U_local_new);

#endif // ALGORITHMS_HPP
