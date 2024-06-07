#ifndef MPI_UTILS_HPP
#define MPI_UTILS_HPP

#include <vector>
#include <string>

class MuparserFun;

void broadcast_parameters(int rank, int& max_it, double& tol_res, std::string& forceFunction, int& gridDim);
void distribute_grid(int rank, int size, int gridDim, MuparserFun& fun, double h, 
                     std::vector<double>& f_local, int& local_nrows, 
                     std::vector<int>& counts, std::vector<int>& displs);
void gather_results(int rank, int size, int gridDim, int local_nrows, 
                    std::vector<double>& U_local, std::vector<double>& U_real_local, 
                    std::vector<double>& U_real_tot, std::vector<int>& counts, 
                    std::vector<int>& displs);

#endif // MPI_UTILS_HPP
