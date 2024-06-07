#ifndef ALGORITHMS_HPP
#define ALGORITHMS_HPP

#include <vector>
#include <string>

void jacobi_parallel( int gridDim, std::string& forceFunction, std::vector<double>& f_tot, 
                      std::vector<double>& f_local, std::vector<int>& counts, 
                      std::vector<int>& displs, std::vector<double>& U_real_local, 
                      std::vector<double>& U_real_tot, int max_it, double tol_res );
    

#endif // ALGORITHMS_HPP
