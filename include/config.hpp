#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <string>

void read_config(const std::string& filename, int& max_it, double& tol_res, std::string& forceFunction, int& gridDim);

#endif // CONFIG_HPP
