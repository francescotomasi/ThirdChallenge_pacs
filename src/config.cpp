#include "config.hpp"
#include <fstream>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

void read_config(const std::string& filename, int& max_it, double& tol_res, std::string& forceFunction, int& gridDim) {
    std::ifstream z(filename);
    json j = json::parse(z);

    max_it = j.value("max_it", 1000);
    tol_res = j.value("tol_res", 1e-7);
    forceFunction = j.value("forceFunction", " ");
    gridDim = j.value("gridDimension", 10);
}
