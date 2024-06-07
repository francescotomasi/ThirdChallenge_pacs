#ifndef GRID_UTILS_HPP
#define GRID_UTILS_HPP

#include <vector>
#include <string>
#include <fstream>

void print_grid(const std::vector<double>& grid, int row, int col);

void write_vtk(const std::string& filename, const std::vector<double>& grid, int row, int col);


#endif // GRID_UTILS_HPP
