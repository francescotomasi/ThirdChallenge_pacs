#include "grid_utils.hpp"
#include <iostream>
#include <vector>

void print_grid(const std::vector<double>& grid, int row, int col) {
    std::cout << "\nGrid:\n";
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            std::cout << grid[i * col + j] << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}


void write_vtk(const std::string& filename, const std::vector<double>& grid, int row, int col) {
    std::ofstream file;
    file.open(filename);

    file << "# vtk DataFile Version 3.0\n";
    file << "Grid data\n";
    file << "ASCII\n";
    file << "DATASET RECTILINEAR_GRID\n";
    file << "DIMENSIONS " << col << " " << row << " 1\n";
    
    file << "X_COORDINATES " << col << " float\n";
    for (int i = 0; i < col; ++i) {
        file << i << " ";
    }
    file << "\n";
    
    file << "Y_COORDINATES " << row << " float\n";
    for (int i = 0; i < row; ++i) {
        file << i << " ";
    }
    file << "\n";
    
    file << "Z_COORDINATES 1 float\n";
    file << "0\n";
    
    file << "POINT_DATA " << row * col << "\n";
    file << "SCALARS solution float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            file << grid[i * col + j] << " ";
        }
        file << "\n";
    }
    
    file.close();
}
