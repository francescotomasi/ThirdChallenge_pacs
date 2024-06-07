
# Third Challenge Pacs

## Project Description

This project implements a Jacobi algorithm for solving linear systems in both parallel and non-parallel modes. The objective is to optimize computational performance by using the MPI library to distribute the workload across multiple processors.

### Compilation
change directory in ThirdChallenge_pacs

```sh
make
```

If you have problem like:
```sh
./main: error while loading shared libraries: libmuparser.so.2: cannot open shared object file: No such file or directory
```

write in your current directory something like that (with your path)
```sh
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/frenci/pacs-examples/Examples/lib 
```

### Esecution

```sh
mpirun -np <numero_di_processori> ./main
```

## Code Structure

### Input Parameters

At the beginning, the following parameters are saved:
- `max_it`: Maximum number of iterations
- `tol_res`: Tolerance for the convergence criterion
- `gridDimension`: Grid dimension
- `function`: Function defined using Muparser

### Data Distribution

1. Calculate `local_nrows` for each rank.
2. `counts`: Counts the points in each rank.
3. `displs`: The starting element for counting.

The initial function values are stored in a vector `f_tot`, which is then scattered into local vectors `f_local` for each rank.

### Jacobi Algorithm

For the algorithm, two new vectors `T_local` and `U_local_new` are created to store the results at iterations \( k \) and \( k+1 \). These vectors have an extra row at the top and bottom to handle boundaries not present in the rank.

1. Execute the Jacobi algorithm.
2. Save the local results in `U_real_local` with the correct dimensions.
3. Gather `U_real_local` into `U_real_tot`, representing the final function.

## Requirements

- MPI (Message Passing Interface)
- Muparser
- C++ Compiler

## Comments

- The compiler will give error but the result and "file vtk" created are correct
