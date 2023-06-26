# heat_conduction
A simple example of the heat conduction equation solved using the PETSc library.

# HowTo
```bash
$ cmake .
$ make
$ ./heat_conduction -d 2 -t 0 -n 101 101
```
where
```bash
    -d - number of dimensions (2 or 3)" << std::endl;
    -t - is problem steady-state (0) or transient (1)" << std::endl;
    -n - number of grid points in each dimension" << std::endl;
```
