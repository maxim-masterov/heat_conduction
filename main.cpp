#include <iostream>
#include <mpi.h>
#include <petsc.h>
#include <string>
#include <algorithm>

#include "equation.h"

struct InputData {
    int dims;
    bool is_transient;
    int nodes[3];
};

bool isNumber(const std::string& s)
{
    return !s.empty() && std::find_if(s.begin(),
                                      s.end(),
                                      [](unsigned char c) { return !std::isdigit(c); }) == s.end();
}

void terminateDueToParserFailure() {
    int pid = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);

    if (pid == 0) {
        std::cerr << "Error! Invalid input parameters." << std::endl;
        std::cerr << "Usage example:" << std::endl;
        std::cerr << "    ./heat_conduction -d 2 -t 0 -n 101 101" << std::endl;
        std::cerr << "Options:" << std::endl;
        std::cerr << "    -d - number of dimensions (2 or 3)" << std::endl;
        std::cerr << "    -t - is problem steady-state (0) or transient (1)" << std::endl;
        std::cerr << "    -n - number of grid points in each dimension" << std::endl;
    }
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
}

InputData parseInput(int argc, char** argv) {
    InputData input;
    bool dims_knowns = false;

    /* Assign the default values first. */
    input.dims = 2;
    input.is_transient = false;
    input.nodes[0] = 51;
    input.nodes[1] = 51;
    input.nodes[2] = -1;

    if (argc < 8 || argc > 9) {
        terminateDueToParserFailure();
    } else {
        for (int n = 1; n < argc; ++n) {
            if (std::string(argv[n]) == "-d") {
                input.dims = atoi(argv[n + 1]);
                if (input.dims != 2 && input.dims != 3) {
                    terminateDueToParserFailure();
                }
                dims_knowns = true;
                ++n;
            } else if (std::string(argv[n]) == "-t") {
                input.is_transient = (bool)atoi(argv[n + 1]);
                ++n;
            } else if (std::string(argv[n]) == "-n") {
                input.nodes[0] = atoi(argv[n + 1]);
                input.nodes[1] = atoi(argv[n + 2]);
                if (dims_knowns) {
                    if (input.dims == 3) {
                        input.nodes[2] = atoi(argv[n + 3]);
                        ++n;
                    }
                }
                else {
                    if (isNumber(argv[n + 3])) {
                        input.nodes[2] = atoi(argv[n + 3]);
                        ++n;
                    }
                }
                n += 2;
            }
        }

        if (input.dims == 3 && input.nodes[2] == -1) {
            terminateDueToParserFailure();
        }
    }

    return input;
}

void run(int dims, bool is_transient, int nodes[3]) {
    Grid *grid;
    Equation *equation;
    int num_procs = 1;
    int pid = 0;

    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);

    if (dims == 2)
        grid = new Grid(0.1, 0.1, nodes[0], nodes[1]);
    else if (dims == 3)
        grid = new Grid(0.1, 0.1, 0.1, nodes[0], nodes[1], nodes[2]);
    else {
        if (pid == 0)
            std::cout << "Incorrect number of dimensions. Should be '2' or '3', not " << dims << "\n";
        return;
    }

    if (is_transient)
        equation = new Equation(grid, 5e-5, 15.);
    else
        equation = new Equation(grid, 15.);

    if (num_procs > 1)
    {
        grid->report();
        equation->solve();
    }
    else {
        std::cout << "This program requires more than one MPI rank" << "\n";
    }

    delete grid;
    delete equation;
}

int main(int argc, char** argv) {
    std::string help = "Transient diffusion problem.\n\n";

    PetscInitialize(&argc, &argv,(char*)0, help.c_str());

    InputData input = parseInput(argc, argv);

    run(input.dims, input.is_transient, input.nodes);

    PetscFinalize();
    return 0;
}
