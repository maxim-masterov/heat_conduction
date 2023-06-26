//
// Created by Maxim Masterov on 24/01/2023.
//
#include <fstream>

#include "equation.h"

Equation::Equation(Grid *_grid, double _dt, double T0) {
    is_transient = true;
    general_setup(_grid, _dt, T0);
}

Equation::Equation(Grid *_grid, double T0) {
    is_transient = false;
    general_setup(_grid, 1., T0);
}

void Equation::general_setup(Grid *_grid, double _dt, double T0) {
    int pid = 0;

    grid = _grid;
    dt = _dt;
    max_t = 50. * dt;

    T_e = 10.;
    T_w = 20.;

    T_s = 30.;
    T_n = 5.;

    T_t = 40.;
    T_b = 15.;

    alpha = 0.19;  // Air //0.082 * 1e-6;       // pine wood

    freq_file_writing = 10;

    solver_name = KSPCG;
    precond_name = "no";

    create_linear_sys();
    create_field(T0);

    MPI_Comm_rank(MPI_COMM_WORLD, &pid);

    if (pid == 0) {
        std::string equation = "";
        std::cout << "\n" << std::endl;
        if (is_transient) {
            std::cout << "Transient heat conduction problem:" << "\n";
            equation += "dT/dt + ";
        }
        else {
            std::cout << "Steady state heat conduction problem:" << "\n";
        }

        equation += "d2T/dx2 + d2T/dy2";
        if (grid->is_2d())
            equation += " = 0";
        else
            equation += " + d2T/dz2 = 0";

        std::cout << "  " << equation << "\n";

        std::cout << "PETSc solver:         " << solver_name << std::endl;
        std::cout << "PETSc preconditioner: " << precond_name << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

void Equation::create_linear_sys() {
    int off_diag_elts = 4;

    if (!grid->is_2d())
        off_diag_elts = 6;

    MatCreate(PETSC_COMM_WORLD, &A);
    VecCreate(PETSC_COMM_WORLD, &x);
    VecCreate(PETSC_COMM_WORLD, &b);

    MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, grid->get_total_num_cells(), grid->get_total_num_cells());
    VecSetSizes(x, PETSC_DECIDE, grid->get_total_num_cells());
    VecSetSizes(b, PETSC_DECIDE, grid->get_total_num_cells());

    MatSetFromOptions(A);
    VecSetFromOptions(x);
    VecSetFromOptions(b);

    MatMPIAIJSetPreallocation(A, off_diag_elts + 1, NULL, off_diag_elts, NULL);
}

void Equation::create_field(double T0) {
    int loc_num_elts;

    VecGetLocalSize(x, &loc_num_elts);
    field.resize(loc_num_elts, T0);
    field_old.resize(loc_num_elts, T0);
}

void Equation::destroy_field() {
    field.clear();
    field_old.clear();
}

void Equation::destroy_linear_sys() {
    MatDestroy(&A);
    VecDestroy(&x);
    VecDestroy(&b);
}

void Equation::spatial_contrib(const Index ind, Coefficients &coeff) {
    double dx2_rec = grid->get_dx2_rec();
    double dy2_rec = grid->get_dy2_rec();
    double dz2_rec = grid->get_dz2_rec();

    if (grid->is_2d())
        coeff.value.name.center += 2. * alpha * dt * (dx2_rec + dy2_rec);
    else
        coeff.value.name.center += 2. * alpha * dt * (dx2_rec + dy2_rec + dz2_rec);
    coeff.exist.name.center = true;

    if (ind.i == 0) {
        // Dirichlet BC
        coeff.value.name.center += alpha * dt * dx2_rec;
        coeff.value.name.rhs += 2. * alpha * dt * dx2_rec * T_w;
        coeff.exist.name.west = false;
    } else {
        coeff.value.name.west -= alpha * dt * dx2_rec;
        coeff.exist.name.west = true;
    }

    if (ind.i == (grid->get_x_cells() - 1)) {
        // Dirichlet BC
        coeff.value.name.center += alpha * dt * dx2_rec;
        coeff.value.name.rhs += 2. * alpha * dt * dx2_rec * T_e;
        coeff.exist.name.east = false;
    } else {
        coeff.value.name.east -= alpha * dt * dx2_rec;
        coeff.exist.name.east = true;
    }

    if (ind.j == 0) {
        // Neumann BC
//        coeff.value.name.center -= alpha * dt * dy2_rec;
        // Dirichlet BC
        coeff.value.name.center += alpha * dt * dy2_rec;
        coeff.value.name.rhs += 2. * alpha * dt * dy2_rec * T_s;
        coeff.exist.name.south = false;
    } else {
        coeff.value.name.south -= alpha * dt * dy2_rec;
        coeff.exist.name.south = true;
    }

    if (ind.j == (grid->get_y_cells() - 1)) {
        // Neumann BC
//        coeff.value.name.center -= alpha * dt * dy2_rec;
        // Dirichlet BC
        coeff.value.name.center += alpha * dt * dy2_rec;
        coeff.value.name.rhs += 2. * alpha * dt * dy2_rec * T_n;
        coeff.exist.name.north = false;
    } else {
        coeff.value.name.north -= alpha * dt * dy2_rec;
        coeff.exist.name.north = true;
    }

    if (!grid->is_2d()) {
        if (ind.k == 0) {
            // Dirichlet BC
            coeff.value.name.center += alpha * dt * dz2_rec;
            coeff.value.name.rhs += 2. * alpha * dt * dz2_rec * T_b;
            coeff.exist.name.bottom = false;
        } else {
            coeff.value.name.bottom -= alpha * dt * dz2_rec;
            coeff.exist.name.bottom = true;
        }

        if (ind.k == (grid->get_z_cells() - 1)) {
            // Dirichlet BC
            coeff.value.name.center += alpha * dt * dz2_rec;
            coeff.value.name.rhs += 2. * alpha * dt * dz2_rec * T_t;
            coeff.exist.name.top = false;
        } else {
            coeff.value.name.top -= alpha * dt * dz2_rec;
            coeff.exist.name.top = true;
        }
    }
}

void Equation::transient_contrib(const int id, Coefficients &coeff) {
    if (is_transient) {
        coeff.value.name.center += 1.;
        coeff.value.name.rhs += field_old[id];
    }
}

void Equation::assemble() {
    Index ind;
    int dimension = 5;
    std::vector<double> values(dimension);  // Array of values in a row (excluding diagonal)
    std::vector<int> indices(dimension);    // Array of column indices (excluding diagonal)
    int rowStart, rowEnd;

    MatGetOwnershipRange(A, &rowStart, &rowEnd);

    for (int row = rowStart; row < rowEnd; ++row) {
        Coefficients coeff;
        Neighbours<int> ids{};
        double x_val = 0.;

        if (grid->is_2d()) {
            ids.name.center = row;
            ids.name.west  = row - grid->get_y_cells();
            ids.name.east  = row + grid->get_y_cells();
            ids.name.south = row - 1;
            ids.name.north = row + 1;

            ind.set_ind(row, grid->get_y_cells());
        }
        else {
            ids.name.center = row;
            ids.name.west   = row - grid->get_y_cells() * grid->get_z_cells();
            ids.name.east   = row + grid->get_y_cells() * grid->get_z_cells();
            ids.name.south  = row - grid->get_z_cells();
            ids.name.north  = row + grid->get_z_cells();
            ids.name.bottom = row - 1;
            ids.name.top    = row + 1;

            ind.set_ind(row, grid->get_y_cells(), grid->get_z_cells());
        }

        spatial_contrib(ind, coeff);
        transient_contrib(row - rowStart, coeff);

        int num_entries = 0;
        for (int n = 0; n < dimension; ++n) {
            if (coeff.exist.data[n]) {
                indices[num_entries] = ids.data[n];
                values[num_entries] = coeff.value.data[n];
                ++num_entries;
            }
        }

        MatSetValues(A, 1, &row, num_entries, indices.data(), values.data(), INSERT_VALUES);
        VecSetValue(x, row, x_val, INSERT_VALUES);
        VecSetValue(b, row, coeff.value.name.rhs, INSERT_VALUES);
    }
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

    VecAssemblyBegin(x);
    VecAssemblyEnd(x);

    VecAssemblyBegin(b);
    VecAssemblyEnd(b);
}

void Equation::solve() {
    KSP solver;
    PC prec;
    PetscReal residual;
    double time_start, time_end;
    double time = 0.;
    int step = 0;
    int pid = 0;

    MPI_Comm_rank(MPI_COMM_WORLD, &pid);

    // create CG solver with ILU(0) preconditioning
    KSPCreate(PETSC_COMM_WORLD, &solver);
    KSPSetOperators(solver, A, A);

    KSPSetType(solver, solver_name.c_str());

//    // setup AMG as a preconditioner
//    KSPGetPC(solver, &prec);
//    PCSetType(prec, PCGAMG);

    KSPSetFromOptions(solver);
    KSPSetTolerances(solver, 1e-8, 1e-8, 1e+2, 1e3);
    KSPSetUp(solver);

    if (pid == 0) {
        if (is_transient)
            std::cout << "\nTime" << "\t" << "Residual" << "\n";
        else
            std::cout << "\nResidual" << "\n";
    }

    // Time loop
    time_start = MPI_Wtime();
    while(1) {
        std::string file_name = "res_" + std::to_string(step) + ".dat";
        if (!(step % freq_file_writing) || !is_transient)
            write_field_to_file(file_name);

        if (time > max_t) break;
        if (step == 1 && !is_transient) break;

        assemble();

        // iterate
        KSPSolve(solver, b, x);
        KSPGetResidualNorm(solver, &residual);

        copy_solution_to_field();
        copy_field_to_old();

        time += dt;
        ++step;

        if (pid == 0) {
            if (is_transient)
                std::cout << time << "\t" << residual << "\n";
            else
                std::cout << residual << "\n";
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    time_end = MPI_Wtime();

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &time_start, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &time_end, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);


    if (pid == 0)
        std::cout << "Solve time: " << time_end - time_start << "\n";
    MPI_Barrier(MPI_COMM_WORLD);
}

void Equation::write_field_to_file(const std::string file_name) {
    Vec tmp;
    IS is1, is2;
    VecScatter ctx = 0;
    PetscScalar *raw_data;
    int num_elts = grid->get_total_num_cells();
    int pid = 0;

    MPI_Comm_rank(PETSC_COMM_WORLD, &pid);

    /* create two index sets */
    ISCreateStride(PETSC_COMM_SELF, num_elts, num_elts * pid, 1, &is1);
    ISCreateStride(PETSC_COMM_SELF, num_elts, 0, 1, &is2);

    VecCreateSeq(PETSC_COMM_SELF, num_elts, &tmp);
    VecScatterCreateToAll(x, &ctx, &tmp);
    VecScatterBegin(ctx, x, tmp, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(ctx, x, tmp, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterDestroy(&ctx);

    VecGetArray(tmp, &raw_data);

    if (pid == 0) {
        int counter = 0;
        std::ofstream out;
        out.open(file_name);

        if (out.is_open()) {
            for (int i = 0; i < grid->get_x_cells(); ++i) {
                double dx = grid->get_dx() * i + grid->get_dx() * 0.5;
                for (int j = 0; j < grid->get_y_cells(); ++j) {
                    double dy = grid->get_dy() * j + grid->get_dy() * 0.5;
                    if (grid->is_2d()) {
                        out << dx << "\t" << dy << "\t" << raw_data[counter] << "\n";
                        ++counter;
                    }
                    else {
                        for (int k = 0; k < grid->get_z_cells(); ++k) {
                            double dz = grid->get_dz() * k + grid->get_dz() * 0.5;
                            out << dx << "\t" << dy << "\t" << dz << "\t" << raw_data[counter] << "\n";
                            ++counter;
                        }
                    }
                }
            }
            out.close();
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    VecRestoreArray(tmp, &raw_data);

    VecDestroy(&tmp);
    ISDestroy(&is1);
    ISDestroy(&is2);
}

void Equation::copy_field_to_old() {
    field_old = field;
}

void Equation::copy_solution_to_field() {
    PetscScalar *raw_data;
    PetscInt size;

    VecGetLocalSize(x, &size);
    VecGetArray(x, &raw_data);
    std::copy(raw_data, raw_data + size, field.data());
    VecRestoreArray(x, &raw_data);
}