//
// Created by Maxim Masterov on 24/01/2023.
//

#ifndef TRANSIENT_DIFFUSION_EQUATION_H
#define TRANSIENT_DIFFUSION_EQUATION_H

#include <vector>
#include <iomanip>
#include <petsc.h>

#include "grid.h"

template <typename T>
union Neighbours {
    struct Names {
        T center;
        T west;
        T east;
        T north;
        T south;
        T top;
        T bottom;
        T rhs;

        void report() {
            std::cout << "\nCoefficients:" << "\n";
            std::cout << "  center: " << center << "\n";
            std::cout << "  west:   " << west   << "\n";
            std::cout << "  east:   " << east   << "\n";
            std::cout << "  north:  " << north  << "\n";
            std::cout << "  south:  " << south  << "\n";
            std::cout << "  top:    " << top    << "\n";
            std::cout << "  bottom: " << bottom << "\n";
            std::cout << "  rhs:    " << rhs    << "\n";
        }
    } name;

    T data[8] = { };
};

class Coefficients {
public:
    Coefficients() {
        value.data[0] = value.data[1] = value.data[2]
                = value.data[3] = value.data[4] = value.data[5]
                    = value.data[6] = value.data[7] = 0.;
        exist.data[0] = exist.data[1] = exist.data[2]
                = exist.data[3] = exist.data[4] = exist.data[5]
                    = exist.data[6] = exist.data[7] = false;
    }
    ~Coefficients() { }

    Neighbours<double> value;
    Neighbours<bool> exist;

    void report() {
        std::cout << "\nCoefficients:" << "\n";
        std::cout << "  center: " << value.name.center << "\n";
        std::cout << "  west:   " << value.name.west   << "\n";
        std::cout << "  east:   " << value.name.east   << "\n";
        std::cout << "  north:  " << value.name.north  << "\n";
        std::cout << "  south:  " << value.name.south  << "\n";
        std::cout << "  top:    " << value.name.top    << "\n";
        std::cout << "  bottom: " << value.name.bottom << "\n";
        std::cout << "  rhs:    " << value.name.rhs    << "\n";
    }
};

class Index {
public:
    Index() { }
    Index(int _i, int _j) : i(_i), j(_j), k(0) { }
    Index(int _i, int _j, int _k) : i(_i), j(_j), k(_k) { }

    int i;
    int j;
    int k;

    int get_id(const int NJ) const {
        return j + i * NJ;
    }

    int get_id(const int NJ, const int NK) const {
        return k + NK * (j + i * NJ);
    }

    void set_ind(const int id, const int NJ) {
        i = floor(id / NJ);
        j = id - i * NJ;
    }

    void set_ind(const int id, const int NJ, const int NK) {
        i = floor(id / (NJ * NK));
        j = floor((id - i * NJ * NK) / NK);
        k = id - j * NK - i * NJ * NK;
    }
};

class Equation {
public:
    /*!
     * @brief Steady state problem setup
     * @param _grid Object of Grid
     * @param T0 Initial temperature
     */
    Equation(Grid *_grid, double T0);

    /*!
     * @brief Transient problem setup
     * @param _grid Object of Grid
     * @param _dt Time step
     * @param T0 Initial temperature
     */
    Equation(Grid *_grid, double _dt, double T0);

    ~Equation() {
        destroy_field();
        destroy_linear_sys();
    }

    void solve();

private:
    /*!
     * @brief General problem setup
     * @param _grid Object of Grid
     * @param _dt Time step
     * @param T0 Initial temperature
     */
    void general_setup(Grid *_grid, double _dt, double T0);

    /*!
     * @brief Populate matrix coefficients using spatial derivatives.
     * @param ind Global indices of the current cell
     * @param coeff Structure of coefficients
     */
    void spatial_contrib(const Index ind, Coefficients &coeff);

    /*!
     * @param id Local index of the current cell
     * @param coeff Structure of coefficients
     */
    void transient_contrib(const int id, Coefficients &coeff);

    /*!
     * @brief Assemble linear system
     */
    void assemble();

    /*!
     * @brief Write solution to the file
     * @param file_name File name
     */
    void write_field_to_file(const std::string file_name);

    /*!
     * @brief Copy new field to the old one
     */
    void copy_field_to_old();

    /*!
     * @brief Copy solution of the linear system to the field
     */
    void copy_solution_to_field();

    /*!
     * @brief Allocate memory for new and old fields
     * Should be called after the @e create_linear_sys() function
     * @param T0 Initial temperature
     */
    void create_field(double T0);

    /*!
     * @brief Deallocate memory used by new and old fields
     */
    void destroy_field();

    /*!
     * @brief Allocate memory for linear system (A, x, b)
     */
    void create_linear_sys();

    /*!
     * @brief Deallocate memory used by liner system (A, x, b)
     */
    void destroy_linear_sys();

private:
    typedef std::vector<double> Field;

    Field field;                // field from the new time step
    Field field_old;            // field from the old time step
    Grid *grid;                 // numerical Grid
    double dt;                  // time step
    double max_t;               // max simulation time (s)
    double T_e, T_w;            // boundary conditions on the west and east walls respectively (C)
    double T_s, T_n;            // boundary conditions on the south and north walls respectively (C)
    double T_t, T_b;            // boundary conditions on the top and bottom walls respectively (C)
    double alpha;               // thermal diffusivity (m2/s)
    bool is_transient;          // "true" if problem is transient, "false" otherwise
    int freq_file_writing;      // frequency of writing files with the field (ignored for steady state problem)

    // linear system
    Mat A;
    Vec x, b;
    std::string solver_name;    // name of the PETSc solver
    std::string precond_name;   // name of the PETSc solver
};


#endif //TRANSIENT_DIFFUSION_EQUATION_H
