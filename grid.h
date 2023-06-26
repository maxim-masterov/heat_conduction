//
// Created by Maxim Masterov on 24/01/2023.
//

#ifndef TRANSIENT_DIFFUSION_GRID_H
#define TRANSIENT_DIFFUSION_GRID_H

#include <iostream>
#include <mpi.h>

class Grid {
public:
    Grid(double _L, double _H, int _NI_nodes, int _NJ_nodes) {
        set2D(_L, _H, _NI_nodes, _NJ_nodes);

        D = 0.;
        z_cells = 0;
        dz = 0.;
        dz2 = 0.;
        dz2_rec = 0.;

        is_grid_2d = true;
    }

    Grid(double _L, double _H, double _D, int _NI_nodes, int _NJ_nodes, int _NK_nodes) {
        set2D(_L, _H, _NI_nodes, _NJ_nodes);

        D = _D;
        z_cells = _NK_nodes;
        dz = D / z_cells;
        dz2 = dz * dz;
        dz2_rec = 1. / dz2;

        is_grid_2d = false;
    }

    ~Grid() { }

    inline double get_dx() const {
        return dx;
    }

    inline double get_dy() const {
        return dy;
    }

    inline double get_dz() const {
        return dz;
    }

    inline double get_dx2_rec() const {
        return dx2_rec;
    }

    inline double get_dy2_rec() const {
        return dy2_rec;
    }

    inline double get_dz2_rec() const {
        return dz2_rec;
    }

    inline int get_x_cells() const {
        return x_cells;
    }

    inline int get_y_cells() const {
        return y_cells;
    }

    inline int get_z_cells() const {
        return z_cells;
    }

    inline int get_total_num_cells() {
        if (is_grid_2d) return x_cells * y_cells;
        else return x_cells * y_cells * z_cells;
    }

    inline bool is_2d() {
        return is_grid_2d;
    }

    inline void report() {
        int pid = 0;

        MPI_Comm_rank(MPI_COMM_WORLD, &pid);

        if (pid == 0) {
            std::cout << "\nGrid data:" << "\n";
            std::cout << "  is 2D:    " << is_grid_2d << "\n";
            std::cout << "  L:        " << L << "\n";
            std::cout << "  H:        " << H << "\n";
            std::cout << "  D:        " << D << "\n";
            std::cout << "  NI_cells: " << x_cells << "\n";
            std::cout << "  NJ_cells: " << y_cells << "\n";
            std::cout << "  NK_cells: " << z_cells << "\n";
            std::cout << "  NI_nodes: " << x_cells + 1 << "\n";
            std::cout << "  NJ_nodes: " << y_cells + 1 << "\n";
            std::cout << "  NK_nodes: " << z_cells + 1 << "\n";
            std::cout << "  dx:       " << dx << "\n";
            std::cout << "  dy:       " << dy << "\n";
            std::cout << "  dz:       " << dz << "\n";
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }

private:
    void set2D(double _L, double _H, int _NI_nodes, int _NJ_nodes) {
        L = _L;
        H = _H;

        x_cells = _NI_nodes - 1;
        y_cells = _NJ_nodes - 1;

        dx = L / x_cells;
        dy = H / y_cells;

        dx2 = dx * dx;
        dy2 = dy * dy;

        dx2_rec = 1. / dx2;
        dy2_rec = 1. / dy2;
    }

private:
    double L;       // domain size in "x" direction (m)
    double H;       // domain size in "y" direction (m)
    double D;       // domain size in "z" direction (m)
    double dx;      // Grid step in "x" direction (m)
    double dy;      // Grid step in "y" direction (m)
    double dz;      // Grid step in "y" direction (m)
    double dx2;     // squared Grid step in "x" direction (m2)
    double dy2;     // squared Grid step in "y" direction (m2)
    double dz2;     // squared Grid step in "y" direction (m2)
    double dx2_rec; // reciprocal squared Grid step in "x" direction (1/m2)
    double dy2_rec; // reciprocal squared Grid step in "y" direction (1/m2)
    double dz2_rec; // reciprocal squared Grid step in "y" direction (1/m2)
    int x_cells;    // number of elements in "x" direction (-)
    int y_cells;    // number of elements in "y" direction (-)
    int z_cells;    // number of elements in "y" direction (-)
    bool is_grid_2d;// "true" if Grid is 2D, "false" if it is 3D
};

#endif //TRANSIENT_DIFFUSION_GRID_H
