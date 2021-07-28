/**
Copyright (c) 2019, Lawrence Livermore National Security, LLC.
Produced at the Lawrence Livermore National Laboratory.
Written by Harsh Bhatia (hbhatia@llnl.gov) and Peer-Timo Bremer (bremer5@llnl.gov)
LLNL-CODE-763493. All rights reserved.

This file is part of MemSurfer, Version 1.0.
Released under GNU General Public License 3.0.
For details, see https://github.com/LLNL/MemSurfer.
*/

/// ----------------------------------------------------------------------------
/// ----------------------------------------------------------------------------

#ifndef _POINTSET_H_
#define _POINTSET_H_

#include <cstdio>
#include <vector>

#include "Types.hpp"

class TriMesh;

/// ---------------------------------------------------------------------------------------
//!
//! \brief This class provides functionality to operate upon a PointSet
//!         possibly periodic in x,y domain
//!
/// ---------------------------------------------------------------------------------------
class PointSet {

    //! points
    std::vector<Point_with_normal> mPoints;

    //! number of given points
    size_t mnPoints;

    //! bounding box
    std::vector<TypeFunction> mBox0, mBox1;

    //! whether normals are available
    bool valid_normals;

public:
    //! constructor
    PointSet(float *_, int n, int d);
    //PointSet(double *_, int n, int d);

    //! set periodicity
    void set_periodic(const std::vector<TypeFunction> &box, const TypeFunction thickness, const bool verbose = false);

    //! get the normals
    std::vector<TypeFunction> get_normals();

    //! compute the normals
    std::vector<TypeFunction> need_normals(const TypeIndex nb_neighbors, const bool verbose = false);

    //! compute an approximate (Poisson) surface
#ifdef CPP_POISSON
    TriMesh* need_approximate_surface(const bool verbose = false);
#endif
};

/// ---------------------------------------------------------------------------------------
#endif  /* _POINTSET_H_ */
