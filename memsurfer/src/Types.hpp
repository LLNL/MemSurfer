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

#ifndef _TYPES_H_
#define _TYPES_H_

/// ----------------------------------------------------------------------------
/// Basic data types
/// ----------------------------------------------------------------------------

#include <cstdint>

typedef uint32_t TypeIndex;
typedef int32_t TypeIndexI;
typedef float TypeFunction;

/// ----------------------------------------------------------------------------
/// Main Geometric Primitives
/// ----------------------------------------------------------------------------

#include "Vec.h"

typedef Vec<3, TypeIndexI> Offset3;
typedef Vec<3, TypeFunction> Vertex;
typedef Vec<2, TypeIndex> Edge;
typedef Vec<3, TypeIndex> Face;
typedef Vertex Normal;

/// ----------------------------------------------------------------------------
/// CGAL data types!
/// ----------------------------------------------------------------------------
//#define CPP_POISSON       // CPP poisson still not ported to cgal 4.13
//#define CPP_REMESHING     // CPP rmeshing still not fully functional

#ifdef CGAL_AVAILABLE
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT          FT;
typedef Kernel::Point_2     Point2;
typedef Kernel::Point_3     Point3;
typedef Kernel::Vector_3    Vector3;
typedef Kernel::Triangle_3  Triangle3;
typedef std::pair<Point2, size_t> Point_with_idx;
//typedef std::pair<Point3, size_t> Point_with_idx;
typedef std::pair<Point3, Vector3> Point_with_normal;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

// Concurrency
#ifdef CGAL_LINKED_WITH_TBB
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif
#endif

/// ----------------------------------------------------------------------------
/// ----------------------------------------------------------------------------
#endif  /* _TYPES_H_ */
