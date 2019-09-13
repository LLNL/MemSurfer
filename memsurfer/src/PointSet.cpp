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

#include <sstream>
#include <stdexcept>

#include "PointSet.hpp"
#include "TriMesh.hpp"

//! ----------------------------------------------------------------------------
//! PointSet constructor
//! ----------------------------------------------------------------------------


Point_with_normal create_point(const TypeFunction &x, const TypeFunction &y, const TypeFunction &z) {
    return Point_with_normal(Point3(x,y,z), Vector3(0,0,0));
}

PointSet::PointSet(float *_, int n, int d) {

    // the dimensionality should be 3
    if (d != 3) {
        std::ostringstream errMsg;
        errMsg << " PointSet::PointSet("<<_<<","<<n<<","<<d<<"): Invalid dimensionality!\n";
        throw std::invalid_argument(errMsg.str());
    }

    valid_normals = false;
    mnPoints = n;
    mPoints.resize(mnPoints);

    for(size_t i = 0; i < mnPoints; i++) {
        mPoints[i] = create_point(_[d*i],  _[d*i+1],  _[d*i+2]);
    }
}
/*
PointSet::PointSet(double *_, int n, int d) {

    // the dimensionality should be 3
    if (d != 3) {
        std::ostringstream errMsg;
        errMsg << " PointSet::PointSet("<<_<<","<<n<<","<<d<<"): Invalid dimensionality!\n";
        throw std::invalid_argument(errMsg.str());
    }

    valid_normals = false;
    mnPoints = n;
    mPoints.resize(mnPoints);

    for(size_t i = 0; i < mnPoints; i++) {
        mPoints[i] = create_point(_[d*i],  _[d*i+1],  _[d*i+2]);
    }
}*/

//! ----------------------------------------------------------------------------
//! get the normals
//! ----------------------------------------------------------------------------

std::vector<TypeFunction> PointSet::get_normals() {

    if (!valid_normals) {
        need_normals(18);
    }

    std::vector<TypeFunction> normals;
    normals.reserve(3*mnPoints);
    size_t i = 0;
    for (auto iter = mPoints.begin(); iter != mPoints.end(); ++iter) {
        const Vector3 &n = iter->second;
        normals.push_back(n[0]);
        normals.push_back(n[1]);
        normals.push_back(n[2]);
        if (++i == mnPoints)
            break;
    }
    return normals;
}

//! ----------------------------------------------------------------------------
//! set periodic box
//! ----------------------------------------------------------------------------

void PointSet::set_periodic(const std::vector<TypeFunction> &box, const TypeFunction thickness, const bool verbose) {

    if(thickness <= 0 || thickness >= 1) {

        std::ostringstream errMsg;
        errMsg << " PointSet::set_periodic(): Invalid thickness of boundary layer (should be 0 < t < 1)! got " << thickness << "!\n";
        throw std::invalid_argument(errMsg.str());
    }

    switch (box.size()) {
        case 3:     mBox0 = std::vector<TypeFunction>({0, 0, 0});                   mBox1 = std::vector<TypeFunction>({box[0], box[1], box[2]});     break;
        case 6:     mBox0 = std::vector<TypeFunction>({box[0], box[1], box[2]});    mBox1 = std::vector<TypeFunction>({box[3], box[4], box[5]});     break;
        default:    {
                    std::ostringstream errMsg;
                    errMsg << " PointSet::set_periodic(): Invalid periodic box (expected [0,0,0 -- b0,b1,b2] or [b0,b1,b2 -- b3,b4,b5])! got " << box.size() << " values!" << std::endl;
                    throw std::invalid_argument(errMsg.str());
                    }
    }

    // width of the box
    const TypeFunction boxw[2] = {mBox1[0]-mBox0[0], mBox1[1]-mBox0[1]};

    // cutoffs for the points to duplicate
    const TypeFunction lcutoff[2] = { mBox0[0] + thickness*boxw[0], mBox0[1] + thickness*boxw[1] };
    const TypeFunction rcutoff[2] = { mBox1[0] - thickness*boxw[0], mBox1[1] - thickness*boxw[1] };

    // start duplicating
    auto eiter = mPoints.end();
    for (auto iter = mPoints.begin(); iter != eiter; ++iter) {
        const Point3 &p = iter->first;
        bool lx = p[0] < lcutoff[0];
        bool ly = p[1] < lcutoff[1];

        bool rx = p[0] > rcutoff[0];
        bool ry = p[1] > rcutoff[1];

        if(!lx && !rx && !ly && !ry)
            continue;

            if (lx) {          mPoints.push_back(create_point(p[0]+boxw[0], p[1], p[2]));          }
       else if (rx) {          mPoints.push_back(create_point(p[0]-boxw[0], p[1], p[2]));          }

            if (ly) {          mPoints.push_back(create_point(p[0], p[1]+boxw[1], p[2]));          }
       else if (ry) {          mPoints.push_back(create_point(p[0], p[1]-boxw[1], p[2]));          }

            if (lx && ly) {    mPoints.push_back(create_point(p[0]+boxw[0], p[1]+boxw[1], p[2]));   }
       else if (lx && ry) {    mPoints.push_back(create_point(p[0]+boxw[0], p[1]-boxw[1], p[2]));   }
       else if (rx && ly) {    mPoints.push_back(create_point(p[0]-boxw[0], p[1]+boxw[1], p[2]));   }
       else if (rx && ry) {    mPoints.push_back(create_point(p[0]-boxw[0], p[1]-boxw[1], p[2]));   }
    }

    if (verbose)
        std::cout << "   > PointSet::set_periodic(" << thickness << ") duplicated " << mPoints.size()-mnPoints << " points!\n";
}

//! ----------------------------------------------------------------------------
//! compute point normals
//! ----------------------------------------------------------------------------

#include <CGAL/pca_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>

std::vector<TypeFunction> PointSet::need_normals(const TypeIndex nb_neighbors, const bool verbose) {

    if (verbose) {
        std::cout << "   > PointSet::need_normals<" << mPoints.size() <<">("<<nb_neighbors<<")...";
        fflush(stdout);
    }

       CGAL::pca_estimate_normals<Concurrency_tag> (mPoints, nb_neighbors,
                                            CGAL::parameters::point_map(CGAL::First_of_pair_property_map<Point_with_normal>()).
                                                              normal_map(CGAL::Second_of_pair_property_map<Point_with_normal>()));

       // Orients normals.
       // Note: mst_orient_normals() requires a range of points
       // as well as property maps to access each point's position and normal.
       //std::list<PointVectorPair>::iterator
       auto unoriented_points_begin =
       CGAL::mst_orient_normals(mPoints, nb_neighbors,
                                CGAL::parameters::point_map(CGAL::First_of_pair_property_map<Point_with_normal>()).
                                                  normal_map(CGAL::Second_of_pair_property_map<Point_with_normal>()));

       // check if any normals failed
    size_t nfailed = std::distance(unoriented_points_begin, mPoints.end());
    if (nfailed > 0) {
        if(verbose)
            std::cout << "\n";

        std::cerr << "   > PointSet::need_normals<" << mPoints.size() <<">("<<nb_neighbors<<") failed for " << nfailed << " points\n";

        if (nfailed >= mPoints.size()-mnPoints) {
            std::cerr << "     > Warning: Some normals may be incorrectly oriented!\n";

            /*
            // Optional: delete points with an unoriented normal
                // if you plan to call a reconstruction algorithm that expects oriented normals.
            mpoints.erase(unoriented_points_begin, mpoints.end());
            return std::vector<TypeFunction> ();*/
        }
    }

    if(verbose)
        std::cout << " Done!\n";

    valid_normals = true;
    return get_normals();
}

//! ----------------------------------------------------------------------------
//! compute Poisson surface
//! ----------------------------------------------------------------------------
#ifdef CPP_POISSON
#include <CGAL/compute_average_spacing.h>
#include <CGAL/Poisson_reconstruction_function.h>

#include <CGAL/Implicit_surface_3.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>

#include <CGAL/make_surface_mesh.h>

//#include <CGAL/Polygon_mesh_processing/distance.h>
//#include <CGAL/IO/Polyhedron_iostream.h>
//#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/output_surface_facets_to_polyhedron.h>

typedef Kernel::Sphere_3 Sphere3;

typedef CGAL::Surface_mesh_default_triangulation_3 STriang3;
typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<STriang3> STriang32;

typedef CGAL::Poisson_reconstruction_function<Kernel> Poisson_reconstruction_function;
typedef CGAL::Implicit_surface_3<Kernel, Poisson_reconstruction_function> ImpSurface_3;

TriMesh* PointSet::need_approximate_surface(bool verbose) {

    if (verbose) {
        std::cout << "   > PointSet::need_approximate_surface()...";
        fflush(stdout);
    }

    // Poisson options
    const FT sm_angle = 20.0;       // Min triangle angle in degrees.
    const FT sm_radius = 30.0;      // Max triangle size w.r.t. point set average spacing.
    const FT sm_distance = 0.375;   // Surface Approximation error w.r.t. point set average spacing.

    const FT radius_factr = 2.0;
    const FT error_factr = 0.001;

    // ------------------------------------------------------------------------------
    // Computes average spacing between points
    FT average_spacing = CGAL::compute_average_spacing<Concurrency_tag>(mPoints.begin(), mPoints.end(), 6 /* knn = 1 ring */);

    // ------------------------------------------------------------------------------
    // Creates implicit function from the read points using the default solver.

    // Note: this method requires an iterator over points
    // + property maps to access each point's position and normal.
    // The position property map can be omitted here as we use iterators over Point_3 elements.
    Poisson_reconstruction_function function(mPoints.begin(), mPoints.end(),
                                             CGAL::make_normal_of_point_with_normal_pmap(std::vector<Point_with_normal>::value_type()) );

    // Computes the Poisson indicator function f()
    // at each vertex of the triangulation.
    if (!function.compute_implicit_function()) {
        std::ostringstream errMsg;
        errMsg << " PointSet::need_approximate_surface(): Failed to compute implicit function!\n";
        throw std::invalid_argument(errMsg.str());
    }

    // ------------------------------------------------------------------------------
    // Gets one point inside the implicit surface
    // and computes implicit function bounding sphere radius.

    Sphere3 bsphere = function.bounding_sphere();
    Point3 inner_point = function.get_inner_point();

    FT radius = std::sqrt(bsphere.squared_radius());

    // Defines the implicit surface: requires defining a
    // conservative bounding sphere centered at inner point.
    FT sm_sphere_radius = radius_factr * radius;
    FT sm_dichotomy_error = sm_distance*average_spacing*error_factr; // Dichotomy error must be << sm_distance

    /*
    std::cout << "avg spacing = " << average_spacing << std::endl;
    std::cout << "bsphere = " << bsphere << std::endl;
    std::cout << "radius = " << radius << std::endl;
    std::cout << "inner point = " << inner_point << std::endl;
    std::cout << "sm_sphere_radius = " << sm_sphere_radius << std::endl;
    std::cout << "sm_dichotomy_error = " << sm_dichotomy_error << std::endl;
    */

    ImpSurface_3 surface(function, Sphere3(inner_point, sm_sphere_radius*sm_sphere_radius),
                                        sm_dichotomy_error/sm_sphere_radius);

    // Defines surface mesh generation criteria
    CGAL::Surface_mesh_default_criteria_3<STriang3> criteria(sm_angle,                      // Min triangle angle (degrees)
                                                             sm_radius*average_spacing,     // Max triangle size
                                                             sm_distance*average_spacing);  // Approximation error

    // Generates surface mesh with manifold option
    STriang3 str3;                  // 3D Delaunay triangulation for surface mesh generation
    STriang32 approx_mesh (str3);   // 2D complex in 3D Delaunay triangulation
    CGAL::make_surface_mesh(approx_mesh,                          // reconstructed mesh
                            surface,                              // implicit surface
                            criteria,                             // meshing criteria
                            CGAL::Manifold_with_boundary_tag());  // require manifold mesh

    if(str3.number_of_vertices() == 0){
        std::ostringstream errMsg;
        errMsg << " PointSet::need_approximate_surface(): Failed to compute approximate surface!\n";
        throw std::invalid_argument(errMsg.str());
    }

    // ------------------------------------------------------------------------------
    Polyhedron surface_mesh;
    CGAL::output_surface_facets_to_polyhedron(approx_mesh, surface_mesh);

    return new TriMesh(surface_mesh);
}
#endif
//! ----------------------------------------------------------------------------
//! ----------------------------------------------------------------------------
