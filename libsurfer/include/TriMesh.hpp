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

#ifndef _TRIMESH_H_
#define _TRIMESH_H_

#include <cstdio>
#include <tuple>
#include <vector>
#include <unordered_map>

#include "Types.hpp"

class DensityKernel;        // kernel for density estimation
class DistanceKernel;

/// ---------------------------------------------------------------------------------------
//!
//! \brief This class provides functionality to operate upon a Triangular Mesh
//!
/// ---------------------------------------------------------------------------------------
class TriMesh {

protected:

    //! Name of the mesh
    std::string mName;

    //! Dimensionality of mesh (planar = 2D, surface = 3D)
    uint8_t mDim;

    //! The set of vertices
    std::vector<Vertex> mVertices;

    //! The set of faces
    std::vector<Face> mFaces;

    //! Other properties
    std::vector<Normal> mPointNormals, mFaceNormals;

    //! A bunch of fields
    std::unordered_map<std::string, std::vector<TypeFunction> > mFields;

    //! For each vertex, all neighboring vertices
    std::vector<std::vector<TypeIndex> > mVNeighbors;

    //! For each vertex, all neighboring faces
    std::vector<std::vector<TypeIndex> > mVAdjFaces;

    //! For each face, the three faces attached to its edges
    //!  (e.g., across_edge[3][2] is the number of the face
    //!   that's touching the edge opposite vertex 2 of face 3)
    std::vector<Offset3> mFAcrossEdge;

    //! Collection of edges on the boundary
    std::vector<Edge> bedges;

    //! geodesic distances

//#define PDIST
#ifndef PDIST
    std::vector<std::vector<TypeFunction>> mgeodesics;
#else
    std::vector<TypeFunction> mgeodesics;
#endif

protected:
    /// ---------------------------------------------------------------------------------------

    //static void compute_geodesics_cgal(const std::vector<Vertex> &mvertices, const std::vector<Face> &mfaces,
    //                                    std::vector<std::vector<TypeFunction>> &mdistances, bool verbose = false);

    static void compute_geodesics_fw(const std::vector<Vertex> &mvertices, const std::vector<Face> &mfaces,
                                     const DistanceKernel &dist, std::vector<std::vector<TypeFunction>> &mdistances, bool verbose = false);

    static Point3 Point2Bary(const Point3 &p, const Point3 &a, const Point3 &b, const Point3 &c);
    static Point2 Bary2Point(const Point3 &bary, const Point2 &a, const Point2 &b, const Point2 &c);

    /// ---------------------------------------------------------------------------------------
    //! linearize data into a single vector for interfacing with Python
    template <uint8_t D,  typename T, typename Tout>
    static std::vector<Tout> linearize(const std::vector<Vec<D,T> > &data, const uint8_t &vdim, size_t sz = 0) {

        static_assert(D == 2 || D == 3, "TriMesh::linearize() expected 2 or 3 dimensional vector");
        std::vector<Tout> ldata;

        if (sz == 0)
            sz = data.size();

        ldata.reserve(vdim*sz);
        for (size_t i = 0; i < sz; i++) {
        for (uint8_t d = 0; d < vdim; d++)
            ldata.push_back(data[i][d]);
        }
        return ldata;
    }

    //! linearize an array into data for interfacing with Python
    template <uint8_t D, typename T, typename Tin>
    static bool delinearize(const Tin *_, const size_t &n, const uint8_t &d,
                            std::vector<Vec<D,T> > &data) {

        static_assert((D == 2 || D == 3), "TriMesh::delinearize() expects 2 or 3 dimensional array");

        data.resize(n);
        if (d == 3) {
            for (size_t i=0; i<n; i++) {
                data[i] = Vec<D,T>(_[3*i], _[3*i+1], _[3*i+2]);
            }
        }
        else if(d == 2) {
            for (size_t i=0; i<n; i++) {
                data[i] = Vec<D,T>(_[2*i], _[2*i+1], 0.0);
            }
        }
        return true;
    }

    /// ---------------------------------------------------------------------------------------
    static void need_normals(const std::vector<Face> &faces, const std::vector<Vertex> &vertices,
                             std::vector<Normal> &fnormals, std::vector<Normal> &pnormals);

    static void need_pointareas(const std::vector<Face> &faces, const std::vector<Vertex> &vertices,
                                std::vector<TypeFunction> &areas);


    /// ---------------------------------------------------------------------------------------
    //! set dimensionalty of the mesh vertices
    bool set_dimensionality(uint8_t _);

    //! sort vertices spatially (using cgal)
    void sort_vertices(std::vector<Point_with_idx> &svertices) const;

    //! Compute connectivity
    void need_neighbors(bool verbose = false);
    void need_adjacentfaces(bool verbose = false);
    void need_across_edge(bool verbose = false);
    std::vector<TypeIndexI> need_boundary(bool verbose = false);

public:

    //! Constructors
    TriMesh() : mDim(0) {
        this->mName = "TriMesh";
    }

    TriMesh(float *_, int n, int d) {
        this->mName = "TriMesh";
        set_dimensionality(d);
        this->delinearize<3,TypeFunction,float>(_,n,d,mVertices);
    }
    /*TriMesh(double *_, int n, int d) {
        this->mName = "TriMesh";
        set_dimensionality(d);
        this->delinearize<3,double,TypeFunction>(_,n,d,mVertices);
    }*/

    TriMesh(const Polyhedron &surface_mesh);

    //! Destructor
    ~TriMesh() {}

    /// ---------------------------------------------------------------------------------------
    //! set vertices and faces
    bool set_faces(const TriMesh &mesh) {
        mFaces = mesh.mFaces;
        return true;
    }
    bool set_faces(uint32_t *_, int n, int d) {
        return this->delinearize<3,TypeIndex,uint32_t>(_,n,d,mFaces);
    }

    //! The number of vertices and faces
    size_t nvertices() const {  return mVertices.size();    }
    size_t nfaces() const {     return mFaces.size();       }

    //! Returns linearized vertices and faces
    std::vector<TypeFunction> get_vertices() const {    return linearize<3,TypeFunction,TypeFunction>(this->mVertices, this->mDim);    }
    std::vector<TypeIndexI> get_faces() const {         return linearize<3,TypeIndex,TypeIndexI>(this->mFaces, 3);          }

    //! Return linearized function
    std::vector<TypeFunction> get_field(const std::string &name) const {
        auto iter = mFields.find(name);
        return (iter != mFields.end()) ? iter->second : std::vector<TypeFunction> ();
    }

    //! Compute vertex normals
    const std::vector<TypeFunction> need_normals(bool verbose = false) {

        // Compute only if point areas are not available
        if (this->mPointNormals.size() != this->mVertices.size()) {

            if (verbose) {
                std::cout << "   > " << this->mName << "::need_normals()...";
                fflush(stdout);
            }

            TriMesh::need_normals(this->mFaces, this->mVertices, this->mFaceNormals, this->mPointNormals);

            if(verbose)
                std::cout << " Done!\n";
        }
        return linearize<3,TypeFunction,TypeFunction>(this->mPointNormals, 3);
    }

    //! Compute per-vertex point areas
    const std::vector<TypeFunction>& need_pointareas(bool verbose = false) {

        // Compute only if point areas are not available
        if (mFields.find("point_areas") == mFields.end()) {

            if (verbose) {
                std::cout << "   > " << this->mName << "::need_pointareas()...";
                fflush(stdout);
            }

            TriMesh::need_pointareas(this->mFaces, this->mVertices, this->mFields["point_areas"]);

            if(verbose)
                std::cout << " Done!\n";
        }
        return mFields["point_areas"];
    }


    //! Compute density
    const std::vector<TypeFunction>& kde(const int &type, const DensityKernel& k, const std::string &name, bool verbose = false) {
        return TriMesh::kde(type, k, name, std::vector<TypeIndexI>(), verbose);
    }
    const std::vector<TypeFunction>& kde(const int &type, const DensityKernel& k, const std::string &name, const std::vector<TypeIndexI> &ids, bool verbose = false);

    //! Compute curvature (using vtk)
    std::vector<TypeFunction> need_curvature(bool verbose = false);             // TriMesh_vtk.cpp

    //! parameterize the surface (using cgal)
    std::vector<TypeFunction> parameterize(bool verbose = false);
    std::vector<TypeFunction> parameterize_xy(bool verbose = false);

    //! project a set of points on the triangulation (using cgal)
    std::vector<TypeFunction> project_on_surface(const std::vector<TypeFunction> &points, bool verbose = false) const;

private:
    std::vector<TypeFunction> project_on_surface(const std::vector<Point3> &points, bool verbose = false) const;

public:
    //! compute the mesh as 2D Delaunay (using cgal)
    std::vector<TypeIndexI> delaunay(bool verbose = false);

    //! compute the distance of "this" mesh from the "other" mesh
    std::vector<TypeFunction> distance_to_other_mesh(const TriMesh &other, bool verbose=false) const;


    void geodesic() const;

#ifdef CPP_REMESHING
    //! Remesh
    void remesh(bool verbose = false);
#endif

    /// ---------------------------------------------------------------------------------------
    //! read/write off format
    static bool read_off(const std::string &fname, std::vector<Vertex> &get_vertices, std::vector<Face> &get_faces, bool verbose = false);
    static bool write_off(const std::string &fname, const std::vector<Vertex> &get_vertices, const std::vector<Face> &get_faces, const uint8_t &dim, bool verbose = false);

    bool read_off(const std::string &fname, bool verbose = false) {
        this->mDim = 3;
        return TriMesh::read_off(fname, mVertices, mFaces, verbose);
    }
    bool write_off(const std::string &fname, bool verbose = false) const {
        return TriMesh::write_off(fname, mVertices, mFaces, mDim, verbose);
    }

    //! write in binary format
    bool write_binary(const std::string &fname);

    //! write in vtp (paraview) format with or without periodic face
    bool write_vtp(const std::string &fname);
};

/// ---------------------------------------------------------------------------------------
//!
//! \brief This class provides functionality to operate upon a periodic Triangular Mesh
//!         periodic in x,y
//!
/// ---------------------------------------------------------------------------------------
class TriMeshPeriodic : public TriMesh {

private:

    //! Bounding box
    Vertex mBox0, mBox1;
    bool bbox_valid;

    //! in addition to the *actual* vertices and faces (i.e., inside the given domain)
    //! we need to store extra information

    //! a periodic vertex is stored as an offset (-1,0,1) with respect to the original vertex
    using periodicVertex = typename std::tuple<TypeIndex, int, int>;

    //! delaunay vertices
    std::vector<std::vector<periodicVertex> > mDelaunayFaces;

    //! the faces that go across the domain (contain original points)
    std::vector<Face> mPeriodicFaces;

    //! the faces that contain duplicate points (to replace periodic triangles)
    std::vector<Face> mTrimmedFaces;

    //! the periodic vertices that map to the original
    std::vector<periodicVertex> mDuplicateVertex_periodic;

    //! the vertices to support duplicate faces
    std::vector<Vertex> mDuplicateVerts;

    //! wrap vertices in xy (dim = 2) or in xyz (dim = 3)
    bool wrap_vertices(uint8_t dim = 2);

    //! trim the periodic faces using the map
    void lift_delaunay(bool verbose = false);


public:

    //! Constructor
    TriMeshPeriodic(float *_, int n, int d) :
        TriMesh(_,n,d) {
        this->mName = "TriMeshPeriodic";
    }
    /*TriMeshPeriodic(double *_, int n, int d) :
        TriMesh(_,n,d) {
        this->mName = "TriMeshPeriodic";
    }*/

    //! Destructor
    ~TriMeshPeriodic() {}

    /// ---------------------------------------------------------------------------------------
    //! set periodic box
    bool set_bbox(float *_, int n);

    //! set faces from a different triangulation
    void lift_delaunay(const TriMeshPeriodic &mesh) {
        this->mDelaunayFaces = mesh.mDelaunayFaces;
        lift_delaunay();
    }

    /// ---------------------------------------------------------------------------------------
    //! access the periodicty data

    std::vector<TypeIndexI> periodic_faces(bool combined = false) const {
        if (!combined)
            return linearize<3,TypeIndex,TypeIndexI>(mPeriodicFaces, 3);
        std::vector<TypeIndexI> f0 = linearize<3,TypeIndex,TypeIndexI>(mFaces, 3);
        std::vector<TypeIndexI> f1 = linearize<3,TypeIndex,TypeIndexI>(mPeriodicFaces, 3);
        f0.insert(f0.end(), f1.begin(), f1.end());
        return f0;
    }
    std::vector<TypeIndexI> trimmed_faces(bool combined = false) const {
        if (!combined)
            return linearize<3,TypeIndex,TypeIndexI>(mTrimmedFaces, 3);
        std::vector<TypeIndexI> f0 = linearize<3,TypeIndex,TypeIndexI>(mFaces, 3);
        std::vector<TypeIndexI> f1 = linearize<3,TypeIndex,TypeIndexI>(mTrimmedFaces, 3);
        f0.insert(f0.end(), f1.begin(), f1.end());
        return f0;
    }
    std::vector<TypeFunction> duplicated_vertices(bool combined = false) const {
        if (!combined)
            return linearize<3,TypeFunction,TypeFunction>(mDuplicateVerts, this->mDim);
        std::vector<TypeFunction> f0 = linearize<3,TypeFunction,TypeFunction>(mVertices, this->mDim);
        std::vector<TypeFunction> f1 = linearize<3,TypeFunction,TypeFunction>(mDuplicateVerts, this->mDim);
        f0.insert(f0.end(), f1.begin(), f1.end());
        return f0;
    }
    std::vector<TypeIndexI> duplicate_ids() const {
        const size_t ndups = mDuplicateVertex_periodic.size();
        std::vector<TypeIndexI> dids (ndups);
        for(size_t i = 0; i < ndups; i++)
            dids[i] = std::get<0>(mDuplicateVertex_periodic[i]);
        return dids;
    }

    /// ---------------------------------------------------------------------------------------
    const std::vector<TypeFunction>& need_pointareas(bool verbose = false) {

        // Compute only if point areas are not available
        if (mFields.find("point_areas") == mFields.end()) {

            if (verbose) {
                std::cout << "   > " << this->mName << "::need_pointareas()...";
                fflush(stdout);
            }

            std::vector<Face> faces = this->mFaces;
            faces.insert(faces.end(), this->mTrimmedFaces.begin(), this->mTrimmedFaces.end());

            std::vector<Vertex> vertices = this->mVertices;
            vertices.insert(vertices.end(), this->mDuplicateVerts.begin(), this->mDuplicateVerts.end());

            TriMesh::need_pointareas(faces, vertices, mFields["point_areas"]);
            mFields["point_areas"].resize(this->mVertices.size());

            if(verbose)
                std::cout << " Done!\n";
        }
        return mFields["point_areas"];
    }

    const std::vector<TypeFunction> need_normals(bool verbose = false) {

        // Compute only if point areas are not available
        if (this->mPointNormals.size() != this->mVertices.size()) {

            if (verbose) {
                std::cout << "   > " << this->mName << "::need_normals()...";
                fflush(stdout);
            }

            std::vector<Face> faces = this->mFaces;
            faces.insert(faces.end(), this->mTrimmedFaces.begin(), this->mTrimmedFaces.end());

            std::vector<Vertex> vertices = this->mVertices;
            vertices.insert(vertices.end(), this->mDuplicateVerts.begin(), this->mDuplicateVerts.end());

            TriMesh::need_normals(faces, vertices, this->mFaceNormals, this->mPointNormals);

            this->mFaceNormals.resize(this->mFaces.size());
            this->mPointNormals.resize(this->mVertices.size());

            if(verbose)
                std::cout << " Done!\n";
        }
        return linearize<3,TypeFunction,TypeFunction>(this->mPointNormals,3);
    }

    /// ---------------------------------------------------------------------------------------
    std::vector<TypeFunction> distance_to_other_mesh(const TriMeshPeriodic &other) const {
        return TriMesh::distance_to_other_mesh(other);
    }

    std::vector<TypeIndexI> delaunay(bool verbose = false);

    const std::vector<TypeFunction>& kde(const int &type, const DensityKernel& k, const std::string &name, bool verbose = false) {
        return TriMeshPeriodic::kde(type, k, name, std::vector<TypeIndexI>(), verbose);
    }
    const std::vector<TypeFunction>& kde(const int &type, const DensityKernel& k, const std::string &name, const std::vector<TypeIndexI> &ids, bool verbose = false);
};

/// ---------------------------------------------------------------------------------------
#endif    /* _TRIMESH_H_ */
