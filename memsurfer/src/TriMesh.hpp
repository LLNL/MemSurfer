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


//! A periodic vertex is stored as an offset (-1,0,1) with respect to the original vertex
struct PeriodicVertex {
    TypeIndex origVidx;
    int8_t offsetx, offsety;

    PeriodicVertex () :
        origVidx(0), offsetx(0), offsety(0) {}

    PeriodicVertex (const TypeIndex v, const int8_t ox, const int8_t oy)
        : origVidx(v), offsetx(ox), offsety(oy) {}

    inline bool operator==(const PeriodicVertex &p) const {
        return origVidx == p.origVidx && offsetx == p.offsetx && offsety == p.offsety;
    }
    inline bool is_original() const {
        return offsetx == 0 and offsety == 0;
    }
};

inline
std::ostream&
operator<<(std::ostream &os, const PeriodicVertex &p) {
    return os << "("<<p.origVidx<<" : ["<<int(p.offsetx)<<","<<int(p.offsety)<<"])";
}

/// ---------------------------------------------------------------------------------------
//!
//! \brief This class provides functionality to operate upon a Triangular Mesh
//!
/// ---------------------------------------------------------------------------------------
class TriMesh {

private:
    //! Dimensionality of mesh (planar = 2D, surface = 3D)
    uint8_t mDim;

    //! Whether this mesh is periodic (in xy) or not
    bool mPeriodic;

    //! -----------------------------------------------------------------------------------
    //! Geometric elements
    //! -----------------------------------------------------------------------------------

    //! The set of vertices
    std::vector<Vertex> mVertices;

    //! The set of faces
    std::vector<Face> mFaces;

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

    //! -----------------------------------------------------------------------------------
    //! Geometric elements to handle periodicity
    //! -----------------------------------------------------------------------------------

    //! Bounding box
    bool bbox_valid;
    Vertex mBox0, mBox1;

    //! In addition to the *actual* vertices and faces (i.e., inside the given domain),
    //! we need to store extra information

    //! Delaunay Faces
    std::vector<std::vector<PeriodicVertex> > mDelaunayFaces;

    //! The faces that go across the domain (contain original points)
    std::vector<Face> mPeriodicFaces;

    //! The faces that contain duplicate points (to replace periodic triangles)
    std::vector<Face> mTrimmedFaces;

    //! The duplicate vertices that map to the original
    std::vector<PeriodicVertex> mDuplicateVertex_periodic;

    //! The vertices to support duplicate faces
    std::vector<Vertex> mDuplicateVerts;

    //! -----------------------------------------------------------------------------------

    //! A bunch of fields defined on the mesh
    std::unordered_map<std::string, std::vector<TypeFunction>> mFields;

    //! Normals
    std::vector<Normal> mPointNormals, mFaceNormals;

    //! Geodesic distances
//#define PDIST
#ifndef PDIST
    std::vector<std::vector<TypeFunction>> mgeodesics;
#else
    std::vector<TypeFunction> mgeodesics;
#endif

    //! -----------------------------------------------------------------------------------
    //! Static methods to compute properties of interest
    //! -----------------------------------------------------------------------------------

    static Point3 Point2Bary(const Point3 &p, const Point3 &a, const Point3 &b, const Point3 &c);
    static Point2 Bary2Point(const Point3 &bary, const Point2 &a, const Point2 &b, const Point2 &c);

    //! -----------------------------------------------------------------------------------
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

    //! -----------------------------------------------------------------------------------
    //! compute normals for a set of vertices and faces
    static void need_normals(const std::vector<Face> &faces, const std::vector<Vertex> &vertices,
                             std::vector<Normal> &fnormals, std::vector<Normal> &pnormals);

    //! compute point areas for a set of vertices
    static void need_pointareas(const std::vector<Face> &faces, const std::vector<Vertex> &vertices,
                                std::vector<TypeFunction> &areas);


    //! compute graph geodesics
    static void compute_geodesics_fw(const std::vector<Vertex> &mvertices, const std::vector<Face> &mfaces,
                                     const DistanceKernel &dist, std::vector<std::vector<TypeFunction>> &mdistances,
                                     bool verbose = false);

    //static void compute_geodesics_cgal(const std::vector<Vertex> &mvertices, const std::vector<Face> &mfaces,
    //                                    std::vector<Field> &mdistances, bool verbose = false);

    //! -----------------------------------------------------------------------------------
    //! nonstatic private methods
    //! -----------------------------------------------------------------------------------
    std::vector<TypeFunction> get_pairwise_distances(uint8_t D) const;

    //! set dimensionalty of the mesh vertices
    bool set_dimensionality(uint8_t _);

    //! sort vertices spatially (using cgal)
    std::vector<Point_with_idx> sort_vertices() const;

    //! wrap vertices in xy (dim = 2) or in xyz (dim = 3)
    bool wrap_vertices(uint8_t dim = 2);

    //! trim the periodic faces generated using periodic Delaunay
    void trim_periodicDelaunay(bool verbose = false);

    //! project a set of points on the surface (using cgal)
    std::vector<TypeFunction> project_on_surface(const std::vector<Point3> &points, bool verbose = false) const;

public:

    //! -----------------------------------------------------------------------------------
    //! API
    //! -----------------------------------------------------------------------------------

    //! Constructors
    TriMesh() : mDim(0) {
        this->mPeriodic = false;
    }
    TriMesh(float *_, int n, int d) {
        this->mPeriodic = false;
        this->set_dimensionality(d);
        this->set_vertices(_,n,d);
    }
    TriMesh(const Polyhedron &surface_mesh);

    //! Destructor
    ~TriMesh() {}

    std::string tag() const {
        return this->mPeriodic ? "TriMeshPeriodic":"TriMesh";
    }

    //! -----------------------------------------------------------------------------------
    //! Set and get vertices and faces
    bool set_vertices(float *_, int n, int d) {
        return this->delinearize<3,TypeFunction,float>(_,n,d,mVertices);
    }
    bool set_faces(uint32_t *_, int n, int d) {
        return this->delinearize<3,TypeIndex,uint32_t>(_,n,d,mFaces);
    }
    bool set_faces(const TriMesh &mesh) {
        mFaces = mesh.mFaces;
        return true;
    }

    //! The number of vertices and faces
    size_t nvertices() const {  return mVertices.size();    }
    size_t nfaces() const {     return mFaces.size();       }

    //! Returns linearized vertices and faces
    std::vector<TypeFunction> get_vertices() const {    return linearize<3,TypeFunction,TypeFunction>(this->mVertices, this->mDim);    }
    std::vector<TypeIndexI> get_faces() const {         return linearize<3,TypeIndex,TypeIndexI>(this->mFaces, 3);          }

    std::vector<TypeIndexI> need_neighbors();

    //! -----------------------------------------------------------------------------------
    //! Set a field with a name
    bool set_field(std::string key, float *_, int n, int d) {

        std::vector<TypeFunction> &v = mFields[key];
        v.resize(n);
        for(int i = 0; i < n; i++){
          v[i] = _[i];
        }
        return true;
    }
    bool set_fields(const TriMesh &mesh, std::string key) {
        for (auto iter = mesh.mFields.begin(); iter != mesh.mFields.end(); iter++) {
            if (iter->first.find(key) == 0) {
                this->mFields[iter->first] = iter->second;
            }
        }
        return true;
    }

    //! Return linearized field
    std::vector<TypeFunction> get_field(const std::string &name) const {
        auto iter = mFields.find(name);
        return (iter != mFields.end()) ? iter->second : std::vector<TypeFunction> ();
    }

    //! -----------------------------------------------------------------------------------
    //! set periodic box
    bool set_periodic() {
        this->mPeriodic = true;
        return true;
    }
    bool set_bbox(float *_, int n);

    //! set faces from a different triangulation and then trim
    void copy_periodicDelaunay(const TriMesh &mesh) {
        this->mDelaunayFaces = mesh.mDelaunayFaces;
        trim_periodicDelaunay();
    }

    //! -----------------------------------------------------------------------------------
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
            dids[i] = mDuplicateVertex_periodic[i].origVidx;
        return dids;
    }

    //! -----------------------------------------------------------------------------------
    //! Compute mesh properties
    //! -----------------------------------------------------------------------------------

    //! Compute vertex normals
    const std::vector<TypeFunction> need_normals(bool verbose = false) {

        // Compute only if point areas are not available
        if (this->mPointNormals.size() != this->mVertices.size()) {

            if (verbose) {
                std::cout << "   > " << tag() << "::need_normals()...";
                fflush(stdout);
            }

            if (!this->mPeriodic) {
                TriMesh::need_normals(this->mFaces, this->mVertices, this->mFaceNormals, this->mPointNormals);
            }
            else {
                std::vector<Face> faces = this->mFaces;
                faces.insert(faces.end(), this->mTrimmedFaces.begin(), this->mTrimmedFaces.end());

                std::vector<Vertex> vertices = this->mVertices;
                vertices.insert(vertices.end(), this->mDuplicateVerts.begin(), this->mDuplicateVerts.end());

                TriMesh::need_normals(faces, vertices, this->mFaceNormals, this->mPointNormals);

                this->mFaceNormals.resize(this->mFaces.size());
                this->mPointNormals.resize(this->mVertices.size());
            }

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
                std::cout << "   > " << tag() << "::need_pointareas()...";
                fflush(stdout);
            }

            if (!this->mPeriodic) {
                TriMesh::need_pointareas(this->mFaces, this->mVertices, this->mFields["point_areas"]);
            }
            else {
                std::vector<Face> faces = this->mFaces;
                faces.insert(faces.end(), this->mTrimmedFaces.begin(), this->mTrimmedFaces.end());

                std::vector<Vertex> vertices = this->mVertices;
                vertices.insert(vertices.end(), this->mDuplicateVerts.begin(), this->mDuplicateVerts.end());

                TriMesh::need_pointareas(faces, vertices, mFields["point_areas"]);
                mFields["point_areas"].resize(this->mVertices.size());
            }

            if(verbose)
                std::cout << " Done!\n";
        }
        return mFields["point_areas"];
    }

    //! compute the mesh as 2D Delaunay (using cgal)
    std::vector<TypeIndexI> delaunay(bool verbose = false);
    std::vector<TypeIndexI> periodicDelaunay(bool verbose = false);

    //! compute the distance of "this" mesh from the "other" mesh
    std::vector<TypeFunction> distance_to_other_mesh(const TriMesh &other, bool verbose=false) const;

    //! parameterize the surface (using cgal)
    std::vector<TypeFunction> parameterize(bool verbose = false);
    std::vector<TypeFunction> parameterize_xy(bool verbose = false);

    //! project a set of points on the triangulation (using cgal)
    std::vector<TypeFunction> project_on_surface(const std::vector<TypeFunction> &points, bool verbose = false) const;

    //! -----------------------------------------------------------------------------------
    //! Compute density
    //! -----------------------------------------------------------------------------------
    const std::vector<TypeFunction>&
        kde(const std::string &name, const int type, const bool get_counts,
            const DensityKernel& dens_kern, const DistanceKernel& dist_kern,
            const bool verbose = false) {

        return TriMesh::kde(name, type, get_counts, dens_kern, dist_kern,
                            std::vector<TypeIndexI>(), verbose);
    }
    const std::vector<TypeFunction>&
        kde(const std::string &name, const int type, const bool get_counts,
            const DensityKernel& dens_kern, const DistanceKernel& dist_kern,
            const std::vector<TypeIndexI> &ids, const bool verbose = false);

    //! -----------------------------------------------------------------------------------

public:

#ifdef CGAL_GEODESIC
    void geodesic() const;
#endif
#ifdef CPP_REMESHING
    //! Remesh
    void remesh(bool verbose = false);
#endif

    /// ---------------------------------------------------------------------------------------
    //! read/write off format
    static bool read_off(const std::string &fname, std::vector<Vertex> &vertices, std::vector<Face> &faces, bool verbose = false);
    static bool write_off(const std::string &fname, const std::vector<Vertex> &vertices, const std::vector<Face> &faces, const uint8_t &dim, bool verbose = false);

    bool read_off(const std::string &fname, bool verbose = false) {
        this->mDim = 3;
        return TriMesh::read_off(fname, mVertices, mFaces, verbose);
    }
    bool write_off(const std::string &fname, bool verbose = false) const {
        return TriMesh::write_off(fname, mVertices, mFaces, mDim, verbose);
    }

    //! write in binary format
    static bool write_binary(const std::string &fname,
                             const std::vector<Vertex> &vertices, const std::vector<Face> &faces,
                             const std::vector<std::string> &field_names,
                             const std::vector<std::vector<TypeFunction>*> &fields,
                             bool verbose = false);

    bool write_binary(const std::string &fname, const std::string &filter_fields="");
};

/// ---------------------------------------------------------------------------------------
#endif    /* _TRIMESH_H_ */
