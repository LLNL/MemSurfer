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

#include <map>
#include <sstream>
#include <stdexcept>
#include <chrono>

#include "TriMesh.hpp"
#include "DensityKernels.hpp"
#include "DistanceKernels.hpp"

#ifdef PDIST
size_t square_to_condensed(const size_t &i, const size_t &j, const size_t &n) {
    if (i == j) {
        std::cerr << " square_to_condensed does not handle i=j!\n";
        exit(1);
    }
    if (i < j) {
        return square_to_condensed(j,i,n);
    }

    //return (j * (n - (j+1)>>1) + i-j-1);
    return (n*j - j*(j+1)/2 + i-j-1);
}
#endif


#ifdef CGAL_GEODESIC
#include "Types.hpp"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Random.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_shortest_path.h>
#include <boost/lexical_cast.hpp>

//typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
//typedef CGAL::Surface_mesh<Kernel::Point_3> SurfaceMesh;

typedef CGAL::Surface_mesh<Kernel::Point_3>                     SurfaceMesh;

typedef boost::graph_traits<SurfaceMesh>::vertex_descriptor     vertex_descriptor;
typedef boost::graph_traits<SurfaceMesh>::face_descriptor       face_descriptor;
typedef boost::graph_traits<SurfaceMesh>::edge_descriptor       edge_descriptor;
typedef boost::graph_traits<SurfaceMesh>::halfedge_descriptor   halfedge_descriptor;


typedef boost::graph_traits<SurfaceMesh> Graph_traits;
typedef Graph_traits::vertex_iterator vertex_iterator;
typedef Graph_traits::face_iterator face_iterator;

typedef CGAL::Surface_mesh_shortest_path_traits<Kernel, SurfaceMesh> Traits;
typedef CGAL::Surface_mesh_shortest_path<Traits> Surface_mesh_shortest_path;

#include <CGAL/boost/graph/dijkstra_shortest_paths.h>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <boost/graph/floyd_warshall_shortest.hpp>
#include <boost/foreach.hpp>

/*
void TriMesh::compute_geodesics_cgal(const std::vector<Vertex> &mvertices, const std::vector<Face> &mfaces,
#ifdef PDIST
                                     std::vector<TypeFunction> &distances,
#else
                                     std::vector<std::vector<TypeFunction>> &distances,
#endif
                                     bool verbose) {

    const size_t nverts = mvertices.size();
    const size_t nfaces = mfaces.size();

    if (verbose) {
        std::cout << " TriMesh::compute_geodesics_cgal() for " << nverts << " vertices and "<< nfaces << " faces!\n";
    }

    using namespace std::chrono;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    // create a surface mesh
    SurfaceMesh smesh;

    std::vector<vertex_descriptor> smvd (nverts);
    for(size_t i = 0; i < nverts; i++) {
        const Vertex &mv = mvertices[i];
        smvd[i] = smesh.add_vertex(Kernel::Point_3(mv[0], mv[1], mv[2]));
    }

    // somehow, CGAL is not able to locate the face based on the vertex descriptor
    // so, to get things to work, let's create an explicit map
    std::map<size_t, std::pair<size_t, uint8_t>> vmap;

    for(size_t i = 0; i < nfaces; i++) {

        const Face &mf = mfaces[i];
        smesh.add_face(smvd[mf[0]], smvd[mf[1]], smvd[mf[2]]);

        for(uint8_t d = 0; d < 3; d++) {
            if (vmap.find(mf[d]) == vmap.end()) {
                vmap[mf[d]] = std::make_pair(i, d);
            }
        }
    }

    if (vmap.size() != nverts) {
        std::cout <<" TriMesh::kde_2m(): could not map all vertices!\n";
        //exit(1);
    }

    // construct a shortest path query object
    Surface_mesh_shortest_path shortest_paths(smesh);

    if (verbose) {
        std::cout << " mesh has " << smesh.num_vertices() << " vertices and " << smesh.number_of_faces() << " faces "
                  << smesh.num_edges() << ", " << smesh.num_halfedges() << "!\n";
    }

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    //std::cout << " step 1 took " << duration_cast<duration<double>>(t2-t1).count() << " seconds.\n";

#ifdef PDIST
    size_t npairs = nverts*(nverts-1)/2;
    distances.resize(npairs, FLT_MAX);
#else
    distances.resize(nverts);
    for(size_t i = 0; i < nverts; i++)
        distances[i].resize(nverts, FLT_MAX);
#endif

    // get distance of the source point from each vertex in the mesh
    for(auto viter = vmap.begin(); viter != vmap.end(); viter++) {

        const size_t vidx = viter->first;           // looking at this vertex
        const size_t fidx = viter->second.first;    //  in this face
        const uint8_t cidx = viter->second.second;  //  at this corner

        t1 = high_resolution_clock::now();

        // represent the vertex as (face, corner)
        face_iterator face_it = faces(smesh).first;
        std::advance(face_it, fidx);

        Traits::Barycentric_coordinates face_location = {{0.0,0.0,0.0}};
        face_location[cidx] = 1.0;

        // add the source point
        shortest_paths.add_source_point(*face_it, face_location);

        // for this vertex, compute the geodesic distance to all other vertices
        size_t v2idx = vidx+1;
        vertex_iterator vert_it = smesh.vertices_begin();
        std::advance(vert_it, v2idx);

        t2 = high_resolution_clock::now();
        //std::cout << vidx << ": " << fidx << ": "<< int(cidx) << " :: " << duration_cast<duration<double>>(t2-t1).count() << " seconds." << std::endl;

        t1 = high_resolution_clock::now();
        size_t cnt = 0;
        for ( ; vert_it != smesh.vertices_end(); ++vert_it, ++v2idx) {

            auto r = shortest_paths.shortest_distance_to_source_points(*vert_it);

#ifndef PDIST
            distances[vidx][v2idx] = r.first;
            distances[v2idx][vidx] = r.first;
#else
            distances[square_to_condensed(vidx, v2idx, nverts)] = r.first;
#endif
        }
        shortest_paths.remove_all_source_points();

        t2 = high_resolution_clock::now();
    }

    if (verbose) {
        std::cout << " TriMesh::compute_geodesics_cgal() done! computed " << distances.size() << " values!\n";
    }
}
*/
#endif

//! Floyd-Warshall algorithm
void TriMesh::compute_geodesics_fw(const std::vector<Vertex> &mvertices, const std::vector<Face> &mfaces,
                                   const DistanceKernel &dist,
#ifdef PDIST
                                    std::vector<TypeFunction> &distances,
#else
                                    std::vector<std::vector<TypeFunction>> &distances,
#endif
                                   bool verbose) {

    using namespace std::chrono;
    const size_t nverts = mvertices.size();
    const size_t nfaces = mfaces.size();
    const TypeFunction threshold = 10000.0;

#ifdef PDIST
    const size_t npairs = nverts*(nverts-1)/2;
    distances.resize(npairs, FLT_MAX);
#else
    distances.resize(nverts);
    for(size_t i = 0; i < nverts; i++)
        distances[i].resize(nverts, FLT_MAX);
#endif

    if (verbose) {
        std::cout << " TriMesh::compute_geodesics_fw(" << nverts << " verts, "<< nfaces << " faces)!\n";
    }

    high_resolution_clock::time_point t1 = high_resolution_clock::now();

#ifndef PDIST
    // initialize the distances for self
    for(size_t i = 0; i < nverts; i++) {
        distances[i][i] = 0;
    }
#endif

    // initialize the distances using immediate neighbors
    for(auto fiter = mfaces.begin(); fiter != mfaces.end(); fiter++) {

        const Face &mf = *fiter;

        for(uint8_t d = 0; d < 3; d++) {

            const size_t &a = mf[d];
            const size_t &b = mf[d==2?0:d+1];

            // distance kernels return square distance
            TypeFunction val = dist(mvertices[a][0], mvertices[a][1], mvertices[a][2],
                                    mvertices[b][0], mvertices[b][1], mvertices[b][2]);

            // FW algorithm needs to "add" distances
            // so, we should take sqrt (convert to actual distance)
            val = std::sqrt(val);

            /*
#ifdef PDIST
            size_t ab = square_to_condensed(a,b, nverts);
            TypeFunction &val = distances[ab];
#else
            TypeFunction &val = distances[a][b];
#endif

            // i do not have a separate structure for edges
            // so use this condition to avoid duplicate work
            if (val < threshold)
                continue;

            val = len(mvertices[a] - mvertices[b]);
            */

            distances[a][b] = val;
#ifndef PDIST
            // not needed, but trying to figure out the issue
            distances[b][a] = distances[a][b];
#endif
            //std::cout << " : " << a << "--"<< b<< " = " << val << " : " << val1 << " -- "
            //          << mvertices[a] << " : "<< mvertices[b] << std::endl;
        }
    }
    //std::cout << "  did " << mfaces.size() << std::endl;

#ifdef PDIST
    size_t ik, kj, ij;
    for (size_t k = 0; k < nverts; k++) {
    for (size_t i = 0; i < nverts; i++)  {

        if (i==k)   continue;
        ik = square_to_condensed(i, k, nverts);
        for (size_t j = 0; j < nverts; j++) {

            if (i==j || j==k)   continue;
            kj = square_to_condensed(k, j, nverts);
            ij = square_to_condensed(i, j, nverts);

            if (distances[ik] + distances[kj] < distances[ij]){
                distances[ij] = distances[ik] + distances[kj];
            }
    }}}
#else
    for (size_t k = 0; k < nverts; k++) {
    for (size_t i = 0; i < nverts; i++) {
    for (size_t j = 0; j < nverts; j++) {

        if (distances[i][k] + distances[k][j] < distances[i][j]) {
            distances[i][j] = distances[i][k] + distances[k][j];
            distances[j][i] = distances[i][j];
        }
    }}}
#endif

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    if (verbose) {
        std::cout << " TriMesh::compute_geodesics_fw() done! computed " << distances.size() << " values!\n";
        std::cout << "\t took " << duration_cast<duration<double>>(t2-t1).count() << " seconds." << std::endl;
    }

#ifndef PDIST
    for(size_t i = 0; i < nverts; i++) {
    for(size_t j = 0; j < nverts; j++) {

        if (distances[i][j] > threshold) {
            std::cout << " invalid: " << i << ", " << j << " = " << distances[i][j] << " : "
                      << mvertices[i] << " : "<< mvertices[j] << std::endl;
            //exit(1);
        }
    }}
#else
    for(size_t i = 0; i < npairs; i++) {

        //std::cout << " : " << distances[i] << std::endl;
        if (distances[i] > threshold) {
            std::cout << " invalid: " << i << " = " << distances[i] << "\n";
            exit(1);
        }
    }
#endif
}

/// -----------------------------------------------------------------------------
//! density estimation (core function)
/// -----------------------------------------------------------------------------

/*void normalize(std::vector<TypeFunction> &_, const size_t n) {
    const TypeFunction sm = std::accumulate(_.begin(), _.end(), TypeFunction(0));
    std::transform(_.begin(), _.end(), _.begin(),
               std::bind(std::multiplies<TypeFunction>(), std::placeholders::_1, TypeFunction(n)/sm));
}*/

void
kde_2d(const std::vector<Vertex> &vertices, const std::vector<TypeIndexI> &ids,
       const DensityKernel& k, const DistanceKernel &dist,
       std::vector<TypeFunction> &density) {

    const size_t nids = ids.size();
    const size_t nverts = vertices.size();
    density.resize(nverts, 0);

    if (nids == 0) {  // compute for all ids
        for (TypeIndex j=0; j<nverts; j++) {
        for (TypeIndex i=0; i<nverts; i++) {
            density[j] += k(dist(vertices[i][0], vertices[i][1],
                                 vertices[j][0], vertices[j][1]));
        }}
    }
    else {              // compute for selected ids
        for (TypeIndex j=0; j<nverts; j++) {
        for (TypeIndex i=0; i<nids;   i++) {
            density[j] += k(dist(vertices[ids[i]][0], vertices[ids[i]][1],
                                 vertices[j][0],      vertices[j][1]));
        }}
    }
}

void
kde_3d(const std::vector<Vertex> &vertices, const std::vector<TypeIndexI> &ids,
       const DensityKernel& k, const DistanceKernel &dist,
       std::vector<TypeFunction> &density) {

    const size_t nids = ids.size();
    const size_t nverts = vertices.size();
    density.resize(nverts, 0);

    if (nids == 0) {  // compute for all ids
        for (TypeIndex j=0; j<nverts; j++) {
        for (TypeIndex i=0; i<nverts; i++) {
            density[j] += k(dist(vertices[i][0], vertices[i][1], vertices[i][2],
                                 vertices[j][0], vertices[j][1], vertices[j][2]));
        }}
    }
    else {              // compute for selected ids
        for (TypeIndex j=0; j<nverts; j++) {
        for (TypeIndex i=0; i<nids;   i++) {
            density[j] += k(dist(vertices[ids[i]][0], vertices[ids[i]][1], vertices[ids[i]][2],
                                 vertices[j][0],      vertices[j][1],      vertices[j][2]));
        }}
    }
}

void kde_2m(const size_t &nverts, const std::vector<TypeIndexI> &ids,
            const DensityKernel& k,
#ifdef PDIST
            const std::vector<TypeFunction> &distances,
#else
            const std::vector<std::vector<TypeFunction>> &distances,
#endif
            std::vector<TypeFunction> &density) {


    const size_t nids = ids.size();
    density.resize(nverts, 0);

    static TypeFunction tmp;

    if (ids.empty()) {  // compute for all ids
        for (TypeIndex j=0; j<nverts; j++) {
        for (TypeIndex i=0; i<nverts; i++) {

#ifdef PDIST
            tmp = distances[square_to_condensed(i,j,nverts)];
#else
            tmp = distances[i][j];
#endif
            // geodesic distances are stored as actual distances,
            // but density kernel requires squared distance
            tmp = k(tmp*tmp);
            density[j] += tmp;
        }}
    }
    else {              // compute for selected ids
        for (TypeIndex j=0; j<nverts; j++) {
        for (TypeIndex i=0; i<nids;   i++) {

#ifdef PDIST
            tmp = distances[square_to_condensed(ids[i],j,nverts)];
#else
            tmp = distances[ids[i]][j];
#endif
            tmp = k(tmp);
            density[j] += tmp;
        }}
    }
}

/// -----------------------------------------------------------------------------
//! density estimation for nonperiodic mesh
/// -----------------------------------------------------------------------------

const std::vector<TypeFunction>&
    TriMesh::kde(const std::string &name, const int type, const bool get_counts,
                 const DensityKernel& dens, const DistanceKernel& dist,
                 const std::vector<TypeIndexI> &ids, const bool verbose) {

    if (type < 1 || type > 3) {
        std::ostringstream errMsg;
        errMsg << "   > " << this->tag() << "::kde(" << type << ">): invalid density type (should be 1, 2, or 3)!\n";
        throw std::invalid_argument(errMsg.str());
    }

    // this name already exists in the map!
    if (mFields.find(name) != mFields.end()) {
        if (verbose){
            std::cout << " " << this->tag() << "::kde("<<name<<") got an already used name. Returning existing field!\n";
        }
        return mFields.at(name);
    }

    if (verbose){
        std::cout << "   > " << this->tag() << "::kde(<" << ids.size() << ">, "<<name<<")...";
        fflush(stdout);
    }

    // we will be using this!
    std::vector<TypeFunction> &density = mFields[name];

    // now, compute the appropriate density!
    if (type == 2) {      kde_2d(mVertices, ids, dens, dist, density); }
    else if (type == 3){  kde_3d(mVertices, ids, dens, dist, density); }

    // geodesic density!
    else {
        // initialize the geodesic graph!
        if (mgeodesics.empty()) {
          if (this->mPeriodic) {
            std::vector<Face> mfaces = this->mFaces;
            mfaces.insert(mfaces.end(), this->mPeriodicFaces.begin(), this->mPeriodicFaces.end());
            compute_geodesics_fw(this->mVertices, mfaces, dist, this->mgeodesics, verbose);
          }
          else {
            compute_geodesics_fw(this->mVertices, this->mFaces, dist, this->mgeodesics, verbose);
          }
        }

        kde_2m(this->mVertices.size(), ids, dens, this->mgeodesics, mFields[name]);
    }

    // -------------------------------------------------------------------------
    // now, normalize the density
    // -------------------------------------------------------------------------
    // need two types of normalizations!
    //    1. we need to divide by nverts
    //        because we added nverts gaussians (one for each vertex)
    //        this is needed by kde (to turn this into a pdf)
    //        if this function was to be sampled dense enough, it would sum to 1
    //    2. if get_counts = True
    //        we want to turn this pdf into number of particles
    //        so we want to make sure the sum is equal to number of ids
    // -------------------------------------------------------------------------

    const size_t ng = this->mVertices.size();           // num of gaussians

    // first normalization
    if (1) {

      const TypeFunction norm = 1.0 / TypeFunction(ng);
      std::transform(density.begin(), density.end(), density.begin(),
                     std::bind(std::multiplies<TypeFunction>(), std::placeholders::_1, norm));
    }

    // second normalization
    if (get_counts) {

      const size_t np = (ids.empty()) ? ng : ids.size();  // num of points counted
      const TypeFunction dsm = std::accumulate(density.begin(), density.end(), TypeFunction(0));
      const TypeFunction norm = (TypeFunction(np) / dsm);

      // now, do the normalization
      std::transform(density.begin(), density.end(), density.begin(),
                     std::bind(std::multiplies<TypeFunction>(), std::placeholders::_1, norm));
    }

    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    if(verbose){
        printf(" Done!\n");
    }
    return mFields.at(name);
}

/// -----------------------------------------------------------------------------
/// -----------------------------------------------------------------------------
