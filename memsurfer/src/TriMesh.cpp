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
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <unordered_set>

#include "TriMesh.hpp"

//! -----------------------------------------------------------------------------
//! static functions
//! -----------------------------------------------------------------------------

Point3 TriMesh::Point2Bary(const Point3 &p, const Point3 &a, const Point3 &b, const Point3 &c) {

    Vector3 v0 = b - a, v1 = c - a, v2 = p - a;
    double d00 = CGAL::scalar_product(v0, v0);
    double d01 = CGAL::scalar_product(v0, v1);
    double d11 = CGAL::scalar_product(v1, v1);
    double d20 = CGAL::scalar_product(v2, v0);
    double d21 = CGAL::scalar_product(v2, v1);
    double denom = d00 * d11 - d01 * d01;

    double v = (d11 * d20 - d01 * d21) / denom;
    double w = (d00 * d21 - d01 * d20) / denom;
    double u = 1.0f - v - w;

    return Point3(u,v,w);
}

Point2 TriMesh::Bary2Point(const Point3 &bary, const Point2 &a, const Point2 &b, const Point2 &c) {

    double x = bary[0]*a[0] + bary[1]*b[0] + bary[2]*c[0];
    double y = bary[0]*a[1] + bary[1]*b[1] + bary[2]*c[1];
    return Point2(x, y);
}

//! -----------------------------------------------------------------------------
//! mesh manipulation and queries
//! -----------------------------------------------------------------------------

bool TriMesh::set_dimensionality(uint8_t _) {

    mDim = _;
    if (mDim != 2 && mDim != 3) {
        std::ostringstream errMsg;
        errMsg << " " << this->tag() << "::set_dimensionality(): Invalid dimensionality of vertices! can only be 2 or 3, but got " << int(mDim) << std::endl;
        throw std::invalid_argument(errMsg.str());
    }
    return true;
}

bool TriMesh::set_bbox(float *_, int n) {

    if (n == 2 && this->mDim == 2) {
        this->mBox0 = Vertex(0, 0, 0);
        this->mBox1 = Vertex(_[0], _[1], 0);
    }
    else if (n == 3 && this->mDim == 3) {
        this->mBox0 = Vertex(0, 0, 0);
        this->mBox1 = Vertex(_[0], _[1], _[2]);
    }
    else if (n == 4 && this->mDim == 2) {
        this->mBox0 = Vertex(_[0], _[1], 0);
        this->mBox1 = Vertex(_[2], _[3], 0);
    }
    else if (n == 6 && this->mDim == 3) {
        this->mBox0 = Vertex(_[0], _[1], _[2]);
        this->mBox1 = Vertex(_[3], _[4], _[5]);
    }
    else {
        std::ostringstream errMsg;
        errMsg << " " << this->tag() << "::set_bbox(): Invalid bounding box! got " << n << " values for " << int(this->mDim) << "D" << std::endl;
        throw std::invalid_argument(errMsg.str());
    }
    bbox_valid = true;
    return true;
}

//! -----------------------------------------------------------------------------
//! compute vertex normals
//! -----------------------------------------------------------------------------

// static
void TriMesh::need_normals(const std::vector<Face> &faces, const std::vector<Vertex> &vertices,
                           std::vector<Normal> &fnormals, std::vector<Normal> &pnormals) {

    const size_t nv = vertices.size();
    const size_t nf = faces.size();

    fnormals.resize(nf, Normal(0,0,0));
    pnormals.resize(nv, Normal(0,0,0));

#pragma omp parallel for
    for (size_t i = 0; i < nf; i++) {

        const Face &f = faces[i];
        const Vertex &p0 = vertices[f[0]];
        const Vertex &p1 = vertices[f[1]];
        const Vertex &p2 = vertices[f[2]];

        Normal a = p0-p1, b = p1-p2, c = p2-p0;

        TypeFunction l2a = len2(a), l2b = len2(b), l2c = len2(c);
        if (!l2a || !l2b || !l2c)
            continue;

        fnormals[i] = a CROSS b;

        pnormals[f[0]] += fnormals[i] * TypeFunction(1.0 / (l2a * l2c));
        pnormals[f[1]] += fnormals[i] * TypeFunction(1.0 / (l2b * l2a));
        pnormals[f[2]] += fnormals[i] * TypeFunction(1.0 / (l2c * l2b));
    }

    // Make them all unit-length
#pragma omp parallel for
    for (size_t i = 0; i < nv; i++)
        normalize(pnormals[i]);
}

//! -----------------------------------------------------------------------------
//! compute per-vertex point areas
//! -----------------------------------------------------------------------------

// static
void TriMesh::need_pointareas(const std::vector<Face> &faces, const std::vector<Vertex> &vertices,
                              std::vector<TypeFunction> &areas) {

    size_t nv = vertices.size();
    size_t nf = faces.size();

    areas.resize(nv, 0.0);
    std::vector<Vertex> cornerareas(nf);

#pragma omp parallel for
    for (size_t i = 0; i < nf; i++) {

        // Edges
        Vertex e[3] = { vertices[faces[i][2]] - vertices[faces[i][1]],
                        vertices[faces[i][0]] - vertices[faces[i][2]],
                        vertices[faces[i][1]] - vertices[faces[i][0]] };

        // Compute corner weights
        TypeFunction area = 0.5f * len(e[0] CROSS e[1]);
        TypeFunction l2[3] = { len2(e[0]), len2(e[1]), len2(e[2]) };
        TypeFunction ew[3] = { l2[0] * (l2[1] + l2[2] - l2[0]),
                               l2[1] * (l2[2] + l2[0] - l2[1]),
                               l2[2] * (l2[0] + l2[1] - l2[2]) };

        if (ew[0] <= 0.0f) {
            cornerareas[i][1] = -0.25f * l2[2] * area / (e[0] DOT e[2]);
            cornerareas[i][2] = -0.25f * l2[1] * area / (e[0] DOT e[1]);
            cornerareas[i][0] = area - cornerareas[i][1] - cornerareas[i][2];
        }
        else if (ew[1] <= 0.0f) {
            cornerareas[i][2] = -0.25f * l2[0] * area / (e[1] DOT e[0]);
            cornerareas[i][0] = -0.25f * l2[2] * area / (e[1] DOT e[2]);
            cornerareas[i][1] = area - cornerareas[i][2] - cornerareas[i][0];
        }
        else if (ew[2] <= 0.0f) {
            cornerareas[i][0] = -0.25f * l2[1] * area / (e[2] DOT e[1]);
            cornerareas[i][1] = -0.25f * l2[0] * area / (e[2] DOT e[0]);
            cornerareas[i][2] = area - cornerareas[i][0] - cornerareas[i][1];
        } else {
            TypeFunction ewscale = 0.5f * area / (ew[0] + ew[1] + ew[2]);
            for (uint8_t j = 0; j < 3; j++)
                cornerareas[i][j] = ewscale * (ew[(j+1)%3] + ew[(j+2)%3]);
        }
#pragma omp atomic
        areas[faces[i][0]] += cornerareas[i][0];
#pragma omp atomic
        areas[faces[i][1]] += cornerareas[i][1];
#pragma omp atomic
        areas[faces[i][2]] += cornerareas[i][2];
    }
}

//! -----------------------------------------------------------------------------
//! compute connectivity
//! -----------------------------------------------------------------------------

//! Find the direct neighbors of each vertex
std::vector<TypeIndexI> TriMesh::need_neighbors() {

    const bool verbose = false;

    // --------------------------------------------------------------------
    size_t nv = mVertices.size(), nf = mFaces.size();

    if (mVNeighbors.empty()) {

        if (verbose) {
            std::cout << "   > TriMesh::need_neighbors()...";
            fflush(stdout);
        }

        std::vector<TypeIndex> numneighbors(nv, 0);
        for (size_t i = 0; i < nf; i++) {
            numneighbors[mFaces[i][0]]++;
            numneighbors[mFaces[i][1]]++;
            numneighbors[mFaces[i][2]]++;
        }


        mVNeighbors.resize(nv);
        for (size_t i = 0; i < nv; i++) {
            mVNeighbors[i].reserve(numneighbors[i]);
            //mVNeighbors[i].reserve(mVNumNeighbors[i]+2); // Slop for boundaries
        }

        for (size_t i = 0; i < nf; i++) {
        for (size_t j = 0; j < 3; j++) {
            std::vector<TypeIndex> &me = mVNeighbors[mFaces[i][j]];
            TypeIndex n1 = mFaces[i][(j+1)%3];
            TypeIndex n2 = mFaces[i][(j+2)%3];
            if (std::find(me.begin(), me.end(), n1) == me.end())    me.push_back(n1);
            if (std::find(me.begin(), me.end(), n2) == me.end())    me.push_back(n2);
        }
        }

        if(verbose)
            std::cout << " Done!\n";
    }

    // --------------------------------------------------------------------
    // now linearize the data
    std::vector<TypeIndexI> linearized_retval;

    // first nv values will store the number of neighbors for each vertex!
        // we are not using numneighbors because they are approximate!
        // the counts can be wrong for boundary!
    linearized_retval.resize(nv, 0);

    for(size_t vidx = 0; vidx < nv; ++vidx) {

        const std::vector<TypeIndex> &vnbrs = mVNeighbors[vidx];
        linearized_retval[vidx] = vnbrs.size();
        std::copy(vnbrs.begin(), vnbrs.end(), std::back_inserter(linearized_retval));
    }
    return linearized_retval;
}

#if 0
//! Find the faces touching each vertex
void TriMesh::need_adjacentfaces(bool verbose) {

    if (!mVAdjFaces.empty())
        return;

    if (verbose) {
        std::cout << "   > TriMesh::need_adjacentfaces()...";
        fflush(stdout);
    }

    size_t nv = mVertices.size(), nf = mFaces.size();

    std::vector<TypeIndex> numadjacentfaces(nv);
    for (size_t i = 0; i < nf; i++) {
        numadjacentfaces[mFaces[i][0]]++;
        numadjacentfaces[mFaces[i][1]]++;
        numadjacentfaces[mFaces[i][2]]++;
    }

    mVAdjFaces.resize(mVertices.size());
    for (size_t i = 0; i < nv; i++)
        mVAdjFaces[i].reserve(numadjacentfaces[i]);

    for (size_t i = 0; i < nf; i++) {
    for (uint8_t j = 0; j < 3; j++)
        mVAdjFaces[mFaces[i][j]].push_back(i);
    }

    if(verbose)
        std::cout << " Done!\n";
}

//! Find the face across each edge from each other face (-1 on boundary)
//! If topology is bad, not necessarily what one would expect...
void TriMesh::need_across_edge(bool verbose) {

    if (!mFAcrossEdge.empty())
        return;

    need_adjacentfaces(verbose);
    if (mVAdjFaces.empty())
        return;

    if (verbose) {
        std::cout << "   > TriMesh::need_across_edge()...";
        fflush(stdout);
    }

    size_t nf = mFaces.size();
    mFAcrossEdge.resize(nf, Offset3(-1,-1,-1));

#pragma omp parallel for
    for (size_t i = 0; i < nf; i++) {
    for (uint8_t j = 0; j < 3; j++) {

        if (mFAcrossEdge[i][j] != -1)
            continue;

        TypeIndex v1 = mFaces[i][(j+1)%3];
        TypeIndex v2 = mFaces[i][(j+2)%3];

        const std::vector<TypeIndex> &a1 = mVAdjFaces[v1];
        const std::vector<TypeIndex> &a2 = mVAdjFaces[v2];

        for (size_t k1 = 0; k1 < a1.size(); k1++) {

            TypeIndex other = a1[k1];
            if (other == i)
                continue;
            std::vector<TypeIndex>::const_iterator it = find(a2.begin(), a2.end(), other);
            if (it == a2.end())
                continue;

            int vidx = (mFaces[other][0] == v1) ? 0 :
                       (mFaces[other][1] == v1) ? 1 :
                       (mFaces[other][2] == v1) ? 2 : -1;

            TypeIndex ind = (vidx+1)%3;
            if (mFaces[other][(ind+1)%3] != v2)
                continue;

            mFAcrossEdge[i][j] = other;
            mFAcrossEdge[other][ind] = i;
            break;
        }
    }
    }

    if(verbose)
        std::cout << " Done!\n";
}

//! Compute boundary edges
std::vector<TypeIndexI> TriMesh::need_boundary(bool verbose) {

    if (bedges.empty()) {

        need_across_edge(verbose);
        if (verbose) {
            std::cout << "   > TriMesh::need_boundary()...";
            fflush(stdout);
        }

        size_t nf = mFaces.size();

        bedges.clear();
        bedges.reserve(mFaces.size());

        for(size_t i = 0; i < nf; i++) {
        for(uint8_t j = 0; j < 3; j++) {

            if (mFAcrossEdge[i][j] == -1) {
                bedges.push_back((Edge( mFaces[i][(j+1)%3], mFaces[i][(j+2)%3] )));
            }
        }
        }

        bedges.shrink_to_fit();

        // orient these edges CCW
        size_t nbedges = bedges.size();
        for(size_t bidx = 0; bidx < nbedges-1; bidx++) {

            bool fnd = false;
            size_t jidx = bidx+1;
            for( ; jidx < nbedges; jidx++) {

                if (bedges[jidx][0] == bedges[bidx][1]) {
                    fnd = true;
                    break;
                }
            }

            if (fnd) {
                std::swap(bedges[bidx+1], bedges[jidx]);
            }

            // could mean that there are more than 1 boundary components
            // for now just carry on, since i eventually only need
            // boundary vertices for paramterization
            // so i dont care about storing the boundary component-wise
            else if (bedges[bidx][1] == bedges[0][0]) {
                continue;
            }
            else {

                //std::cout <<" failed for " << bidx << " : " <<bedges[bidx] << std::endl;

                std::ostringstream errMsg;
                errMsg << " " << this->tag() << "::need_boundary(): Could not CCW orient edges! returning possibly unoriented edges! (does the boundary have multiple components?)\n";
                std::cerr << errMsg.str();
                //throw std::invalid_argument(errMsg.str());

                if (verbose) {
                    std::cout << "   > TriMesh::need_boundary()...";
                }
                break;
            }
        }

        if(verbose)
            std::cout << " Done!\n";
    }

    return linearize<2,TypeIndex,TypeIndexI>(bedges, 2);
}
#endif


//! -----------------------------------------------------------------------------
//! Compute pairwise distances of ll vertuces
//! -----------------------------------------------------------------------------
std::vector<TypeFunction> TriMesh::get_pairwise_distances(uint8_t dim) const {

    if (dim < 1 || dim > this->mDim) {
        std::ostringstream errMsg;
        errMsg << " " << this->tag() << "::wrap_vertices("<<int(dim)<<"): Invalid dim specified for " << int(mDim) << "D vertices!" << std::endl;
        throw std::logic_error(errMsg.str());
    }

    const size_t nverts = this->mVertices.size();
    std::vector<TypeFunction> distances;
    distances.reserve(nverts*nverts/2);

    for (size_t i = 0; i < nverts; i++) {
    for (size_t j = i+1; j < nverts; j++) {

        TypeFunction dst = 0;

        auto a = this->mVertices[i];
        auto b = this->mVertices[j];
        for (size_t d = 0; d < dim; d++) {
            dst += (a[d]-b[d])*(a[d]-b[d]);
        }
        distances.push_back(std::sqrt(dst));
    }}
    return distances;
}

//! -----------------------------------------------------------------------------
//! Periodic Mesh
//! -----------------------------------------------------------------------------

bool TriMesh::wrap_vertices(uint8_t dim) {

    if (dim < 1 || dim > this->mDim) {
        std::ostringstream errMsg;
        errMsg << " " << this->tag() << "::wrap_vertices("<<int(dim)<<"): Invalid dim specified for " << int(mDim) << "D vertices!" << std::endl;
        throw std::logic_error(errMsg.str());
    }
    if (!this->mPeriodic) {
        std::ostringstream errMsg;
        errMsg << " " << this->tag() << "::wrap_vertices(): Mesh is not periodic!" << std::endl;
        throw std::invalid_argument(errMsg.str());
    }
    if (!this->bbox_valid) {
        std::ostringstream errMsg;
        errMsg << " " << this->tag() << "::wrap_vertices("<<int(dim)<<"): Bounding box not available!" << std::endl;
        throw std::logic_error(errMsg.str());
    }

    const Vertex boxw = this->mBox1 - this->mBox0;

    const size_t nverts = mVertices.size();
    for(size_t i=0; i<nverts; i++) {
    for(uint8_t d=0; d<dim; d++) {
        TypeFunction &x = mVertices[i][d];
        if (x <  this->mBox0[d])  x += boxw[d];
        if (x >= this->mBox1[d])  x -= boxw[d];
    }
    }
    return true;
}

void TriMesh::trim_periodicDelaunay(bool verbose) {

    if (!this->mPeriodic) {
        std::ostringstream errMsg;
        errMsg << " " << this->tag() << "::lift_delaunay(): Mesh is not periodic!" << std::endl;
        throw std::invalid_argument(errMsg.str());
    }

    this->mFaces.clear();
    this->mPeriodicFaces.clear();
    this->mTrimmedFaces.clear();
    this->mDuplicateVerts.clear();
    this->mDuplicateVertex_periodic.clear();

    if (verbose) {
        std::cout << "   > " << tag() << "::lift_delaunay()...";
        fflush(stdout);
    }

    const size_t noverts = this->mVertices.size();
    const size_t ndfaces = this->mDelaunayFaces.size();
    const Vertex boxw = mBox1 - mBox0;

    for(size_t i = 0; i < ndfaces; i++) {

        const std::vector<PeriodicVertex> &dface = this->mDelaunayFaces[i];

        Face face, tface;
        uint8_t num_orig_verts = 0;

        for(uint8_t d = 0; d < 3; d++) {

            const PeriodicVertex &pvertex = dface[d];

            const TypeIndex orig_vid = pvertex.origVidx;

            face[d] = orig_vid;
            tface[d] = orig_vid;

            // is an original vertex
            if (pvertex.is_original()) {
                num_orig_verts++;
                continue;
            }

            // needs to be duplicated
            auto iter = std::find(mDuplicateVertex_periodic.begin(), mDuplicateVertex_periodic.end(), pvertex);
            size_t ndupid = 0;

            if (iter != mDuplicateVertex_periodic.end()) {
                ndupid = noverts + std::distance(mDuplicateVertex_periodic.begin(), iter);
            }
            else {
                const int offx = pvertex.offsetx;
                const int offy = pvertex.offsety;

                // otherwise, let's duplicate
                Vertex dv (mVertices[orig_vid]);
                dv[0] += (offx == 0) ? 0.0 : float(offx) *boxw[0];
                dv[1] += (offy == 0) ? 0.0 : float(offy) *boxw[1];

                ndupid = noverts + mDuplicateVerts.size();
                mDuplicateVertex_periodic.push_back(pvertex);
                mDuplicateVerts.push_back(dv);
            }

            tface[d] = ndupid;
        }

        if (num_orig_verts == 3) {
            this->mFaces.push_back(face);
        }
        else {
            this->mPeriodicFaces.push_back(face);
            this->mTrimmedFaces.push_back(tface);
        }
    }

    if (verbose) {
        std::cout << " Done! created [" << mFaces.size() << ", " << mPeriodicFaces.size() << ", "<< mTrimmedFaces.size() << "] triangles "
                  << " with [" << mVertices.size() << ", " << mDuplicateVerts.size() << "] vertices!\n";
    }
}

//! -----------------------------------------------------------------------------
//! read/write off format
//! -----------------------------------------------------------------------------

bool TriMesh::read_off(const std::string &filename, std::vector<Vertex> &vertices, std::vector<Face> &faces, bool verbose) {

    if (verbose){
        std::cout << " TriMesh::read_off(" << filename << ")...";
        fflush(stdout);
    }

    std::ifstream file(filename);
    std::string str;

    std::getline(file, str);        // comment!

    int nverts, nfaces, c;
    std::getline(file, str);
    sscanf (str.c_str(), "%d %d %d", &nverts, &nfaces, &c);

    vertices.resize(nverts);
    for(int i = 0; i < nverts; i++) {
        std::getline(file, str);
        sscanf (str.c_str(), "%f %f %f", &vertices[i][0], &vertices[i][1], &vertices[i][2]);
    }

    int tmp;
    faces.resize(nfaces);
    for(int i = 0; i < nfaces; i++) {
        std::getline(file, str);
        sscanf (str.c_str(), "%d %d %d %d", &tmp, &faces[i][0], &faces[i][1], &faces[i][2]);
    }

    if (verbose) {
        std::cout << " Done! Read " << nverts << " points and " << nfaces << " faces!\n";
    }
    return true;
}

bool TriMesh::write_off(const std::string &filename, const std::vector<Vertex> &vertices, const std::vector<Face> &faces, const uint8_t &dim, bool verbose) {

    if (verbose){
        std::cout << " TriMesh::write_off(" << filename << ")...";
        fflush(stdout);
    }

    std::ofstream file(filename);

    size_t nverts = vertices.size();
    size_t nfaces = faces.size();

    file << "OFF\n";
    file << nverts <<" " << nfaces << " 0\n";

    if (dim == 2) {
        for(size_t i = 0; i < nverts; i++) {
            file << vertices[i][0] << " " << vertices[i][1] << " 0.0\n";
        }
    }
    else {
        for(size_t i = 0; i < nverts; i++) {
            file << vertices[i][0] << " " << vertices[i][1] << " " << vertices[i][2] << "\n";
        }
    }

    for(size_t i = 0; i < nfaces; i++) {
        file << "3 " << faces[i][0] << " " << faces[i][1] << " " << faces[i][2] << "\n";
    }

    if (verbose) {
        std::cout << " Done!\n";
    }
    return true;
}

//! -----------------------------------------------------------------------------
//! output binary format
//! -----------------------------------------------------------------------------

//static
bool TriMesh::write_binary(const std::string &fname,
                           const std::vector<Vertex> &vertices, const std::vector<Face> &faces,
                           const std::vector<std::string> &field_names,
                           const std::vector<std::vector<TypeFunction>*> &fields,
                           bool verbose) {

    const uint32_t nfields = fields.size();
    if (nfields != field_names.size()) {
        std::ostringstream errMsg;
        errMsg << " TriMesh::write_binary(): Got " << nfields << " fields, but " << field_names.size() << " field_names!\n";
        throw std::invalid_argument(errMsg.str());
    }

    // -------------------------------------------------------------------------
    if (verbose) {
        std::cout << "   > TriMesh::write_binary("<<fname<<")...";
        fflush(stdout);
    }

    FILE* outfile = fopen(fname.c_str(), "wb");

    uint32_t index_size = sizeof(uint32_t);
    uint32_t function_size = sizeof(float);
    uint32_t dummy_dimensions[3] = {1,1,1};

    fwrite(&index_size, sizeof(index_size), 1, outfile);
    fwrite(&function_size, sizeof(function_size), 1, outfile);
    fwrite(dummy_dimensions, sizeof(dummy_dimensions[0]), 3, outfile);

    // find the number of faces for each vertex
    std::vector<uint8_t> counts(vertices.size(),0);
    std::vector<uint32_t> first(vertices.size(), (uint32_t)-1);

    uint32_t count = 0;
    for (auto iter=faces.begin(); iter!=faces.end(); iter++) {

        const Face &f = *iter;

        counts[f[0]]++;
        counts[f[1]]++;
        counts[f[2]]++;

        first[f[0]] = std::min(first[f[0]], count);
        first[f[1]] = std::min(first[f[1]], count);
        first[f[2]] = std::min(first[f[2]], count);
        count++;
    }

    // now, actually write
    char token;
    TypeIndex path[2];

    count = 0;
    for (auto iter=faces.begin(); iter!=faces.end(); iter++) {

        const Face &face = *iter;

        token = 'v';

        // write vertices
        for (uint8_t k=0; k<3; k++) {

            const TypeIndex &vert = face[k];

            // if this is the first occurance of this vertex
            if (count == first[vert]) {

                fwrite(&token, sizeof(token), 1, outfile);
                fwrite(&vert, sizeof(TypeIndex), 1, outfile);
                fwrite(vertices[vert].data(), sizeof(TypeFunction), 3, outfile);

                for (uint8_t f = 0; f < nfields; f++) {
                    const std::vector<TypeFunction> &data = *(fields[f]);
                    fwrite(&(data[vert]), sizeof(TypeFunction), 1, outfile);
                }
            }
        }

        // now, write edges
        token = 'e';

        path[0] = face[0];      path[1] = face[1];
        fwrite(&token, sizeof(token), 1, outfile);
        fwrite(path, sizeof(TypeIndex), 2, outfile);

        path[0] = face[1];      path[1] = face[2];
        fwrite(&token, sizeof(token), 1, outfile);
        fwrite(path, sizeof(TypeIndex), 2, outfile);

        path[0] = face[2];      path[1] = face[0];
        fwrite(&token, sizeof(token), 1, outfile);
        fwrite(path, sizeof(TypeIndex), 2, outfile);

        // finally, write the face
        token = 'f';

        for (uint8_t k=0; k<3; k++) {

            counts[face[k]]--;
            if (counts[face[k]] == 0) {
                fwrite(&token, sizeof(token), 1, outfile);
                fwrite(&(face[k]), sizeof(TypeIndex), 1, outfile);
            }
        }

        count++;
    }
    fflush(outfile);
    fclose(outfile);
    if (verbose) {
        std::cout << " Done! Wrote " << vertices.size() << " vertices, "
                                     << faces.size() << " faces, and "
                                     << fields.size() << " fields!\n";

        if (fields.size() > 0){
            std::cout << "     > ";
            for(uint8_t i = 0; i < fields.size(); i++) {
                const std::vector<TypeFunction> &data = *fields[i];
                auto p = std::minmax_element(data.begin(), data.end());
                std::cout << "["<<field_names[i] << " : " << *p.first << ", " << *p.second << "] ";
            }
            std::cout << std::endl;
        }
    }
    return true;
}

bool TriMesh::write_binary(const std::string &fname, const std::string &filter_fields) {

    // -------------------------------------------------------------------------
    // addition to write only the relevant fields!
    std::vector<std::string> field_names;
    std::vector<std::vector<TypeFunction>*> fields;

    if (filter_fields.empty()) {
        for (auto iter=mFields.begin(); iter!=mFields.end(); ++iter) {
            field_names.push_back(iter->first);
            fields.push_back(&(iter->second));
        }
    }
    else {

        // hardcoded! also write point areas
        field_names.push_back("point_areas");
        fields.push_back(&(mFields["point_areas"]));

        for (auto iter=mFields.begin(); iter!=mFields.end(); ++iter) {
            if (iter->first.find(filter_fields) != std::string::npos) {
                field_names.push_back(iter->first);
                fields.push_back(&(iter->second));
            }
        }
    }

    // -------------------------------------------------------------------------
    if (!this->mPeriodic) {
        return TriMesh::write_binary(fname, this->mVertices, this->mFaces, field_names, fields, true);
    }
    else {
        std::vector<Face> mfaces = this->mFaces;
        mfaces.insert(mfaces.end(), this->mPeriodicFaces.begin(), this->mPeriodicFaces.end());
        return TriMesh::write_binary(fname, this->mVertices, mfaces, field_names, fields, true);
    }
}
//! -----------------------------------------------------------------------------
//! -----------------------------------------------------------------------------
