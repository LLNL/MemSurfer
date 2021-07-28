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

#ifndef _DISTANCE_KERNELS_H_
#define _DISTANCE_KERNELS_H_

#include <cmath>
#include <vector>

#include "Types.hpp"

//! ----------------------------------------------------------------------------
//!
//! \brief This file provides kernels to compute distances
//!
//! ----------------------------------------------------------------------------

//! ----------------------------------------------------------------------------
//! the (abstract) base class for distance kernel
//! ----------------------------------------------------------------------------
class DistanceKernel {
public:
    DistanceKernel() {}
    virtual ~DistanceKernel() {}

    virtual TypeFunction operator()(TypeFunction ax, TypeFunction ay, TypeFunction bx, TypeFunction by) const = 0;
    virtual TypeFunction operator()(TypeFunction ax, TypeFunction ay, TypeFunction az, TypeFunction bx, TypeFunction by, TypeFunction bz) const = 0;
    TypeFunction operator()(const TypeFunction *a, const TypeFunction *b, const uint8_t &dim) const {
        return (dim == 2) ? this->operator ()(a[0], a[1], b[0], b[1]) :
               (dim == 3) ? this->operator ()(a[0], a[1], a[2], b[0], b[1], b[2]) : -1;
    }

    //! ------------------------------------------------------------------------
protected:
    static bool parse_bbox(float *_, int n, Vertex &bbox0, Vertex &bbox1) {

      if (n == 2) {
          bbox0 = Vertex(0, 0, 0);
          bbox1 = Vertex(_[0], _[1], 0);
      }
      else if (n == 3) {
          bbox0 = Vertex(0, 0, 0);
          bbox1= Vertex(_[0], _[1], _[2]);
      }
      else if (n == 4) {
          bbox0 = Vertex(_[0], _[1], 0);
          bbox1 = Vertex(_[2], _[3], 0);
      }
      else if (n == 6) {
          bbox0 = Vertex(_[0], _[1], _[2]);
          bbox1 = Vertex(_[3], _[4], _[5]);
      }
      else {
        return false;
      }
      return true;
    }
    //! ------------------------------------------------------------------------
};

//! ----------------------------------------------------------------------------
//! Square euclidean distance
//! ----------------------------------------------------------------------------
class DistanceSquared : public DistanceKernel {

public:
    DistanceSquared() {}
    TypeFunction operator()(TypeFunction ax, TypeFunction ay, TypeFunction bx, TypeFunction by) const {
        return (ax - bx)*(ax - bx) + (ay - by)*(ay - by);
    }
    TypeFunction operator()(TypeFunction ax, TypeFunction ay, TypeFunction az, TypeFunction bx, TypeFunction by, TypeFunction bz) const {
        return this->operator ()(ax, ay, bx, by) + (az - bz)*(az - bz);
    }
};

//! ----------------------------------------------------------------------------
//! Square euclidean distance on periodic domain in xy
//! ----------------------------------------------------------------------------
class DistancePeriodicXYSquared : public DistanceKernel {

    Vertex mBox0, mBox1, mBoxw;

public:
    DistancePeriodicXYSquared(const Vertex& bbox0, const Vertex& bbox1) {
        mBox0 = bbox0;  mBox1 = bbox1;
        mBoxw = mBox1 - mBox0;
        return;
        std::cout << " Initializing DistancePeriodicXYSquared. Domain = "
                  << " [" << mBox0[0] << ", " << mBox0[1] << "] -- "
                  << " [" << mBox1[0] << ", " << mBox1[1] << "], width = "
                  << " [" << mBoxw[0] << ", " << mBoxw[1] << "]\n";
    }

    DistancePeriodicXYSquared(float *_, int n) {
      if (!DistanceKernel::parse_bbox(_, n, mBox0, mBox1)) {
          std::ostringstream errMsg;
          errMsg << " DistancePeriodicXYSquared(): Invalid bounding box!" << std::endl;
          throw std::invalid_argument(errMsg.str());
      }
      this->mBoxw = mBox1 - mBox0;
      return;

      std::cout << " Initializing DistancePeriodicXYSquared. Domain = "
                << " [" << mBox0[0] << ", " << mBox0[1] << ", " << mBox0[2] << "] -- "
                << " [" << mBox1[0] << ", " << mBox1[1] << ", " << mBox1[2] << "], width = "
                << " [" << mBoxw[0] << ", " << mBoxw[1] << ", " << mBoxw[2] << "]\n";
    }

    TypeFunction operator()(TypeFunction ax, TypeFunction ay, TypeFunction bx, TypeFunction by) const {

        if (fabs(ax - bx) >= 0.5*mBoxw[0]) {
            if (ax < bx)    ax += mBoxw[0];
            else            bx += mBoxw[0];
        }

        if (fabs(ay - by) >= 0.5*mBoxw[1]) {
            if (ay < by)    ay += mBoxw[1];
            else            by += mBoxw[1];
        }
        return (ax - bx)*(ax - bx) + (ay - by)*(ay - by);
    }
    TypeFunction operator()(TypeFunction ax, TypeFunction ay, TypeFunction az, TypeFunction bx, TypeFunction by, TypeFunction bz) const {
        return this->operator ()(ax, ay, bx, by) + (az - bz)*(az - bz);
    }
};

//! ---------------------------------------------------------------------------------------------
//! Square euclidean distance on periodic domain in xyz
class DistancePeriodicSquared : public DistanceKernel {

    Vertex mBox0, mBox1, mBoxw;

public:
    DistancePeriodicSquared(const Vertex& bbox0, const Vertex& bbox1) {
        mBox0 = bbox0;  mBox1 = bbox1;
        mBoxw = mBox1 - mBox0;
    }
    DistancePeriodicSquared(float *_, int n) {
      if (!DistanceKernel::parse_bbox(_, n, mBox0, mBox1)) {
          std::ostringstream errMsg;
          errMsg << " DistancePeriodicSquared(): Invalid bounding box!" << std::endl;
          throw std::invalid_argument(errMsg.str());
      }
      this->mBoxw = mBox1 - mBox0;
    }
    TypeFunction operator()(TypeFunction ax, TypeFunction ay, TypeFunction bx, TypeFunction by) const {

        if (fabs(ax - bx) >= 0.5*mBoxw[0]) {
            if (ax < bx)    ax += mBoxw[0];
            else            bx += mBoxw[0];
        }

        if (fabs(ay - by) >= 0.5*mBoxw[1]) {
            if (ay < by)    ay += mBoxw[1];
            else            by += mBoxw[1];
        }

        return (ax - bx)*(ax - bx) + (ay - by)*(ay - by);
    }
    TypeFunction operator()(TypeFunction ax, TypeFunction ay, TypeFunction az, TypeFunction bx, TypeFunction by, TypeFunction bz) const {

        if (fabs(az - bz) >= 0.5*mBoxw[2]) {
            if (az < bz)    az += mBoxw[2];
            else            bz += mBoxw[2];
        }

        return this->operator ()(ax, ay, bx, by) + (az - bz)*(az - bz);
    }
};
//! ---------------------------------------------------------------------------------------------

#endif /* _DISTANCE_KERNELS_H_ */
