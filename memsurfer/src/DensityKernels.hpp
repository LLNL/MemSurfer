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

#ifndef _DENSITY_KERNELS_H_
#define _DENSITY_KERNELS_H_

#include <cmath>
#include <vector>

#include "Types.hpp"

//! ----------------------------------------------------------------------------
//!
//! \brief This file provides kernels to compute density
//!
//! ----------------------------------------------------------------------------

//! ----------------------------------------------------------------------------
//! the (abstract) base class for density kernel
//! ----------------------------------------------------------------------------
class DensityKernel {

public:
    DensityKernel() {}
    virtual ~DensityKernel() {}

    virtual TypeFunction operator()(const TypeFunction &) const = 0;
};

//! ----------------------------------------------------------------------------
//! a generic (abstract) gaussian kernel
//! ----------------------------------------------------------------------------
class GaussianKernel : public DensityKernel {

    const TypeFunction efactor;
    const TypeFunction sfactor;

public:
    GaussianKernel(const TypeFunction e, const TypeFunction s) :
        efactor(e), sfactor(s) {}

    TypeFunction operator()(const TypeFunction &xsquared) const {
      return sfactor * exp(xsquared*efactor);
    }
};

//! ----------------------------------------------------------------------------
//! 1D 2D ad 3D implementations of Gaussian kernels
//! ----------------------------------------------------------------------------
class GaussianKernel1D : public GaussianKernel {
public:
    GaussianKernel1D(const TypeFunction &sigma) :
        GaussianKernel(-1.0/(2.0*sigma*sigma),
                        1.0/(std::sqrt(2.0*M_PI)*sigma))  // 1 / (sqrt(2pi) s)
    {}
};

class GaussianKernel2D : public GaussianKernel {
public:
    GaussianKernel2D(const TypeFunction &sigma) :
        GaussianKernel(-1.0/(2.0*sigma*sigma),
                        1.0/(2.0*M_PI*sigma*sigma))     // 1 / (2pi s^2)    == 1 / (sqrt(2pi) s)^2
    {}
};

class GaussianKernel3D : public GaussianKernel {

public:
    GaussianKernel3D(const TypeFunction &sigma) :
        GaussianKernel(-1.0/(2.0*sigma*sigma),
                        1.0/std::pow(std::sqrt(2.0*M_PI)*sigma, 3.0))  // 1 / (sqrt(2pi) s)^3
    {}
};

//! ----------------------------------------------------------------------------
//! ----------------------------------------------------------------------------
/*
// commented these kernels on Mar 22, 2019
// and instead, implemented correct Gaussian kernels for different dimensions
#if 1
class DensityKernel {

    const TypeFunction efactor;
    const TypeFunction sfactor;

public:

    //! note: this is a 2D kernel!
    //! this should be parameterized to support 1d and 3d
    DensityKernel(const TypeFunction &sigma) :
        efactor(-1.0/(2.0*sigma*sigma)),
        sfactor( 1.0/(2.0*M_PI*sigma*sigma)) {}

    TypeFunction operator()(const TypeFunction &xsquared) const {
        return sfactor * exp(xsquared*efactor);
    }
};
#else
//! ----------------------------------------------------------------------------
//! On Feb 05, 2020
//! Harsh added this simple kernel to reproduce the functionality
//! of the code for Brain Paper
//!
//! see, e.g., Timo's code
//! https://lc.llnl.gov/bitbucket/projects/KRAS/repos/densityestimation/browse/src/DensityEstimation.h
//! ----------------------------------------------------------------------------
class DensityKernel {

public:
  DensityKernel(TypeFunction sigma, bool fake=true) : mSigmaSquare(sigma*sigma) {}
  virtual TypeFunction operator()(TypeFunction x) const {return exp(-x / mSigmaSquare);}

private:
  TypeFunction mSigmaSquare;
};
#endif
//! ----------------------------------------------------------------------------
*/

//! ---------------------------------------------------------------------------------------------

#endif /* _DENSITY_KERNELS_H_ */
