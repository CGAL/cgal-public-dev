// Copyright (c) 2007-2008  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Tong Zhao

#ifndef CGAL_COVARIANCE_MATRIX_3_H
#define CGAL_COVARIANCE_MATRIX_3_H

#include <CGAL/license/Implicit_surface_reconstruction_3.h>


#include <cmath>
#include <CGAL/array.h>
#ifdef CGAL_EIGEN3_ENABLED
#include <CGAL/Eigen_diagonalize_traits.h>
#else
#include <CGAL/Diagonalize_traits.h>
#endif
#include <cassert>

namespace CGAL {


/// \internal
/// The Covariance class represents a 3D covariance matrix.
/// The purpose of this class is to create a tensor for each point 
///
/// \cgalModels `Kernel::Vector_3`
///
/// @tparam Kernel   A traits class in the CGAL algorithms and data structures
template<class Kernel>
class Covariance_matrix_3
{
// Public types
public:

    typedef typename Kernel::FT       FT;
    typedef typename Kernel::Point_3  Point;
    typedef typename Kernel::Vector_3 Vector; ///< Kernel's Vector_3 class.
    typedef typename Kernel::Plane_3  Plane;
    typedef CGAL::cpp11::array<FT, 6>   Eigen_matrix;
    typedef CGAL::cpp11::array<FT, 3>   Eigen_vector;
    typedef CGAL::cpp11::array<FT, 9>   Eigen_three_vectors;

    #ifdef CGAL_EIGEN3_ENABLED
    typedef CGAL::Eigen_diagonalize_traits<FT, 3>           Diagonalize_traits;
    #else
    typedef CGAL::Diagonalize_traits<FT, 3>                 Diagonalize_traits;
    #endif

// Private variables
private:

     Eigen_matrix           m_tensor;
     Eigen_vector           m_eigen_values;
     Eigen_three_vectors    m_eigen_vects;

// Public methods
public:

    /// m_tensor is I_3 by default
    Covariance_matrix_3()
    {
      set_id();
    }

    Covariance_matrix_3(const Vector& normal, const FT anisotropy)
    {
      if(normal == CGAL::NULL_VECTOR)
        this->set_id();
      else{
        m_tensor[0] = normal[0] * normal[0];
        m_tensor[1] = normal[0] * normal[1];
        m_tensor[2] = normal[0] * normal[2];
        m_tensor[3] = normal[1] * normal[1];
        m_tensor[4] = normal[1] * normal[2];
        m_tensor[5] = normal[2] * normal[2];

        diagonalize();

        // FT inv_anisotropy = 1. / std::sqrt(anisotropy);
        FT inv_anisotropy = 1;
        build_from_eigen(eigen_vect(0), eigen_vect(1), eigen_vect(2), inv_anisotropy, inv_anisotropy, anisotropy);
      }
    }

    Covariance_matrix_3(const Point& p, const Vector& normal, const FT anisotropy)
    {
      if(normal == CGAL::NULL_VECTOR || std::sqrt(normal * normal) < 0.9)
        this->set_id();
      else{
        Plane tangent_plane(p, normal);
        Vector vmin = tangent_plane.base1();
        Vector vmid = tangent_plane.base2();
        Vector vmax = normal;

        // Normalize vectors
        vmin = vmin / std::sqrt(vmin * vmin);
        vmid = vmid / std::sqrt(vmid * vmid);
        vmax = vmax / std::sqrt(vmax * vmax);

        //FT inv_anisotropy = 1. / std::sqrt(anisotropy);
        FT inv_anisotropy = 1.;
        build_from_eigen(vmin, vmid, vmax, inv_anisotropy, inv_anisotropy, anisotropy);
      } 
    }

    Covariance_matrix_3(const Covariance_matrix_3& covariance)
    {
      m_tensor = covariance.tensors();
      m_eigen_values = covariance.eigen_values();
      m_eigen_vects = covariance.eigen_vects();
    }

    Covariance_matrix_3(const Covariance_matrix_3& c1, const Covariance_matrix_3& c2, const bool convert = true)
    {
      if(!convert){
        for(unsigned int i = 0; i < 6; i++)
          m_tensor[i] = 0.5 * (c1.tensor(i) + c2.tensor(i));
        diagonalize();
      }
      else{
        Covariance_matrix_3 log_c1, log_c2;
        log_c1.build_from_eigen(c1.eigen_vect(0), c1.eigen_vect(1), c1.eigen_vect(2),
             log(c1.eigen_value(0)), log(c1.eigen_value(1)), log(c1.eigen_value(2)));
        log_c2.build_from_eigen(c2.eigen_vect(0), c2.eigen_vect(1), c2.eigen_vect(2),
             log(c2.eigen_value(0)), log(c2.eigen_value(1)), log(c2.eigen_value(2)));
            
        Covariance_matrix_3 log_c12;
        for(unsigned int i = 0; i < 6; i++)
          log_c12.tensor(i) = 0.5 * (log_c1.tensor(i) + log_c2.tensor(i));
        log_c12.diagonalize();

        build_from_eigen(log_c12.eigen_vect(0), log_c12.eigen_vect(1), log_c12.eigen_vect(2),
              exp(log_c12.eigen_value(0)), exp(log_c12.eigen_value(1)), exp(log_c12.eigen_value(2)));
      }
    }

    Covariance_matrix_3& operator+(const Covariance_matrix_3& other)
    {
      if(this != &other){
        Covariance_matrix_3 log_this, log_other;
        log_this.build_from_eigen(this->eigen_vect(0), this->eigen_vect(1), this->eigen_vect(2),
                log(this->eigen_value(0)), log(this->eigen_value(1)), log(this->eigen_value(2)));
        log_other.build_from_eigen(other.eigen_vect(0), other.eigen_vect(1), other.eigen_vect(2),
                log(other.eigen_value(0)), log(other.eigen_value(1)), log(other.eigen_value(2)));

        Covariance_matrix_3 log_sum;
        for(unsigned int i = 0; i < 6; i++)
          log_sum.tensor(i) = 0.5 * (log_this.tensor(i) + log_other.tensor(i));
        log_sum.diagonalize();

        this->build_from_eigen(log_sum.eigen_vect(0), log_sum.eigen_vect(1), log_sum.eigen_vect(2),
              exp(log_sum.eigen_value(0)), exp(log_sum.eigen_value(1)), exp(log_sum.eigen_value(2)));
      }

      return *this;
    }

    ~Covariance_matrix_3()
    {
    }

    void set_id()
    {
      m_tensor = make_array(1., 0., 0., 1., 0., 1.);
      m_eigen_values = make_array(1., 1., 1.);
      // make_array takes at most 6 
      m_eigen_vects[0] = 0.;
      m_eigen_vects[1] = 0.;
      m_eigen_vects[2] = 1.;
      m_eigen_vects[3] = 0.;
      m_eigen_vects[4] = 1.;
      m_eigen_vects[5] = 0.;
      m_eigen_vects[6] = 1.;
      m_eigen_vects[7] = 0.;
      m_eigen_vects[8] = 0.;
    }

    bool diagonalize()
    {
      if(!(Diagonalize_traits::diagonalize_selfadjoint_covariance_matrix(m_tensor,
                                                                         m_eigen_values,
                                                                         m_eigen_vects)))
      {
        std::cerr << "Error: cannot diagonalize matrix" << std::endl;
        return false;
      }

      return true;
    }

    FT& tensor(const unsigned int i) 
    {
      if(i < 6)
        return m_tensor[i]; 
      else {
        std::cerr << "Error: index " << i << " is out of range." << std::endl;
        return m_tensor[0];
      }
    }

    const FT& tensor(const unsigned int i) const
    {
      if(i < 6)
        return m_tensor[i]; 
      else {
        std::cerr << "Error: index " << i << " is out of range." << std::endl;
        return m_tensor[0];
      }
    }

    FT& eigen_value(const unsigned int i) 
    {
      if(i < 3)
        return m_eigen_values[i]; 
      else {
        std::cerr << "Error: index " << i << " is out of range." << std::endl;
        return m_eigen_values[0];
      }
    }

    const FT& eigen_value(const unsigned int i) const
    {
      if(i < 3)
        return m_eigen_values[i]; 
      else {
        std::cerr << "Error: index " << i << " is out of range." << std::endl;
        return m_eigen_values[0];
      }
    }

    const Vector eigen_vect(const unsigned int i) const
    {
      Vector m_vect;

      if(i < 3)
        m_vect = Vector(m_eigen_vects[i * 3], m_eigen_vects[i * 3 + 1], m_eigen_vects[i * 3 + 2]);
      else {
        std::cerr << "Error: index " << i << " is out of range." << std::endl;
        m_vect = Vector(-1., -1., -1.);
      }
      return m_vect;
    }

    Eigen_matrix& tensors() { return m_tensor;}
    const Eigen_matrix& tensors() const { return m_tensor;}

    Eigen_vector& eigen_values() { return m_eigen_values;}
    const Eigen_vector& eigen_values() const { return m_eigen_values;}

    Eigen_three_vectors& eigen_vects() { return m_eigen_vects;}
    const Eigen_three_vectors& eigen_vects() const { return m_eigen_vects;}

    FT isotropy() { return m_eigen_values[0] / m_eigen_values[2];}
    FT anisotropy() { return 1.0 - isotropy();}
    FT trace() { return m_eigen_values[0] * m_eigen_values[1] * m_eigen_values[2];}
    bool isotropic() { return isotropy() == 1.;}

    void normalize(const FT factor)
    {
      double iso = isotropy();
      assert(iso >= 0.0 && iso <= 1.0);
      if(iso != 1.0)
      {
        double emin = iso;
        double emax = factor - (factor - 1.0) * iso;
        double emid = m_eigen_values[1] / m_eigen_values[0] * emax;

        Vector vmin = eigen_vect(0);
        Vector vmid = eigen_vect(1);
        Vector vmax = eigen_vect(2);
        
        build_from_eigen(vmin, vmid, vmax, emin, emid, emax);
      }
    }

    void normalize(const FT fmin, const FT fmid, const FT fmax)
    {
      Vector vmin = eigen_vect(0);
      Vector vmid = eigen_vect(1);
      Vector vmax = eigen_vect(2);
      build_from_eigen(vmin, vmid, vmax, fmin, fmid, fmax);
    }

    void build_from_eigen(const Vector& vmin, const Vector& vmid, const Vector& vmax,
                          const FT fmin, const FT fmid, const FT fmax)
    {
      m_tensor[0] = fmin * vmin.x() * vmin.x() + fmid * vmid.x() * vmid.x() + fmax * vmax.x() * vmax.x();
      m_tensor[1] = fmin * vmin.x() * vmin.y() + fmid * vmid.x() * vmid.y() + fmax * vmax.x() * vmax.y();
      m_tensor[2] = fmin * vmin.x() * vmin.z() + fmid * vmid.x() * vmid.z() + fmax * vmax.x() * vmax.z();
      m_tensor[3] = fmin * vmin.y() * vmin.y() + fmid * vmid.y() * vmid.y() + fmax * vmax.y() * vmax.y();
      m_tensor[4] = fmin * vmin.y() * vmin.z() + fmid * vmid.y() * vmid.z() + fmax * vmax.y() * vmax.z();
      m_tensor[5] = fmin * vmin.z() * vmin.z() + fmid * vmid.z() * vmid.z() + fmax * vmax.z() * vmax.z();
      
		  diagonalize();
    }

    FT ut_c_v(const Vector& u, const Vector& v)
    {
      Vector ut_c = Vector(u.x() * m_tensor[0] + u.y() * m_tensor[1] + u.z() * m_tensor[2],
                           u.x() * m_tensor[1] + u.y() * m_tensor[3] + u.z() * m_tensor[4], 
                           u.x() * m_tensor[2] + u.y() * m_tensor[4] + u.z() * m_tensor[5]);

      FT dot = ut_c.x() * v.x() + ut_c.y() * v.y() + ut_c.z() * v.z();

      return dot;
    }

};


} //namespace CGAL

#endif //CGAL_COVARIANCE_MATRIX_3_H
