// ============================================================================
//
// Copyright (c) 2001-2006 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of EXACUS (http://www.mpi-inf.mpg.de/projects/EXACUS/);
// you may redistribute it under the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with EXACUS.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// ----------------------------------------------------------------------------
//
// Library       : QdX
// File          : demos/xtrinagulate/include/math_routines.h
// QdX_release   : $Name:  $
// Revision      : $Revision: 1.2 $
// Revision_date : $Date: 2009-07-24 13:21:30 $
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

#ifndef MATH_ROUTINES_H 
#define MATH_ROUTINES_H

#include "includes_common.h"
#include <CGAL/Handle_for.h>

#include <math.h>
#include <cmath>

template <class NT>
class Vector_3_rep {

public:
    Vector_3_rep() {
        _m_data[0] = NT(0);
        _m_data[1] = NT(0);
        _m_data[2] = NT(0);
    }

    Vector_3_rep(const NT& x, const NT& y, const NT& z) {
        _m_data[0] = x;
        _m_data[1] = y;
        _m_data[2] = z;
    }

    NT _m_data[3];
};

struct Quaternion;

template <class NT_>
class Vector_3 :
    public CGAL::Handle_for< Vector_3_rep< NT_ > > {

public:
    
    //! this instance's template parameter
    typedef NT_ NT;
    
    //! represenation type
    typedef Vector_3_rep< NT > Rep;
    
    //! myself
    typedef Vector_3<NT> Self;

    //! the handle superclass
    typedef ::CGAL::Handle_for<Rep> Base;
    
    Vector_3() : 
        Base(Rep()) {  
    }

    Vector_3(const Self& p) : 
        Base(static_cast<const Base&>(p)) {  
    }

    Vector_3(const boost::array< NT, 3 >& v) :
        Base(Rep(v[0], v[1], v[2])) {
    }

    Vector_3(const NT& x, const NT& y, const NT& z) :
        Base(Rep(x, y, z)) {
    }

    inline const NT& operator[](unsigned i) const {
        CGAL_precondition(i <= 2);
        return this->ptr()->_m_data[i];
    }

    inline NT& operator[](unsigned i) {
        CGAL_precondition(i <= 2);
        return this->ptr()->_m_data[i];
    }
           
    inline const NT& x() const {
        return this->ptr()->_m_data[0];
    }

    inline const NT& y() const {
        return this->ptr()->_m_data[1];
    } 

    inline const NT& z() const {
        return this->ptr()->_m_data[2];
    }

    NT square_magnitude() const {
        return (x()*x() + y()*y() + z()*z());
    }

    Self normalize() {
        NT mag = square_magnitude();
        if(mag == NT(0))    
            return *this;
        mag = sqrt(mag);
        (*this)[0] /= mag;
        (*this)[1] /= mag;
        (*this)[2] /= mag;
        return *this;
    }
           
    //double distance(const Vector3d& vec) const;
    NT dot_product(const Self& v) const {
        return (x()*v.x() + y()*v.y() + z()*v.z());
    }

    Self cross_product(const Self& v) const {
        Self res;        
        res[0] = y()*v.z() - z()*v.y();
        res[1] = z()*v.x() - x()*v.z();
        res[2] = x()*v.y() - y()*v.x();
        return res;
    }

    friend Self normalize(const Self& v) {
        Self tmp(v);
        tmp.copy_on_write();
        tmp.normalize();
        return tmp;
    }

    friend Self operator +(const Self& v1, const Self& v2) {
        Self res;
        res[0] = v1.x() + v2.x();
        res[1] = v1.y() + v2.y();
        res[2] = v1.z() + v2.z();
        return res;
    } 

    friend Self operator -(const Self& v1, const Self& v2) {
        Self res;
        res[0] = v1.x() - v2.x();
        res[1] = v1.y() - v2.y();
        res[2] = v1.z() - v2.z();
        return res;
    } 

    Self operator *(const NT& scalar) const {
        Self res;
        res[0] = x() * scalar;
        res[1] = y() * scalar;
        res[2] = z() * scalar;
        return res;
    }

    Self operator *=(const NT& scalar) {
        (*this)[0] *= scalar;
        (*this)[1] *= scalar;
        (*this)[2] *= scalar;
        return *this;
    }

    Self operator /=(const NT& scalar) {
        CGAL_precondition(scalar != NT(0));
        (*this)[0] /= scalar;
        (*this)[1] /= scalar;
        (*this)[2] /= scalar;
        return *this;
    }

    Self operator-() {
        (*this)[0] = -x();
        (*this)[1] = -y();
        (*this)[2] = -z();
        return *this;    
    }

    bool operator ==(const Self& v) const {
        return (x() == v.x() && y() == v.y() && z() == v.z());
    }

    bool operator !=(const Self& v) const {
        return !(operator ==(v));
    }

    friend std::ostream& operator <<(std::ostream& os, const Self& v) {
        os << "[" << v[0] << "; " << v[1] << "; " << v[2] << "]";
        return os;
    }

    friend struct Quaternion;
};

typedef Vector_3< float > Vector_3f;
typedef Vector_3< double > Vector_3d;

struct Quaternion {
    Vector_3d v;
    double w;

    void copy_on_write() {
        v.copy_on_write();
    }
};

void quaternion_from_axis(Quaternion& q, const Vector_3d& axis, double angle);
    
void quaternion_trackball(Quaternion& q, const Point_2d& p1,
        const Point_2d& p2);

void quaternion_2_matrix(const Quaternion& q, double *mat);

void quaternion_multiply(Quaternion& res, const Quaternion& p, 
        const Quaternion& q);

void quaternion_normalize(Quaternion& q);

#endif // MATH_ROUTINES_H

