
// Copyright (c) 2010  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: $
// $Id: $
//
//
// Author(s)     : Asaf Porat          <asafpor1@post.tau.ac.il>

#ifndef LINES_THROUGH_SEGMENTS_POINT_ADAPT_H
#define LINES_THROUGH_SEGMENTS_POINT_ADAPT_H

#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_rational_function_traits_2.h>
#include <CGAL/Algebraic_kernel_d_1.h>
// #include <CGAL/Algebraic_kernel_2_1.h>
#include <CGAL/Lazy_exact_nt.h>
/*! \file
 *
 * This class adapts the representation of Point_2 at the algebraic traits.
 *
 */

namespace CGAL {


template <typename Traits_3_, typename Point_, typename Number_type_>
class Lines_through_segments_point_adapt_2
{
public:
   typedef Traits_3_ Traits_3;

   typedef enum {
      LTS_VERTEX_INTERIOR = 0,
      LTS_VERTEX_X_MINUS_INFINITY = 1,
      LTS_VERTEX_X_PLUS_INFINITY = 2,
      LTS_VERTEX_Y_MINUS_INFINITY = 4,
      LTS_VERTEX_Y_PLUS_INFINITY = 8
   } Point_type;

private:
  typedef typename Traits_3::Rational_kernel            Rational_kernel;
  typedef typename Traits_3::Alg_kernel                 Alg_kernel;

public:
  typedef typename Rational_kernel::FT                  Rational;
  typedef typename Alg_kernel::FT                       Algebraic;

private:
  typedef Lines_through_segments_point_adapt_2          Self;

  typedef typename Rational_kernel::Segment_3           Rational_segment_3;



  /* Specific typedefs for conic arc traits. */
  typedef CGAL::CORE_algebraic_number_traits            Nt_traits;
  typedef CGAL::Arr_conic_traits_2<Rational_kernel, Alg_kernel, Nt_traits>
  Conic_traits_arr_on_plane_2_;
  typedef CGAL::Arr_consolidated_curve_data_traits_2<Conic_traits_arr_on_plane_2_,const Rational_segment_3*>
  Conic_traits_arr_on_plane_2;
  typedef CORE::BigInt                                  Integer;
   /* Specific typedefs for Rational arc traits . */
#if USE_SQRT_TRAITS
   typedef CGAL::Algebraic_kernel_2_1<Rational>	   AK1;
#else
  typedef CGAL::Algebraic_kernel_d_1<Integer>	   AK1;
#endif
  typedef CGAL::Arr_rational_function_traits_2<AK1>  Rational_arc_traits_arr_on_plane_2_;
  typedef CGAL::Arr_consolidated_curve_data_traits_2<Rational_arc_traits_arr_on_plane_2_,const Rational_segment_3*>
  Rational_arc_traits_arr_on_plane_2;

  typedef typename Rational_kernel::Point_2         Rational_point_2;

protected:
  typedef Number_type_                                  Number_type;


protected:
  /* If the point was created with rational representation keep it. */
  Rational m_rat_x;
  Rational m_rat_y;
  bool m_rat_rep_exists;

  Number_type m_x;
  Number_type m_y;
  Point_ m_orig_point;
  bool m_original_point_created; /* True if m_orig_point was intialized. */

  Point_type m_point_type;

public:

  typedef Point_ Point;
  Lines_through_segments_point_adapt_2()
  {
    m_rat_rep_exists = false;
    m_point_type = LTS_VERTEX_INTERIOR;
  }

  Lines_through_segments_point_adapt_2(const typename Rational_arc_traits_arr_on_plane_2::Point_2& p,
                                       const Rational& rat_x,
                                       const Rational& rat_y)
  {
    Lines_through_segments_traits_on_plane_adapt<
    Traits_3>  traits_2_adapt;
    m_rat_rep_exists = true;
    m_point_type = LTS_VERTEX_INTERIOR;

    m_rat_x = rat_x;
    m_rat_y = rat_y;
    m_orig_point = p;
    traits_2_adapt.convert_point(p,m_x,m_y);
    m_original_point_created = true;
  }

  Lines_through_segments_point_adapt_2(const typename Conic_traits_arr_on_plane_2::Point_2& p,
                                       const Rational& rat_x,
                                       const Rational& rat_y)
  {
    m_rat_rep_exists = true;
    m_point_type = LTS_VERTEX_INTERIOR;

    m_rat_x = rat_x;
    m_rat_y = rat_y;

    m_orig_point = p;
    m_x = p.x();
    m_y = p.y();
    m_original_point_created = true;
  }


  Lines_through_segments_point_adapt_2(const typename Rational_arc_traits_arr_on_plane_2::Point_2& p)
  {
    Lines_through_segments_traits_on_plane_adapt<
    Traits_3>   traits_2_adapt;
    m_rat_rep_exists = false;
    m_point_type = LTS_VERTEX_INTERIOR;

    m_rat_x = 0;
    m_rat_y = 0;
    m_orig_point = p;
    traits_2_adapt.convert_point(p,m_x,m_y);
    m_original_point_created = true;
  }

  Lines_through_segments_point_adapt_2(const typename Conic_traits_arr_on_plane_2::Point_2& p)
  {
    m_rat_rep_exists = false;
    m_point_type = LTS_VERTEX_INTERIOR;

    m_rat_x = 0;
    m_rat_y = 0;

    m_orig_point = p;
    m_x = p.x();
    m_y = p.y();
    m_original_point_created = true;
  }


  Lines_through_segments_point_adapt_2(const Self& p)
  {
    m_rat_rep_exists = p.m_rat_rep_exists;
    m_point_type = p.m_point_type;

    m_rat_y = p.m_rat_y;
    m_rat_x = p.m_rat_x;
    m_orig_point = p.m_orig_point;
    m_x = p.m_x;
    m_y = p.m_y;
    m_original_point_created = p.m_original_point_created;
  }

  template <typename N_x,typename N_y>
  Lines_through_segments_point_adapt_2(const N_x& x,const N_y& y)
  {
    m_rat_rep_exists = false;
    m_point_type = LTS_VERTEX_INTERIOR;

    m_x = x;
    m_y = y;
  }

  Lines_through_segments_point_adapt_2(const typename Rational_arc_traits_arr_on_plane_2::Algebraic_real_1& number,
                                       Point_type ver_type)
  {
    Lines_through_segments_traits_on_plane_adapt<
    Traits_3>  traits_2_adapt;

    m_rat_rep_exists = false;
    m_point_type = ver_type;
    switch (ver_type)
    {
     case LTS_VERTEX_X_PLUS_INFINITY:
     case LTS_VERTEX_X_MINUS_INFINITY:
      {
       m_y = traits_2_adapt.convert_real_to_algebraic(number);

      }
      break;
     case LTS_VERTEX_Y_PLUS_INFINITY:
     case LTS_VERTEX_Y_MINUS_INFINITY:
      {
       m_x = traits_2_adapt.convert_real_to_algebraic(number);
      }
      break;
     default:
      CGAL_error_msg("LTS_ERROR");
    }
  }


  Lines_through_segments_point_adapt_2(const Rational& x, const Rational& y)
  {
    m_x = x;
    m_y = y;
    m_rat_x = x;
    m_rat_y = y;
    m_rat_rep_exists = true;
    m_point_type = LTS_VERTEX_INTERIOR;
    m_original_point_created = true;
    m_orig_point = Point(x,y);
  }

  Number_type x() const
  {
    CGAL_assertion(m_point_type != LTS_VERTEX_X_MINUS_INFINITY &&
                   m_point_type != LTS_VERTEX_X_PLUS_INFINITY);

    return m_x;
  }

  Point_type type() const
  {
    return m_point_type;
  }

  Number_type y() const
  {
    CGAL_assertion(m_point_type != LTS_VERTEX_Y_MINUS_INFINITY &&
                   m_point_type != LTS_VERTEX_Y_PLUS_INFINITY);

    return m_y;
  }

  void get_original_point(Point& point) const
  {

    CGAL_assertion(m_original_point_created);
    point = m_orig_point;
  }

  Rational_point_2 get_rational_point() const
  {

    CGAL_assertion(m_rat_rep_exists);
    return Rational_point_2(m_rat_x,m_rat_y);
  }
private:
  bool operator==(const Self& point_2)
  {
    switch (point_2.m_point_type)
    {
     case LTS_VERTEX_INTERIOR:
      if (this->m_point_type == LTS_VERTEX_INTERIOR &&
          this->m_x == point_2.m_x && this->m_y == point_2.m_y)
        return true;
      break;

     case LTS_VERTEX_X_MINUS_INFINITY:
      if (this->m_point_type == LTS_VERTEX_X_MINUS_INFINITY &&
          this->m_y == point_2.m_y)
        return true;
      break;

     case LTS_VERTEX_X_PLUS_INFINITY:
      if (this->m_point_type == LTS_VERTEX_X_PLUS_INFINITY &&
          this->m_y == point_2.m_y)
        return true;
      break;

     case LTS_VERTEX_Y_MINUS_INFINITY:
      if (this->m_point_type == LTS_VERTEX_Y_MINUS_INFINITY &&
          this->m_x == point_2.m_x)
        return true;
      break;

     case LTS_VERTEX_Y_PLUS_INFINITY:
      if (this->m_point_type == LTS_VERTEX_Y_PLUS_INFINITY &&
          this->m_x == point_2.m_x)
        return true;
      break;

     default:
      break;
    }
    return false;
  }

public:
  bool operator<(const Self& point_2) const
  {
    switch (point_2.m_point_type)
    {
     case LTS_VERTEX_INTERIOR:
     case LTS_VERTEX_Y_MINUS_INFINITY:
     case LTS_VERTEX_Y_PLUS_INFINITY:
      if (this->m_point_type == LTS_VERTEX_INTERIOR ||
          this->m_point_type == LTS_VERTEX_Y_MINUS_INFINITY ||
          this->m_point_type == LTS_VERTEX_Y_PLUS_INFINITY)
      {
        if (this->m_x < point_2.m_x)
          return true;

        if (this->m_x == point_2.m_x)
        {
          if (this->m_point_type == LTS_VERTEX_INTERIOR &&
              point_2.m_point_type == LTS_VERTEX_INTERIOR &&
              this->m_y < point_2.m_y)
            return true;

          if (point_2.m_point_type == LTS_VERTEX_Y_PLUS_INFINITY &&
              this->m_point_type != LTS_VERTEX_Y_PLUS_INFINITY)
            return true;

          if (point_2.m_point_type == LTS_VERTEX_INTERIOR &&
              this->m_point_type == LTS_VERTEX_Y_MINUS_INFINITY)
            return true;
        }
      }

      if (this->m_point_type == LTS_VERTEX_X_MINUS_INFINITY)
        return true;
      break;

     case LTS_VERTEX_X_MINUS_INFINITY:
      if (this->m_point_type == LTS_VERTEX_X_MINUS_INFINITY &&
          this->m_y < point_2.m_y)
        return true;
      break;

     case LTS_VERTEX_X_PLUS_INFINITY:
      if (this->m_point_type != LTS_VERTEX_X_PLUS_INFINITY ||
          this->m_y < point_2.m_y)
        return true;
      break;

     default:
      break;
    }
    return false;
  }

  virtual std::string to_string() const
  {
    std::ostringstream o;
    o << "Point (";

    if (m_point_type == LTS_VERTEX_X_MINUS_INFINITY)
      o << "MINUS INFINITY";
    else if (m_point_type == LTS_VERTEX_X_PLUS_INFINITY)
      o << "PLUS INFINITY";
    else
      o << this->m_x;

    o << ",";
    if (m_point_type == LTS_VERTEX_Y_MINUS_INFINITY)
      o << "MINUS INFINITY";
    else if (m_point_type == LTS_VERTEX_Y_PLUS_INFINITY)
      o << "PLUS INFINITY";
    else
      o << this->m_y;

    o <<  ")" << std::endl;
    return o.str();
  }

};

template <typename Traits_3_, typename Point, typename Number_type>
class Lines_through_segments_point_adapt_3 :
    public Lines_through_segments_point_adapt_2<Traits_3_, Point, Number_type>
{
private:
  typedef Lines_through_segments_point_adapt_3                  Self;
  typedef Traits_3_                                             Traits_3;

  Number_type m_z;

public:
  Lines_through_segments_point_adapt_3()
  {
  }

  Lines_through_segments_point_adapt_3(const Point& p)
  {
    this->m_orig_point = p;
    this->m_x = p.dx();
    this->m_y = p.dy();
    m_z = p.dz();
    this->m_original_point_created = true;
  }

  Lines_through_segments_point_adapt_3(const Self& p) :
    Lines_through_segments_point_adapt_2<Traits_3, Point, Number_type>(p)
  {
    m_z = p.m_z;
  }

  template <typename N_x,typename N_y,typename N_z>
  Lines_through_segments_point_adapt_3(N_x x,N_y y,N_z z)
  {
    this->m_x = x;
    this->m_y = y;
    m_z = z;
    this->m_original_point_created = false;
  }

  Number_type z() const
  {
    return this->m_z;
  }

  Number_type dx() const
  {
    return this->m_x;
  }
  Number_type dy() const
  {
    return this->m_y;
  }
  Number_type dz() const
  {
    return this->m_z;
  }

  Point get_original_point() const
  {
    if (this->m_original_point_created)
      return this->m_orig_point;
    else
      return Point(this->m_x,this->m_y,this->m_z);
  }

  void get_original_point(Point& point) const
  {
    if (this->m_original_point_created)
      point = this->m_orig_point;
    else
      point = Point(this->m_x,this->m_y,this->m_z);
  }


  std::string to_string() const
  {
    std::ostringstream o;
    if (this->m_original_point_created)
    {
      o << this->m_orig_point << std::endl;
    }
    else
    {
      o << Point(this->m_x,this->m_y,this->m_z) << std::endl;
    }

    return o.str();
  }
};


} //namespace CGAL

#endif /*LINES_THROUGH_SEGMENTS_POINT_ADAPT_H*/
