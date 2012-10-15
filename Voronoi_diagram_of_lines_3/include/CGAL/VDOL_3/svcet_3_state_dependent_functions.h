// Copyright (c) 2005-2009 Tel-Aviv University (Israel).
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
// $Id:
// 
//
// Author(s): Ophir Setter          <ophirset@post.tau.ac.il>
//            Michael Hemmer        <hemmer@mpi-inf.mpg.de>
//

// SVCET_3 = Single_voronoi_cell_envelope_traits_3

// TODO: Change use of svcet_3 pointer to const& 

#ifndef CGAL_VDL3_SVCET_33_STATE_DEPENDENT_FUNCTIONS
#define CGAL_VDL3_SVCET_33_STATE_DEPENDENT_FUNCTIONS

#include <CGAL/config.h>
#include <CGAL/VDOL_3/vol_tools.h>
#include <CGAL/VDOL_3/CGAL_SNAP_SVCET_3_TYPEDEFS.h>
#include <CGAL/VDOL_3/vol_tools.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Root_of_traits.h>

namespace CGAL { 
namespace VDOL_3 {

template <typename SVCET_3>
typename SVCET_3::Poly_int_3
construct_bisector_3(
    const SVCET_3 * svcet_3,
    const typename SVCET_3::Line_3& l1){

  CGAL_SNAP_SVCET_3_TYPEDEFS;
  
  Poly_rat_3 bisector_rat = 
    construct_bisector_3<Linear_kernel>(svcet_3->base_line(),l1);

  typedef typename CGAL::Fraction_traits<Poly_rat_3> FRT;
  typename FRT::Decompose decompose;
  typename FRT::Denominator_type dummy;  
  typename SVCET_3::Poly_int_3 bisector_int;
  
  decompose(bisector_rat,bisector_int,dummy);
  return CGAL::canonicalize(bisector_int);
}


template <class SVCET_3>
typename SVCET_3::Poly_int_3 
construct_transformed_bisector_3(
    const SVCET_3 *svcet_3,
    const typename SVCET_3::Line_3& line){ 
  CGAL_SNAP_SVCET_3_TYPEDEFS;
  Poly_int_3 bisector = construct_bisector_3(svcet_3,line); 
  return construct_transformed_bisector_3(svcet_3,bisector);
}


template <class SVCET_3>
typename SVCET_3::Poly_int_3 
construct_transformed_bisector_3(
    const SVCET_3 *svcet_3,
    const typename SVCET_3::Poly_int_3& bisecetor){
  
  CGAL_SNAP_SVCET_3_TYPEDEFS;
  
  Poly_rat_3 rat_bisector = 
    typename Coercion_traits<Poly_rat_3,Poly_int_3>::Cast()(bisecetor);

  typename  PT_rat_3::Shift shift_rat_3;
  Poly_rat_3 u = shift_rat_3(Poly_rat_3(1),1,0); // u^1
  Poly_rat_3 v = shift_rat_3(Poly_rat_3(1),1,1); // v^1
  Poly_rat_3 r = shift_rat_3(Poly_rat_3(1),1,2); // r^1
  Poly_rat_3 vr = v*r; 
  
  typename PT_rat_3::Construct_polynomial construct_polynomial;

  typename Linear_kernel::Compute_x_3 compute_x_3 = 
    svcet_3->linear_kernel().compute_x_3_object();
  typename Linear_kernel::Compute_y_3 compute_y_3 = 
    svcet_3->linear_kernel().compute_y_3_object();
  typename Linear_kernel::Compute_z_3 compute_z_3 = 
    svcet_3->linear_kernel().compute_z_3_object();
  
  std::vector<Poly_rat_3> subtitution_vec;
  FT v1_x = compute_x_3(svcet_3->v1());
  FT v2_x = compute_x_3(svcet_3->v2());
  FT v3_x = compute_x_3(svcet_3->v3());
  FT origin_x = compute_x_3(svcet_3->origin());
  Poly_rat_3 px = construct_polynomial(v1_x) * u + 
    construct_polynomial(v2_x) * vr + 
    construct_polynomial(v3_x) * r + construct_polynomial(origin_x);
  FT v1_y = compute_y_3(svcet_3->v1());
  FT v2_y = compute_y_3(svcet_3->v2());
  FT v3_y = compute_y_3(svcet_3->v3());
  FT origin_y = compute_y_3(svcet_3->origin());
  Poly_rat_3 py = construct_polynomial(v1_y) * u + 
    construct_polynomial(v2_y) * vr + 
    construct_polynomial(v3_y) * r + construct_polynomial(origin_y);
  FT v1_z = compute_z_3(svcet_3->v1());
  FT v2_z = compute_z_3(svcet_3->v2());
  FT v3_z = compute_z_3(svcet_3->v3());
  FT origin_z = compute_z_3(svcet_3->origin());
  Poly_rat_3 pz = construct_polynomial(v1_z) * u + 
    construct_polynomial(v2_z) * vr + 
    construct_polynomial(v3_z) * r + construct_polynomial(origin_z);
  
  subtitution_vec.push_back(px);
  subtitution_vec.push_back(py);
  subtitution_vec.push_back(pz);

  Poly_rat_3 rat_transformed_bisector = typename PT_rat_3::Substitute()
    (rat_bisector, subtitution_vec.begin(), subtitution_vec.end());

  typedef typename CGAL::Fraction_traits<Poly_rat_3> FRT;
  typename FRT::Decompose decompose;
  typename FRT::Denominator_type dummy;  
  Poly_int_3 transformed_bisector;
  decompose(rat_transformed_bisector,transformed_bisector,dummy);
  return CGAL::canonicalize(transformed_bisector);
}


template <class SVCET_3>
typename SVCET_3::Poly_int_2 
intersect_with_hstar(
    const SVCET_3 *svcet_3,
    const typename SVCET_3::Poly_int_3& bisecetor){
  
  CGAL_SNAP_SVCET_3_TYPEDEFS;
  
  Poly_rat_3 rat_bisector = 
    typename Coercion_traits<Poly_rat_3,Poly_int_3>::Cast()(bisecetor);

  // typename  PT_rat_2::Shift shift_rat_2;
  Poly_rat_2 u = CGAL::shift(Poly_rat_2(1),1,0); // u^1
  Poly_rat_2 v = CGAL::shift(Poly_rat_2(1),1,1); // v^1
  
  typename PT_rat_2::Construct_polynomial construct_polynomial;

  typename Linear_kernel::Compute_x_3 compute_x_3 = 
    svcet_3->linear_kernel().compute_x_3_object();
  typename Linear_kernel::Compute_y_3 compute_y_3 = 
    svcet_3->linear_kernel().compute_y_3_object();
  typename Linear_kernel::Compute_z_3 compute_z_3 = 
    svcet_3->linear_kernel().compute_z_3_object();
  
  std::vector<Poly_rat_2> subtitution_vec;
  FT v1_x = compute_x_3(svcet_3->v1());
  FT v2_x = compute_x_3(svcet_3->v2());
  FT origin_x = compute_x_3(svcet_3->origin());
  Poly_rat_2 px = construct_polynomial(v1_x) * u + 
    construct_polynomial(v2_x) * v + construct_polynomial(origin_x);
  FT v1_y = compute_y_3(svcet_3->v1());
  FT v2_y = compute_y_3(svcet_3->v2());
  FT origin_y = compute_y_3(svcet_3->origin());
  Poly_rat_2 py = construct_polynomial(v1_y) * u + 
    construct_polynomial(v2_y) * v +  construct_polynomial(origin_y);
  FT v1_z = compute_z_3(svcet_3->v1());
  FT v2_z = compute_z_3(svcet_3->v2());
  FT origin_z = compute_z_3(svcet_3->origin());
  Poly_rat_2 pz = construct_polynomial(v1_z) * u + 
    construct_polynomial(v2_z) * v +  construct_polynomial(origin_z);
  
  subtitution_vec.push_back(px);
  subtitution_vec.push_back(py);
  subtitution_vec.push_back(pz);

  Poly_rat_2 rat_intersection = typename PT_rat_3::Substitute()
    (rat_bisector, subtitution_vec.begin(), subtitution_vec.end());

  typedef typename CGAL::Fraction_traits<Poly_rat_2> FRT;
  typename FRT::Decompose decompose;
  typename FRT::Denominator_type dummy;  
  Poly_int_2 intersection;
  decompose(rat_intersection,intersection,dummy);
  return CGAL::canonicalize(intersection);
}

template <class SVCET_3>
typename SVCET_3::Poly_int_1 shoot_ray_at(
    const SVCET_3 *svcet_3, 
    const typename SVCET_3::Line_3& l1,
    const std::pair<typename SVCET_3::FT, typename SVCET_3::FT>& point){
  CGAL_SNAP_SVCET_3_TYPEDEFS;

  Poly_int_3 s1 = construct_transformed_bisector_3(svcet_3,l1);
  
  std::vector<Poly_rat_1> rat_polys;
  typename PT_rat_1::Construct_polynomial construct_poly;
  rat_polys.push_back(construct_poly(point.first));
  rat_polys.push_back(construct_poly(point.second));
  rat_polys.push_back(construct_poly(FT(0), FT(1)));

  typename PT_int_3::Substitute substitue;
  Poly_rat_1 rat_poly1 = substitue(s1, rat_polys.begin(), rat_polys.end());
  Poly_int_1 poly1;
  Integer dummy_int;
  typename Fraction_traits<Poly_rat_1>::Decompose() 
    (rat_poly1, poly1, dummy_int);
  return CGAL::canonicalize(poly1);
}

#if 1
template <class SVCET_3>
static Comparison_result compare_bisectors_at(
    const SVCET_3* svcet_3,
    const typename SVCET_3::Xy_monotone_surface_3& s1,
    const typename SVCET_3::Xy_monotone_surface_3& s2,
    const std::pair<typename SVCET_3::FT, typename SVCET_3::FT>& point)
{
  typedef typename SVCET_3::Poly_lazy_rat_3 Poly_lazy_rat_3; 
  typedef typename SVCET_3::Poly_int_3 Poly_int_3; 
  typedef typename SVCET_3::FT FT;
  typedef CGAL::Interval_nt<true> Interval; 
  typename CGAL::Coercion_traits<Poly_lazy_rat_3,Interval>::Cast cast; 
  try
    {
      return 
        compare_transformed_bisectors_at(
            cast(s1.lazy_rat_transformed_bisector()), 
            cast(s2.lazy_rat_transformed_bisector()), 
            std::make_pair(
                Interval(CGAL::to_interval(point.first)),
                Interval(CGAL::to_interval(point.second))));
    }
  catch (...)
    {      
      typename CGAL::Coercion_traits<Poly_int_3,FT>::Cast cast; 
      return 
        compare_transformed_bisectors_at(
            cast(s1.transformed_bisector()), 
            cast(s2.transformed_bisector()), 
            point); 
      
    }
}


template <class Polynomial_3, class FT>
static Comparison_result compare_transformed_bisectors_at(
    const Polynomial_3& f,
    const Polynomial_3& g,
    const std::pair<FT, FT>& point)
{
  // assuming that FT is Innermost_coefficient
  
  typedef CGAL::Polynomial_traits_d<Polynomial_3>      PT_3; 
  typedef typename PT_3:: template Rebind<FT,1>::Other PT_1;
  typedef typename PT_1::Polynomial_d Poly_1; 

  // prepare for substitue of  x and y 
  std::vector<Poly_1> polys;
  typename PT_1::Construct_polynomial construct_npoly_1;
  polys.push_back(construct_npoly_1((point.first)));
  polys.push_back(construct_npoly_1((point.second)));
  polys.push_back(construct_npoly_1(FT(0), FT(1)));

  // substitue of  x and y 
  typename PT_3::Substitute substitue;
  Poly_1 f_1 = substitue(f, polys.begin(), polys.end());
  Poly_1 g_1 = substitue(g, polys.begin(), polys.end());
  CGAL_postcondition(CGAL::degree(f_1) <= 2);
  CGAL_postcondition(CGAL::degree(g_1) <= 2);
  
  typename CGAL::Root_of_traits<FT>::Root_of_2 rf(
    make_root_of_2(
        CGAL::get_coefficient(f_1,2),
        CGAL::get_coefficient(f_1,1),
        CGAL::get_coefficient(f_1,0),
        false)); // get the larger root

  typename CGAL::Root_of_traits<FT>::Root_of_2 rg(
  make_root_of_2(
        CGAL::get_coefficient(g_1,2),
        CGAL::get_coefficient(g_1,1),
        CGAL::get_coefficient(g_1,0),
        false)); // get the larger root  

  // of there is no positive root the surface is at infinity  
  CGAL::Sign srf = CGAL::sign(rf);
  CGAL::Sign srg = CGAL::sign(rg);
  if(srf == CGAL::NEGATIVE && srg == CGAL::NEGATIVE) return CGAL::EQUAL; 
  if(srf == CGAL::NEGATIVE) return CGAL::LARGER;  // f is at infinity 
  if(srg == CGAL::NEGATIVE) return CGAL::SMALLER; // g is at infinity 
  CGAL::Comparison_result result = CGAL::compare(rf,rg); 
  return result; 
  
}

#else 
template <class SVCET_3>
static Comparison_result compare_bisectors_at(
    const SVCET_3* svcet_3,
    const typename SVCET_3::Xy_monotone_surface_3& s1,
    const typename SVCET_3::Xy_monotone_surface_3& s2,
    const std::pair<typename SVCET_3::FT, typename SVCET_3::FT>& point)
{
  return compare_transformed_bisectors_at(
      svcet_3, s1.transformed_bisector(), s2.transformed_bisector(),point);


}

template <class SVCET_3>
static Comparison_result compare_transformed_bisectors_at(
    const SVCET_3* svcet_3,
    const typename SVCET_3::Poly_int_3& s1,
    const typename SVCET_3::Poly_int_3& s2,
    const std::pair<typename SVCET_3::FT, typename SVCET_3::FT>& point)
{
  // TODO use STATE 
  CGAL_SNAP_SVCET_3_TYPEDEFS;
  
  std::vector<Poly_rat_1> rat_polys;
  typename PT_rat_1::Construct_polynomial construct_poly;
  rat_polys.push_back(construct_poly(point.first));
  rat_polys.push_back(construct_poly(point.second));
  rat_polys.push_back(construct_poly(FT(0), FT(1)));
  
  typename PT_int_3::Substitute substitue;
  Poly_rat_1 rat_poly1 = substitue(s1, rat_polys.begin(), rat_polys.end());
  Poly_rat_1 rat_poly2 = substitue(s2, rat_polys.begin(), rat_polys.end());
  Poly_int_1 poly1, poly2;
  Integer dummy_int;
  typename Fraction_traits<Poly_rat_1>::Decompose() 
    (rat_poly1, poly1, dummy_int);
  typename Fraction_traits<Poly_rat_1>::Decompose() 
    (rat_poly2, poly2, dummy_int);
  
  Algebraic_kernel_d_1 ak; // should be obtained from traits 
  return CGAL::VDOL_3::compare_smallest_nonnegative_roots(ak, poly1, poly2);
}
#endif 

template <class SVCET_3>
static Comparison_result compare_lines_at(
    const SVCET_3* svcet_3,
    const typename SVCET_3::Line_3& l1,
    const typename SVCET_3::Line_3& l2,
    const std::pair<typename SVCET_3::FT, typename SVCET_3::FT>& point)
{
  CGAL_SNAP_SVCET_3_TYPEDEFS;
  Poly_int_1 poly_1 = shoot_ray_at(svcet_3,l1,point);
  Poly_int_1 poly_2 = shoot_ray_at(svcet_3,l2,point);
  // std::cout << poly_1 << std::endl;
  // std::cout << poly_2 << std::endl;
  Algebraic_kernel_d_2 ak_2 = svcet_3->kernel();
  CGAL::Comparison_result result = 
    CGAL::VDOL_3::compare_smallest_nonnegative_roots(ak_2, poly_1, poly_2);
  // std::cout << result << std::endl;
  return result; 
}

template <class SVCET_3>
static bool is_at_infinity(
    const SVCET_3* svcet_3,
    const typename SVCET_3::Line_3& l1,
    const std::pair<typename SVCET_3::FT, typename SVCET_3::FT>& point)
{
  CGAL_SNAP_SVCET_3_TYPEDEFS;
  Poly_int_1 poly_1 = shoot_ray_at(svcet_3, l1, point);
  
  // In this case the surface contains the ray
  if(CGAL::is_zero(poly_1)) return false;
  
  // case of intersecting lines
  // the ray is parallel to the singular line of the intersecting planes
  // hence it hits no surface in both directions 
  if(CGAL::degree(poly_1)==0) 
    return true; 
  // the ray hits the surface twice, i.e. no root at infinity
  if(CGAL::degree(poly_1)==2) return false; 
  // the ray hits the surface only once (+ at infinity)
  if( CGAL::sign(CGAL::get_coefficient(poly_1,0))*
      CGAL::sign(CGAL::get_coefficient(poly_1,1))== 
      CGAL::NEGATIVE){
    return false; // a positive root 
  }else{
    return true; // since the existing root is positive
  }
}

template <class SVCET_3> 
bool is_singular_point_of(
    const SVCET_3* svcet_3, 
    const typename SVCET_3::Point_2& p, 
    const typename SVCET_3::Xy_monotone_surface_3& surface){
  return  
    surface.boundary().size()==2 && 
    svcet_3->is_on_2_object()(p,surface.boundary().front().first) && 
    svcet_3->is_on_2_object()(p,surface.boundary().back().first);
}
 

template <class SVCET_3>
bool is_undefined_at(
    const SVCET_3* svcet_3, 
    const typename SVCET_3::Point_2& p, 
    const typename SVCET_3::Xy_monotone_surface_3& surface){
  
  if(CGAL::total_degree(surface.bisector()) == 2) 
    return false; 
  
  CGAL_assertion(CGAL::total_degree(surface.bisector()) == 1 );
  CGAL_assertion(surface.boundary().size() == 1 );
  
  //if(this->svcet_3()->is_on_2_object()(p,surface.boundary().front().first))
  //  return false;
    
  typename SVCET_3::Compare_y_at_x_2
    compare_y_at_x_2(svcet_3->compare_y_at_x_2_object());
  
  if( compare_y_at_x_2(p,surface.boundary().front().first) == - surface.boundary().front().second )
    return true; 
  return false;  // on boundary or inside 
}

template <class SVCET_3>
bool has_infinit_distance_at(
    const SVCET_3* svcet_3, 
    const typename SVCET_3::Point_2& p, 
    const typename SVCET_3::Xy_monotone_surface_3& surface){
  CGAL_precondition(!is_undefined_at(svcet_3,p,surface));
  return 
    surface.boundary().size() >= 1 && 
    svcet_3->is_on_2_object()(p,surface.boundary().front().first); 
}

template <class SVCET_3>
bool has_zero_distance_at(
    const SVCET_3* svcet_3, 
    const typename SVCET_3::Point_2& p, 
    const typename SVCET_3::Xy_monotone_surface_3& surface){
  CGAL_precondition(!is_undefined_at(svcet_3,p,surface));
  return 
    surface.boundary().size()==2 &&  
    svcet_3->is_on_2_object()(p,surface.boundary().back().first);
}

template <class SVCET_3>
bool is_singular_line_at(
    const SVCET_3* svcet_3, 
    const typename SVCET_3::Point_2& p, 
    const typename SVCET_3::Xy_monotone_surface_3& surface){
  CGAL_precondition(!is_undefined_at(svcet_3,p,surface));
  return  
    surface.boundary().size()==2 && 
    svcet_3->is_on_2_object()(p,surface.boundary().front().first) && 
    svcet_3->is_on_2_object()(p,surface.boundary().back().first);
}


} // namespace VDOL_3 
} // namespace CGAL

#endif // CGAL_VDL3_SVCET_33_STATE_DEPENDENT_FUNCTIONS
