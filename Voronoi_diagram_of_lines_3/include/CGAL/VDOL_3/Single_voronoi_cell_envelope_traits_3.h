// Copyright (c) 2005-2007 Tel-Aviv University (Israel).
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
//

#ifndef CGAL_VDL3_SINGLE_CELL_ENVELOPE_TRAITS_3_H
#define CGAL_VDL3_SINGLE_CELL_ENVELOPE_TRAITS_3_H

#include <CGAL/basic.h>
#include <CGAL/Fraction_traits.h>
#include <CGAL/Lazy_exact_nt.h>

#include <CGAL/VDOL_3/vol_tools.h>
#include <CGAL/VDOL_3/svcet_3_state_dependent_functions.h>
#include <CGAL/VDOL_3/exceptions.h>
#include <CGAL/VDOL_3/Cache_container.h>
#include <CGAL/VDOL_3/Bisector.h>

#include <CGAL/VDOL_3/CGAL_SNAP_SVCET_3_TYPEDEFS.h>
#include <CGAL/VDOL_3/SVCET_3_functors/Intersect_2.h>
#include <CGAL/VDOL_3/SVCET_3_functors/Make_x_monotone_2.h>
#include <CGAL/VDOL_3/SVCET_3_functors/Make_xy_monotone_3.h>
#include <CGAL/VDOL_3/SVCET_3_functors/Construct_projected_boundary_2.h>
#include <CGAL/VDOL_3/SVCET_3_functors/Construct_projected_intersections_2.h>
#include <CGAL/VDOL_3/SVCET_3_functors/Compare_z_at_xy_3.h>
#include <CGAL/VDOL_3/SVCET_3_functors/Compare_z_at_xy_above_3.h>
#include <CGAL/VDOL_3/SVCET_3_functors/Compare_z_at_xy_below_3.h>

#include <boost/static_assert.hpp>
#include <boost/optional.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

namespace CGAL {

namespace VDOL_3
{

/*! \class Single_voronoi_cell_envelope_traits_3 
 *  \brief This function is used to compute a single cell of a voronoi 
 *         diagram of lines in space. It computes the lower envelope of
 *         the bisector surfaces of the lines and the base line over
 *         the base line.
 * \tparam T_LinearKernel Linear kernel which its Line_3 objects are used
 *                        to describe the lines.
 * \tparam T_CurvedKernel A Curved kernel class 
 *                        (like Curved_kernel_via_analysis_2). The curved
 *                        kernel should model arrangment traits concept.
 * \todo   Drop the inheritence and use the member curved kernel.
 */
template <class T_LinearKernel_2, class T_CurvedKernel_2>
class Single_voronoi_cell_envelope_traits_3 : public T_CurvedKernel_2
{
public:
  typedef T_LinearKernel_2                                Linear_kernel;
  typedef T_CurvedKernel_2                                Curved_kernel_2;
  typedef typename Curved_kernel_2::Curve_kernel_2        Algebraic_kernel_d_2;
  typedef typename Algebraic_kernel_d_2::Algebraic_kernel_d_1 Algebraic_kernel_d_1;
  typedef Single_voronoi_cell_envelope_traits_3
  < Linear_kernel, Curved_kernel_2 >                      Self;
  typedef Curved_kernel_2                                 Base;

  typedef typename Linear_kernel::FT                      FT;
  typedef typename Linear_kernel::FT                      Rational;
  
  // making sure that both kernels have the same rational number.
  // The FT of the kernel should be of the same type as the type
  // used to isolate the polynomial roots.
  BOOST_MPL_ASSERT((boost::is_same<FT,                                  \
          typename Curved_kernel_2::Curve_kernel_2::Boundary>));

  typedef typename Linear_kernel::Line_3                 Line_3;
  typedef typename Linear_kernel::Vector_3               Vector_3;
  typedef typename Linear_kernel::Point_3                Point_3;
  typedef typename Linear_kernel::Plane_3                Plane_3;

  typedef typename Base::Point_2                         Point_2;
  typedef typename Base::Curve_2                         Curve_2;
  typedef typename Base::X_monotone_curve_2              X_monotone_curve_2;
  typedef typename Base::Multiplicity                    Multiplicity;
  
  typedef typename Algebraic_kernel_d_2::Coordinate_1  Coordinate_1;
  
  // the bisector are distance functions, so if we shoot a ray from the
  // base line to any direction we get at most one intersection with the
  // quadric. This is the reason we refer to quadrics as xy_monotone_surface_3.

  typedef typename Fraction_traits<FT>::Numerator_type   Integer;
  typedef CGAL::Polynomial<Integer>                      Poly_int_1;
  typedef CGAL::Polynomial<Poly_int_1>                   Poly_int_2;
  typedef CGAL::Polynomial<Poly_int_2>                   Poly_int_3;
  typedef CGAL::Polynomial_traits_d<Poly_int_1>          PT_int_1;
  typedef CGAL::Polynomial_traits_d<Poly_int_2>          PT_int_2;
  typedef CGAL::Polynomial_traits_d<Poly_int_3>          PT_int_3;

  typedef CGAL::Polynomial<FT>                           Poly_rat_1;
  typedef CGAL::Polynomial<Poly_rat_1>                   Poly_rat_2;
  typedef CGAL::Polynomial<Poly_rat_2>                   Poly_rat_3;
  typedef CGAL::Polynomial_traits_d<Poly_rat_1>          PT_rat_1;
  typedef CGAL::Polynomial_traits_d<Poly_rat_2>          PT_rat_2;
  typedef CGAL::Polynomial_traits_d<Poly_rat_3>          PT_rat_3;

  typedef CGAL::Lazy_exact_nt<FT>                        Lazy_rat; 
  typedef CGAL::Polynomial<Lazy_rat>                     Poly_lazy_rat_1;
  typedef CGAL::Polynomial<Poly_lazy_rat_1>              Poly_lazy_rat_2;
  typedef CGAL::Polynomial<Poly_lazy_rat_2>              Poly_lazy_rat_3;
  typedef CGAL::Polynomial_traits_d<Poly_lazy_rat_1>     PT_lazy_rat_1;
  typedef CGAL::Polynomial_traits_d<Poly_lazy_rat_2>     PT_lazy_rat_2;
  typedef CGAL::Polynomial_traits_d<Poly_lazy_rat_3>     PT_lazy_rat_3;



  typedef std::pair<CGAL::Exponent_vector, FT>           Monomial;

  typedef std::list<CGAL::Object>                        Obj_list;
 
  typedef CGAL::VDOL_3::Bisector<Self> Bisector;
  typedef Bisector                                         Xy_monotone_surface_3;
  typedef Line_3                                           Surface_3;

  
  /*! The cache is for not computing all the curves again. We can have
      several $xy$-monotone of the same surface so this is quite useful.
    \todo The key of the trisector should contain also a boolean flag for 
          the positive/negative plane.
    \todo Why is it CGAL::Object? 
          Because it is the required value type of OutputIterator in 
          Construct_projected_intersections_2
    pair<xcurve, multiplicity>.
   */

  typedef CGAL::VDOL_3::Cache_container<Self> Cache_container;
  typedef typename Cache_container::Trisector_cache Trisector_cache; 

 
protected:
  
  Line_3 m_base_line;
  Linear_kernel m_linear_kernel;
  Curved_kernel_2 m_curved_kernel_2;

  // local space parameterization:
  Vector_3 m_v1;
  Vector_3 m_v2;
  Vector_3 m_v3;
  Point_3 m_origin;
  CGAL::Sign m_state;

  // pointer to common cache object
  boost::shared_ptr<Cache_container> m_shared_caches ;

public:

  Single_voronoi_cell_envelope_traits_3(){};

  Single_voronoi_cell_envelope_traits_3
  //< Linear_kernel, Curved_kernel_2 >
  (const Line_3& base_line, int seed = 0 , CGAL::Sign state_ = CGAL::POSITIVE )
    : m_base_line(base_line), m_state(state_)
  {
    m_origin =      
      m_linear_kernel.construct_point_on_3_object() (m_base_line, 0);
  
    Plane_3 pl = 
      m_linear_kernel.construct_perpendicular_plane_3_object() 
      (m_base_line, m_origin);
    
    m_v1 = m_linear_kernel.construct_vector_3_object() (m_base_line);
    m_v2 = m_linear_kernel.construct_base_vector_3_object() (pl, 1);
    m_v3 = m_linear_kernel.construct_base_vector_3_object() (pl, 2);
    
    m_v2 = m_v2 + seed * m_v3;
    m_v3 = m_linear_kernel.construct_cross_product_vector_3_object()(m_v1,m_v2);
    
    if(m_state == CGAL::NEGATIVE){
      m_v2 = -m_v2;
      m_v3 = -m_v3;
    }
     
    m_shared_caches = boost::make_shared<Cache_container>();
  }
  
  const Line_3& base_line() const { return m_base_line; }  
  const Line_3& base_line_3() const { return m_base_line; }  
  const Linear_kernel& linear_kernel() const { return m_linear_kernel; }  
  const Curved_kernel_2& curved_kernel_2() const { return m_curved_kernel_2; }  
  const Vector_3& v1() const { return m_v1; }  
  const Vector_3& v2() const { return m_v2; }  
  const Vector_3& v3() const { return m_v3; }  
  const Point_3& origin() const { return m_origin; }  
  CGAL::Sign state() const { return m_state; }  
  boost::shared_ptr<Cache_container> caches_ptr() const { return m_shared_caches ; }  
  
  Trisector_cache& trisector_cache() const 
  { return caches_ptr()->m_trisector_cache; }

  Self mirror() const{
    Self result = *this; 
    result.m_state = -(result.m_state);
    result.m_v2    = -(result.m_v2);
    result.m_v3    = -(result.m_v3);
    return result;
  }
  
#define SVCET_3_snap_functor(_Name,_get_object)         \
  typedef SVCET_3_functors::_Name<Self> _Name;          \
    _Name _get_object() const  { return _Name(this); }  \
    
  SVCET_3_snap_functor(
      Intersect_2,
      intersect_2_object);  

  SVCET_3_snap_functor(
      Make_x_monotone_2,
      make_x_monotone_2_object);

  SVCET_3_snap_functor(
      Make_xy_monotone_3,
      make_xy_monotone_3_object);
  SVCET_3_snap_functor(
      Construct_projected_boundary_2,
      construct_projected_boundary_2_object);
  SVCET_3_snap_functor(
      Construct_projected_intersections_2,
      construct_projected_intersections_2_object);
  SVCET_3_snap_functor(
      Compare_z_at_xy_3,
      compare_z_at_xy_3_object);
  SVCET_3_snap_functor(
      Compare_z_at_xy_above_3,
      compare_z_at_xy_above_3_object);
  SVCET_3_snap_functor(
      Compare_z_at_xy_below_3,
      compare_z_at_xy_below_3_object);
#undef SVCET_3_snap_functor

public: 
  const Poly_int_2& projection(const Bisector& s1, const Bisector& s2) const {
    CGAL_assertion(s1.bisector() == CGAL::canonicalize(s1.bisector()));
    CGAL_assertion(s2.bisector() == CGAL::canonicalize(s2.bisector()));
    
    if(s1.bisector() < s2.bisector()) return projection(s2,s1);

    // this is a boost::optional<Poly_int_2> 
    typedef typename Cache_container::Projection_cache_data Projection_cache_data;
    Projection_cache_data& projection = 
      (caches_ptr()->m_projection_cache)[std::make_pair(s1.bisector(),s2.bisector())];
    if(!projection.is_initialized()){
      projection = Projection_cache_data(
          CGAL::canonicalize(
              CGAL::resultant(
                  s1.transformed_bisector(), s2.transformed_bisector())));
      //also check for generic projection 
      if (total_degree(projection) < 
          2*total_degree(s1.bisector())*total_degree(s2.bisector())){
        // compute gcd if degree of projection seems too low
        Poly_int_2 g = CGAL::gcd(
            VDOL_3::intersect_with_hstar(this,s1.bisector()),
            VDOL_3::intersect_with_hstar(this,s2.bisector()));
        if(CGAL::total_degree(g)!=0){
          throw Event_at_discontinuity_exception(
              " 1-dimensional component of trisector in H^* ");
        }
      }

    }else
    CGAL_assertion(projection.is_initialized());
    CGAL_assertion( projection.get() == 
        CGAL::canonicalize(
            CGAL::resultant(
                s1.transformed_bisector(), s2.transformed_bisector()))); 
    
    return projection.get(); 
  }
  
  template<class OutputIterator>
  OutputIterator projection(const Bisector& s1, const Bisector& s2, OutputIterator oit) const {
    CGAL_assertion(s1.bisector() == CGAL::canonicalize(s1.bisector()));
    CGAL_assertion(s2.bisector() == CGAL::canonicalize(s2.bisector()));
    if(s1.bisector() < s2.bisector()) return projection(s2,s1,oit);
    
    // get a reference to (maybe empty) data 
    typename Cache_container::Projection_sqff_cache_data& data = 
      (caches_ptr()->m_projection_sqff_cache) [ std::make_pair(s1.bisector(),s2.bisector()) ];
    
    // std::cout << data.is_initialized() << std::flush; 
    if(!data.is_initialized()){
      data = std::list<std::pair<Poly_int_2,int> >();
      this->kernel().square_free_factorize_2_object()(projection(s1,s2),std::back_inserter(*data));
    }
    CGAL_postcondition(data.is_initialized());
    std::copy(data->begin(),data->end(),oit);
    return oit;
  }
  
  const Poly_int_2& swap_2(const Poly_int_2& f) const {
    boost::optional<Poly_int_2>& data = 
      caches_ptr()->m_swap_2_cache[f]; 
    if(!data.is_initialized()){data= CGAL::swap(f,0,1);}
    CGAL_postcondition(data.is_initialized());
    return *data;
  }

  template <typename OutputIterator>
  OutputIterator solve_1(const Poly_int_1& f, OutputIterator oit) const{
    CGAL_precondition(CGAL::is_square_free(f));
    typename Cache_container::Solve_1_cache_data& data =
      this->caches_ptr()->m_solve_1_cache[f];
    // std::cout << data.is_initialized() << std::flush; 
    if(!data.is_initialized()){
      data = std::list<Coordinate_1>(); 
      this->kernel().solve_1_object()(f,std::back_inserter(*data),false);
    }
    CGAL_postcondition(data.is_initialized());
    return std::copy(data->begin(),data->end(),oit);     
  }

  // PRIVATE MEMBER FUNCTION 
  template <typename InputIterator, typename OutputIterator>
  OutputIterator 
  split_curves(InputIterator begin, InputIterator end, OutputIterator oit) const {
    // split curves into non-intersection x-monotone curves 
    return CGAL::compute_subcurves (begin,end, oit, false, *this);
  }
  
};

} // end of namespace VDOL_3

} //namespace CGAL

#endif // CGAL_VDL3_SINGLE_CELL_ENVELOPE_TRAITS_3_H
