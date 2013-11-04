#ifndef CGAL_DUMMY_BEZIER_TRAITS_H_
#define CGAL_DUMMY_BEZIER_TRAITS_H_


#include <CGAL/tags.h>
#include <CGAL/Arr_tags.h>

//#include "Bezier_curve_2.h"
//#include "Bezier_point_2.h"
//#include "Bezier_x_monotone_2.h"
//#include <CGAL/Arr_geometry_traits/Bezier_curve_2.h>
//#include <CGAL/Arr_geometry_traits/Bezier_point_2.h>
//#include <CGAL/Arr_geometry_traits/Bezier_x_monotone_2.h>

#include <CGAL/Arr_geometry_traits/Bezier_bounding_rational_traits.h>

// New Header files added
#include <CGAL/Algebraic_kernel_d_1.h>
#include "../Bezier_traits_cache.h"



namespace CGAL {

/*! \class
 * A traits class for maintaining an arrangement of Bezier curves with
 * rational control points.
 *
 * The class is templated with four parameters: 
 * Rat_kernel A kernel that defines the type of control points.
 * Alg_kernel A geometric kernel, where Alg_kernel::FT is the number type
 *            for the coordinates of arrangement vertices and is used to
 *            represent algebraic numbers.
 * Nt_traits A number-type traits class. This class defines the Rational
 *           number type (should be the same as Rat_kernel::FT) and the
 *           Algebraic number type (should be the same as Alg_kernel::FT)
 *           and supports various operations on them.
 * Bounding_traits A traits class for filtering the exact computations.
 *                 By default we use the rational bounding traits.
 */

template <class AlgebraicKernel_d_1_>
//, class RatKernel_, 
//				class AlgKernel_, class NtTraits_,
//				class BoundingTraits_ = Bezier_bounding_rational_traits<RatKernel_> >
class Arr_Bezier_curve_traits_2 
{
public:


  // Algebraic kernel to remove Nt_traits and Alg_kernel
  typedef AlgebraicKernel_d_1_			 Algebraic_kernel_d_1;

/*
  typedef RatKernel_                             Rat_kernel;
  typedef AlgKernel_                             Alg_kernel;
  typedef NtTraits_                              Nt_traits;
  typedef BoundingTraits_                        Bounding_traits;

  typedef Arr_Bezier_curve_traits_2<Algebraic_kernel_d_1,
				    Rat_kernel,
                                    Alg_kernel,
                                    Nt_traits,
                                    Bounding_traits>   Self;
 
  typedef typename Nt_traits::Integer            Integer;
  typedef typename Rat_kernel::FT                Rational;
  typedef typename Alg_kernel::FT                Algebraic;

  typedef typename Rat_kernel::Point_2           Rat_point_2;
  typedef typename Alg_kernel::Point_2           Alg_point_2;
  
  // Category tags:
  typedef Tag_true                               Has_left_category;
  typedef Tag_true                               Has_merge_category;
  typedef Tag_false                              Has_do_intersect_category;

  typedef Arr_oblivious_side_tag                 Left_side_category;
  typedef Arr_oblivious_side_tag                 Bottom_side_category;
  typedef Arr_oblivious_side_tag                 Top_side_category;
  typedef Arr_oblivious_side_tag                 Right_side_category;

  // Traits-class types:
  typedef _Bezier_curve_2<Rat_kernel,
                          Alg_kernel,
                          Nt_traits,
                          Bounding_traits>             Curve_2;

  typedef _Bezier_x_monotone_2<Rat_kernel,
                               Alg_kernel,
                               Nt_traits,
                               Bounding_traits>        X_monotone_curve_2;

  typedef _Bezier_point_2<Rat_kernel,
                          Alg_kernel,
                          Nt_traits,
                          Bounding_traits>             Point_2;

  typedef typename X_monotone_curve_2::Multiplicity    Multiplicity;

*/
  // Type definition for the vertical-tangnecy and intersection point cache.
  // Changed to use the algebraic kernel
  typedef _Bezier_cache<Algebraic_kernel_d_1>		Bezier_cache;


private:

  // Type definition for the bounded intersection points mapping.
//  typedef typename X_monotone_curve_2::Intersection_map   Intersection_map;

  // Data members:
  mutable Bezier_cache * p_cache;         /*!< Caches vertical tangency points
                                           * and intersection points that have
                                           * been computed in an exact manner.
                                           */
//  mutable Intersection_map * p_inter_map; /*!< Maps curve pairs to their
//                                           * intersection points.
//                                           */
//  bool m_owner;                           /*!< Does this instance own its cache
 //                                          * and map structures.
//                                           */

public:

  /// \name Construction.
  //@{

  /*! Default constructor. */
  Arr_Bezier_curve_traits_2 ()
  {
    p_cache = new Bezier_cache;
//    p_inter_map = new Intersection_map;
//    m_owner = true;
  }


  //@}
};

} //namespace CGAL


#endif
