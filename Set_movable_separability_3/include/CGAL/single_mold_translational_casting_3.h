// Copyright (c) 2005,2006,2007,2008,2009,2010,2011 Tel-Aviv University (Israel).
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
// Author(s): Shahar    <shasha94@gmail.com>
//            Efi Fogel <efif@gmail.com>

#ifndef CGAL_SINGLE_MOLD_TRANSLATIONAL_CASTING_3_H
#define CGAL_SINGLE_MOLD_TRANSLATIONAL_CASTING_3_H

#include <iostream>
#include <list>
#include <boost/type_traits/is_same.hpp>

#include <CGAL/Kernel_traits.h>
#include <CGAL/enum.h>
#include <CGAL/value_type_traits.h>
#include <CGAL/property_map.h>

#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>

namespace CGAL {

enum geom_traits_t { geom_traits };
enum all_default_t { all_default };

namespace Set_movable_separability_3 {

template<typename Polyhedron, typename NamedParameters>
class Get_kernel {
  typedef typename boost::property_map<Polyhedron,
                                       boost::vertex_point_t>::const_type
    Propert_map;
  typedef typename boost::property_traits<Propert_map>::value_type      Point;

  typedef typename CGAL::Kernel_traits<Point>::Kernel
    Default_kernel;

public:
  typedef typename boost::lookup_named_param_def<CGAL::geom_traits_t,
                                                 NamedParameters,
                                                 Default_kernel>::type
    type;
};

template <typename Polyhedron,  typename NamedParameters,
          typename OutputIterator>
OutputIterator
single_mold_translational_casting_3_impl(const Polyhedron& polyhedron,
                                         const NamedParameters& np,
                                         OutputIterator oi, boost::false_type)
{
  return oi;
}

template <typename Polyhedron,  typename NamedParameters,
          typename OutputIterator>
OutputIterator
single_mold_translational_casting_3_impl(const Polyhedron& polyhedron,
                                         const NamedParameters& np,
                                         OutputIterator oi, boost::true_type)
{
  return oi;
}

/*! \fn OutputIterator find_single_mold_translational_casting_3(const Polyhedron& polyhedron, OutputIterator oi)
 * \param[in] polyhedron the input polyhedron.
 * \param[out] oi the output iterator. Its value type is a pair, where
 *             (i) the first element in the pair identifies a valid top face
 *                 represented by its index the type of which is convertible to
                   `boost::graph_traits<Polyhedron>::face_descriptor`, and
 *             (ii) the second element is a closed spherical patch of pull-out
 *                  3D directions represented as a sequence of the extreme
 *                  directions in the patch of type `Kernel::Direction_3`.
 * \return the past-the-end iterator of the output container.
 * \pre `polyhedron` must be non-degenerate (has at least 4 vertices and 6
 *      edges), simple, and does not have neighboring coplanar facets.
 */
template <typename Polyhedron,  typename NamedParameters,
          typename OutputIterator, typename DirectionType>
OutputIterator
single_mold_translational_casting_3(const Polyhedron& polyhedron,
                                    const NamedParameters& np,
                                    OutputIterator oi)
{
  typedef typename Get_kernel<Polyhedron, NamedParameters>::type Kernel;
  typedef typename Kernel::Direction_3                           Direction_3;

  //! \todo consider using CGAL::is_same_or_derived instaed of boost::is_same
  typedef typename boost::is_same<DirectionType, Direction_3>::type
    Is_direction;

  return single_mold_translational_casting_3_impl(polyhedron, np, oi,
                                                  Is_direction());
}
#define PRINT_PLANE(p) do{ if((p).a()!=0)\
std::cout<<(p).a()<<'x'<<(((p).b())>0?"+":"");\
if((p).b()!=0)\
std::cout<<(p).b()<<'y'<<(((p).c()>0)?"+":"");\
if((p).c()!=0)\
std::cout<<(p).c()<<'z';\
std::cout<<'='<<(p).d()<<std::endl;}while(0)

// program and solution types
typedef CGAL::Linear_program_from_iterators
<int**,                                                // for A
 int*,                                                 // for b
 CGAL::Const_oneset_iterator<CGAL::Comparison_result>, // for r
 bool*,                                                // for fl
 int*,                                                 // for l
 bool*,                                                // for fu
 int*,                                                 // for u
 int*>                                                 // for c
Program;
template <typename Polyhedron,  typename NamedParameters,
          typename OutputIterator>
OutputIterator
single_mold_translational_casting_3(const Polyhedron& polyhedron,
                                    const NamedParameters& np,
                                    OutputIterator oi)
{
  typedef typename value_type_traits<OutputIterator>::type       Value_type;
  typedef typename Value_type::second_type                       Direction_type;

  typedef typename Get_kernel<Polyhedron, NamedParameters>::type Kernel;
  typedef typename Kernel::Direction_3                           Direction_3;

  //! \todo consider using CGAL::is_same_or_derived instaed of boost::is_same
  typedef typename boost::is_same<Direction_type, Direction_3>::type
    Is_direction;

	typedef typename Kernel::FT                       FT;

	typedef CGAL::Quadratic_program_solution<FT> Solution;
	for(auto it = polyhedron.planes_begin(); it!=polyhedron.planes_end();it++)
	{
	    PRINT_PLANE(*it);
		//std::cout<<it->a()<<'x'<<((it->b()>0)?"+":"")<<it->b()<<it->c()<<it->d()<<std::endl;

//		const CGAL::Plane_3<CGAL::Epick> plane = *it;
//		Polyhedron::PolyhedronTraits_3 b;
//		int a=b;
//		std::cout<<*it<<std::endl;
	}
	  int  Ax[] = {1, -1 ,0, 0};                        // column for x

			  int  Ay[] = {1, 1 , -1, -1};                        // column for y
			  int*  A[] = {Ax, Ay};                       // A comes columnwise
			  int   b[] = {0, 0,-1 , -7 };                         // right-hand side
/*
 * x+y<=0
 * -x+y<=0
 * -y<=1
 *
 */

			  CGAL::Const_oneset_iterator<CGAL::Comparison_result>
			        r(    CGAL::SMALLER);                 // constraints are "<="
			  bool fl[] = {false, false};                   // both x, y are lower-bounded
			  int   l[] = {0, 0};
			  bool fu[] = {false, false};                  // only y is upper-bounded
			  int   u[] = {0, 4};                         // x's u-entry is ignored
			  int   c[] = {0, -32};
			  int  c0   = 64;                             // constant term
			  // now construct the linear program; the first two parameters are
			  // the number of variables and the number of constraints (rows of A)
			  Program lp (2, 4, A, b, r, fl, l, fu, u, c, c0);
			  // solve the program, using ET as the exact type
			  Solution s = CGAL::solve_linear_program(lp, FT());
//			  CGAL::Quadratic_program_solution<double> a;
			  if(s.is_infeasible())
			    {
				for(auto it= s.infeasibility_certificate_begin(); it!=s.infeasibility_certificate_end();it++)
				  {
				    std::cout<<*it<<std::endl;
				  }
			    }
			 std::cout<< s.is_infeasible()<<std::endl;
  return single_mold_translational_casting_3_impl(polyhedron, np, oi,
                                                  Is_direction());
}

} // end of namespace Set_movable_separability_3
} // end of namespace CGAL

#endif
