// Copyright (c) 2014-2015 GeometryFactory (France).
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
// 
// Author(s) : Jane Tournois
//

#ifndef CGAL_LLOYD_OPTIMIZE_MESH_2_H
#define CGAL_LLOYD_OPTIMIZE_MESH_2_H

#include <CGAL/Mesh_2/Mesh_global_optimizer_2.h>
#include <CGAL/Mesh_2/Lloyd_move_2.h>
#include <CGAL/Mesh_2/Mesh_sizing_field.h>
#include <CGAL/Mesh_optimization_return_code.h>
#include <CGAL/iterator.h>

#include <fstream>

#define BOOST_PARAMETER_MAX_ARITY 12
#include <boost/parameter.hpp>
#include <boost/parameter/name.hpp>

namespace CGAL
{
namespace parameters
{
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
#pragma clang diagnostic ignored "-Wunused-parameter"
#endif

BOOST_PARAMETER_NAME( cdt )
BOOST_PARAMETER_NAME( (max_iteration_number, tag) max_iteration_number_ )
BOOST_PARAMETER_NAME( (convergence, tag) convergence_)
BOOST_PARAMETER_NAME( (time_limit, tag) time_limit_ )
BOOST_PARAMETER_NAME( (freeze_bound, tag) freeze_bound_)
BOOST_PARAMETER_NAME( (seeds_begin, tag) seeds_begin_)
BOOST_PARAMETER_NAME( (seeds_end, tag) seeds_end_)
BOOST_PARAMETER_NAME( (mark, tag) mark_)

//CGAL_PRAGMA_DIAG_POP
#ifdef __clang__
#pragma clang diagnostic pop
#endif

}//end namespace parameters
}//end namespace CGAL


namespace CGAL
{
  using namespace parameters;

  BOOST_PARAMETER_FUNCTION(
  (Mesh_optimization_return_code),
  lloyd_optimize_mesh_2,
  tag,
  (required (in_out(cdt),*))
  (optional
    (max_iteration_number_, *, 0 )
    (convergence_, *, 0.001 )
    (time_limit_, *, 0. )
    (freeze_bound_, *, 0.001 )
    (seeds_begin_, *, CGAL::Emptyset_iterator())//see comments below
    (seeds_end_, *, CGAL::Emptyset_iterator())//see comments below
    (mark_, *, false) //if "false", seeds indicate "outside" regions
  )
  )
  {
    return lloyd_optimize_mesh_2_impl(cdt,
                                      max_iteration_number_,
                                      convergence_,
                                      freeze_bound_,
                                      time_limit_,
                                      seeds_begin_,
                                      seeds_end_,
                                      mark_);
  }

  /**
  * this partial specialization is a workaround
  * to avoid compilation errors when seeds_begin and seeds_end are
  * not initialized. Indeed, there is no way to have a
  * "default empty iterator" for these named parameters.
  * Emptyset_iterator implements OutputIterator, 
  * but stands here for "any empty input iterator"
  * (and any other type could).
  */
  template<typename CDT>
  Mesh_optimization_return_code
  lloyd_optimize_mesh_2_impl(CDT& cdt,
                             const int max_iterations,
                             const double convergence_ratio,
                             const double freeze_bound,
                             const double time_limit,
                             CGAL::Emptyset_iterator,
                             CGAL::Emptyset_iterator,
                             const bool mark)
  {
    std::list<typename CDT::Point> seeds;
    return lloyd_optimize_mesh_2_impl(cdt, max_iterations, convergence_ratio,
      freeze_bound, time_limit, seeds.begin(), seeds.end(), mark);
  }

  template<typename CDT, typename InputIterator>
  Mesh_optimization_return_code
  lloyd_optimize_mesh_2_impl(CDT& cdt,
                             const int max_iterations,
                             const double convergence_ratio,
                             const double freeze_bound,
                             const double time_limit,
                             InputIterator seeds_begin,
                             InputIterator seeds_end,
                             const bool mark)
  {
    typedef Mesh_2::Mesh_sizing_field<CDT>           Sizing;
    typedef Mesh_2::Lloyd_move_2<CDT, Sizing>        Mv;
    typedef Mesh_2::Mesh_global_optimizer_2<CDT, Mv> Optimizer;

    Optimizer lloyd(cdt,
                    convergence_ratio,
                    freeze_bound);
    lloyd.set_time_limit(time_limit);
    lloyd.set_seeds(seeds_begin, seeds_end, mark);

#ifdef CGAL_MESH_2_OPTIMIZERS_DEBUG
    std::ofstream os("before_lloyd.angles.txt");
    lloyd.output_angles_histogram(os);
    os.close();
#endif

    // 1000 iteration max to avoid infinite loop
    int nb_iterations = (0 == max_iterations)
      ? 1000
      : max_iterations;

    //run optimization
    Mesh_optimization_return_code rc = lloyd(nb_iterations);

#ifdef CGAL_MESH_2_OPTIMIZERS_DEBUG
    std::ofstream os2("after_lloyd.angles.txt");
    lloyd.output_angles_histogram(os2);
    os2.close();
#endif

    return rc;
  }

} //end namespace CGAL

#endif
