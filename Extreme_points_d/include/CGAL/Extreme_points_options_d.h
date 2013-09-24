// Copyright (c) 2010 ETH Zurich (Switzerland)
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
// $Id:  $
// 
//
// Author(s)     : Christian Helbling

#ifndef CGAL_EXTREME_POINTS_OPTIONS_D_H
#define CGAL_EXTREME_POINTS_OPTIONS_D_H

#include <CGAL/QP_models.h>

namespace CGAL {

/// \ingroup PkgExtremePointsDEnum
/// This is an enumeration type used to specify an extreme point algorithm in `Extreme_points_options_d`.
/// \sa `CGAL::Extreme_points_options_d`
enum Extreme_point_algorithm_d {
    /// This is the default value of the algorithm in `Extreme_points_options_d`, and it lets the implementation 
    /// choose the algorithm that it thinks is most appropriate for the current situation.
    EP_CHOOSE_APPROPRIATE,
    /// This is the straightforward algorithm `CGAL::extreme_points_d_simple`.\ If most of the input points are extreme points, 
    /// this algorithm can be faster than `EP_DULA_HELGASON`.
    EP_SIMPLE,
    /// This is the output-sensitive algorithm from Dulá and Helgason implemented in `CGAL::extreme_points_d_dula_helgason`.\ If a 
    /// small fraction of the input points are extreme points, this algorithm performs significantly better than `EP_SIMPLE`.
    EP_DULA_HELGASON
};

/// \ingroup PkgExtremePointsDClasses
/// @{
/*! 
   The class `Extreme_points_options_d` is used for passing options to the class `Extreme_points_d<Traits>`.

   \sa `CGAL::Extreme_points_d<Traits>`
   \sa `CGAL::Extreme_point_algorithm_d`
*/
class Extreme_points_options_d {
private:
    Extreme_point_algorithm_d algo_;
    bool deletion_, anti_cycling_;
    Quadratic_program_options qp_options_;
    
public:
    // default constructor
    /// Constructs an instance of `Extreme_points_options_d`. If `algo` is specified it is set as the chosen extreme point algorithm.
    Extreme_points_options_d(Extreme_point_algorithm_d algo = EP_CHOOSE_APPROPRIATE, bool deletion=false, 
                Quadratic_program_options qp_options = Quadratic_program_options())
      : algo_(algo), deletion_(deletion), qp_options_(qp_options) {}
    
    // set/get algorithm
    // ------------------------
    /// Returns the algorithm used for extreme point computations.
    Extreme_point_algorithm_d get_algorithm() const
    {
        return algo_;
    }
    
    /// Sets the algorithm used for extreme point computations to `algo`. For more information see `Extreme_point_algorithm_d`.
    void set_algorithm (Extreme_point_algorithm_d algo)
    {
        algo_ = algo;
    }
    
    // set/get deletion
    // ------------------------
    /// Returns whether or not deletion is permitted
    bool get_deletion() const
    {
        return deletion_;
    }
    
    /// Sets whether or not deletion is permitted
    void set_deletion (bool deletion)
    {
        deletion_ = deletion;
    }
    
    // set/get anti-cycling
    // ------------------------
    /// Returns whether or not anti-cycling is activated
    bool get_anti_cycling() const
    {
        return anti_cycling_;
    }
    
    /// Sets whether or not anti-cycling is activated
    void set_anti_cycling (bool anti_cycling)
    {
        anti_cycling_ = anti_cycling;
        if (anti_cycling_)
          qp_options_.set_pricing_strategy(CGAL::QP_BLAND);
    }

    #ifndef DOXYGEN_RUNNING
    Quadratic_program_options get_qp_options() const
    {
        return qp_options_;
    }
    #endif

};
/// @}


} //namespace CGAL

#endif // CGAL_EXTREME_POINTS_OPTIONS_D_H
