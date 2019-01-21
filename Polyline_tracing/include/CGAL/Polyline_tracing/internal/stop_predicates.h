// Copyright (c) 2017, 2018 GeometryFactory (France).
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
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_POLYLINE_TRACING_STOP_PREDICATE_H
#define CGAL_POLYLINE_TRACING_STOP_PREDICATE_H

#include <CGAL/assertions.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <iostream>

namespace CGAL {

namespace Polyline_tracing {

namespace internal {

template <typename MotorcycleGraph>
struct Default_stop_predicate
{
  typedef bool result_type;

  // By default, it should always stops when encountering a crash
  bool operator()(const int /*mc_id*/, const MotorcycleGraph& /*motorcycle_graph*/) const { return true; }
};

template <typename MotorcycleGraph>
struct Complete_zombifier
{
  typedef bool result_type;

  // This is for complete madness (it will stop at borders still)
  bool operator()(const int /*mc_id*/, const MotorcycleGraph& /*motorcycle_graph*/) const
  {
    return false;
  }
};

} // namespace internal

} // namespace Polyline_tracing

} // namespace CGAL

#endif // CGAL_POLYLINE_TRACING_STOP_PREDICATE_H
