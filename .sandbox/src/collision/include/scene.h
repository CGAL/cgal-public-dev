// Copyright (C) 2001-2005 by Computer Graphics Group, RWTH Aachen
// Copyright (C) 2011 by Graphics & Geometry Group, Bielefeld University
// Copyright (C) 2014 GeometryFactory
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//

#ifndef CGAL_SCENE_H
#define CGAL_SCENE_H

#include <CGAL/license/Surface_mesh.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Surface_mesh/IO.h>
#include <CGAL/Surface_mesh/Surface_mesh_fwd.h>
#include <CGAL/Surface_mesh/Properties.h>

#include <CGAL/assertions.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/circulator.h>
#include <CGAL/Handle_hash_function.h>
#include <CGAL/IO/Verbose_ostream.h>
#include <CGAL/Iterator_range.h>
#include <CGAL/property_map.h>

#include <boost/cstdint.hpp>
#include <boost/array.hpp>
#include <boost/iterator/iterator_facade.hpp>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <iterator>
#include <iostream>
#include <sstream>
#include <string>
#include <typeinfo>
#include <utility>
#include <vector>

namespace CGAL {

template <typename R>
class Scene
{
    private:
        typedef R::Point_3          Point;
        typedef Surface_mesh<Point> Mesh;
    
    public:
        void add_mesh(Mesh& mesh);
        void update(const UpdateScene& update_scene );

        Mesh[] 



}

template<typename P>
void Scene::add_mesh(Surface_mesh<P>& mesh) {
    return;
}

} // namespace CGAL


#include <CGAL/enable_warnings.h>

#endif /* CGAL_SCENE_H */
