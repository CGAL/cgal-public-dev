// Copyright (c) 2023
// INRIA Sophia-Antipolis (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jeffrey Cochran

#ifndef COLLISION_TYPE_H
#define COLLISION_TYPE_H

namespace CGAL {
namespace Collisions {
namespace internal {

enum class COLLISION_TYPE {
    EDGE_EDGE, 
    POINT_TRIANGLE
};

}
}
}

#endif