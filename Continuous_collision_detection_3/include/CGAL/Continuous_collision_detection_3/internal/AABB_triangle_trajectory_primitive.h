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

#ifndef AABB_TRIANGLE_TRAJECTORY_PRIMITVE_H
#define AABB_TRIANGLE_TRAJECTORY_PRIMITVE_H

#include <iostream>
#include <list>
#include <vector>
#include <utility>

#include <CGAL/Trajectories.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>


// custom triangle type with
// three pointers to points


namespace CGAL {
namespace Collisions {
namespace internal {



template<class K_, class Index_>
struct Triangle_trajectory_observer {

    typedef          K_                          K;
    typedef          Index_                      Index;
    typedef typename K_::FT                      FT;
    typedef typename K_::Point_3                 Point;
    typedef typename K_::Iso_cuboid_3            Iso_cuboid_3;
    typedef          Point_3_trajectory<K_>      Point_trajectory;
    typedef          Triangle_3_trajectory<K_>   Triangle_trajectory;

    struct Extrema {
        Point min_point;
        Point max_point;

        Extrema() {}
        Extrema(const Point & p) : min_point{p}, max_point{p} {}

        void update( const Point & p) {
            FT x_min = p.x() < min_point.x() ? p.x() : min_point.x();
            FT x_max = p.x() > max_point.x() ? p.x() : max_point.x();

            FT y_min = p.y() < min_point.y() ? p.y() : min_point.y();
            FT y_max = p.y() > max_point.y() ? p.y() : max_point.y();

            FT z_min = p.z() < min_point.z() ? p.z() : min_point.z();
            FT z_max = p.z() > max_point.z() ? p.z() : max_point.z();

            min_point = Point(x_min, y_min, z_min);
            max_point = Point(x_max, y_max, z_max);
        }
    };

    // Bounding box
    Iso_cuboid_3 bounding_iso_cuboid;
    const Point * pa;
    const Point * pb;
    const Point * pc;
    const Point * next_pa;
    const Point * next_pb;
    const Point * next_pc;
    Index_ index;

    Triangle_trajectory_observer(
        const Point * current_position_a,
        const Point * current_position_b,
        const Point * current_position_c,
        const Point * next_position_a,
        const Point * next_position_b,
        const Point * next_position_c,
        Index_ index
    ) : pa{current_position_a}, pb{current_position_b}, pc{current_position_c}, next_pa{next_position_a}, next_pb{next_position_b}, next_pc{next_position_c}, index{index} {
        bounding_iso_cuboid = compute_bounding_iso_cuboid();
    }

    Iso_cuboid_3 compute_bounding_iso_cuboid() {
        Extrema extrema(*pa);
        extrema.update(*pb);
        extrema.update(*pc);
        extrema.update(*next_pa);
        extrema.update(*next_pb);
        extrema.update(*next_pc);
        return Iso_cuboid_3(
            extrema.min_point, extrema.max_point
        );
    }

    void update(){ bounding_iso_cuboid = compute_bounding_iso_cuboid(); }
};

template <class K, class Index>
Triangle_3_trajectory<K> to_Triangle_3_trajectory(Triangle_trajectory_observer<K, Index> observer) {
    return Triangle_3_trajectory<K>(
        Point_3_trajectory<K>(*observer.pa, *observer.next_pa),
        Point_3_trajectory<K>(*observer.pb, *observer.next_pb),
        Point_3_trajectory<K>(*observer.pc, *observer.next_pc)
    );
}

// The following primitive provides the conversion facilities between
// the custom triangle and point types and the CGAL ones
template <class K, class Index>
struct AABB_Triangle_trajectory_primitive {
public:

    // CGAL types returned
    typedef typename K::FT                                      FT;
    typedef typename K::Point_3                                 Point; // CGAL 3D point type
    typedef typename K::Iso_cuboid_3                            Datum; // CGAL 3D triangle type
    typedef typename K::Aff_transformation_3                    Transform;
    typedef          Triangle_trajectory_observer<K, Index>     Trajectory;
    typedef typename std::vector<Trajectory>::const_iterator    Iterator;

    // this is the type of data that the queries returns. For this example
    // we imagine that, for some reasons, we do not want to store the iterators
    // of the vector, but raw pointers. This is to show that the Id type
    // does not have to be the same as the one of the input parameter of the
    // constructor.
    typedef const Trajectory* Id;

private:
    Id m_pt; // this is what the AABB tree stores internally

public:
    AABB_Triangle_trajectory_primitive() {} // default constructor needed

    // the following constructor is the one that receives the iterators from the
    // iterator range given as input to the AABB_tree
    AABB_Triangle_trajectory_primitive(Iterator it)
        : m_pt(&(*it)) {}

    const Id& id() const { return m_pt; }

    // on the fly conversion from the internal data to the CGAL types
    Datum datum() const
    {
        return m_pt->bounding_iso_cuboid;
    }

    // returns a reference point which must be on the primitive
    Point reference_point() const
    { return *(m_pt->pa); }
};


}
}
} // end CGAL

#endif

