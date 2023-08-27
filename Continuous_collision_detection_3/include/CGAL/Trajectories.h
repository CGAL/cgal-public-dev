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

#ifndef TRAJECTORIES_H
#define TRAJECTORIES_H

#include <tuple>

namespace CGAL {



template <class K>
class Point_3_trajectory : public std::tuple<typename K::Point_3, typename K::Point_3> {

    typedef typename K::Point_3               Point;
    typedef          std::tuple<Point, Point> Base;

    public:
        Point_3_trajectory() {}

        Point_3_trajectory(const Point& current, const Point& next)
        : Base(current, next)
        {}

        Point_3_trajectory(const Base& point_pair)
        : Base(point_pair)
        {}

        friend std::ostream& operator<<(std::ostream& os, Point_3_trajectory const& point_trajectory)
        {
            return (os << std::get<0>(point_trajectory) << " -> " << std::get<1>(point_trajectory));
        }

        Point current() {
            return std::get<0>(*this);
        }

        Point next() {
            return std::get<1>(*this);
        }
};

template <class K>
class Segment_3_trajectory : public std::tuple<Point_3_trajectory<K>, Point_3_trajectory<K>> {

    using Segment            = typename K::Segment_3;
    using Point_trajectory = Point_3_trajectory<K>;
    using Base               = std::tuple<Point_trajectory, Point_trajectory>;

    public:
        Segment_3_trajectory() {}

        Segment_3_trajectory(const Point_trajectory& source, const Point_trajectory& target)
        : Base(source, target)
        {}

        Segment_3_trajectory(Point_trajectory&& source, Point_trajectory&& target)
        : Base(std::move(source), std::move(target))
        {}

        Segment_3_trajectory(const Base& point_trajectory_pair)
        : Base(point_trajectory_pair)
        {}

        friend std::ostream& operator<<(std::ostream& os, Segment_3_trajectory const& segment_trajectory)
        {
            return (os << std::get<0>(segment_trajectory) << "\n" << std::get<1>(segment_trajectory));
        }

        Segment current() {
            return Segment(
                std::get<0>(*this).current(),
                std::get<1>(*this).current()
            );
        }

        Segment next() {
            return Segment(
                std::get<0>(*this).next(),
                std::get<1>(*this).next()
            );
        }
};

template <class K>
class Triangle_3_trajectory : public std::tuple<Point_3_trajectory<K>, Point_3_trajectory<K>, Point_3_trajectory<K>> {

    using Triangle           = typename K::Triangle_3;
    using Point_trajectory   = Point_3_trajectory<K>;
    using Base               = std::tuple<Point_trajectory, Point_trajectory, Point_trajectory>;

    public:
        Triangle_3_trajectory() {}

        Triangle_3_trajectory(const Point_trajectory& v0, const Point_trajectory& v1, const Point_trajectory& v2)
        : Base(v0, v1, v2)
        {}

        Triangle_3_trajectory(Point_trajectory&& v0, Point_trajectory&& v1, Point_trajectory&& v2)
        : Base(std::move(v0), std::move(v1), std::move(v2))
        {}

        Triangle_3_trajectory(const Base& point_trajectory_triplet)
        : Base(point_trajectory_triplet)
        {}

        friend std::ostream& operator<<(std::ostream& os, Triangle_3_trajectory const& triangle_trajectory)
        {
            return (os << std::get<0>(triangle_trajectory) << "\n" << std::get<1>(triangle_trajectory) << "\n" << std::get<2>(triangle_trajectory));
        }

        Triangle current() {
            return Triangle(
                std::get<0>(*this).current(),
                std::get<1>(*this).current(),
                std::get<2>(*this).current()
            );
        }

        Triangle next() {
            return Triangle(
                std::get<0>(*this).next(),
                std::get<1>(*this).next(),
                std::get<2>(*this).next()
            );
        }
};



}

#endif