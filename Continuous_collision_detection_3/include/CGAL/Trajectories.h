// Copyright (c) 2023 GeometryFactory (France).
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jeffrey Cochran

#ifndef TRAJECTORIES_H
#define TRAJECTORIES_H

#include <tuple>

namespace CGAL {

/// \ingroup PkgCollisions3Classes
/// @{

/// \brief The class `Point_3_trajectory` serves as a container for the current and next positions of a `Point_3` object
/// during the time of interest for a collision query.
template <class K>
class Point_3_trajectory : public std::tuple<typename K::Point_3, typename K::Point_3> {

    typedef typename K::Point_3               Point;
    typedef          std::tuple<Point, Point> Base;

    public:
        /// \name Creation
        /// @{

        /// @brief Creates an empty trajectory.
        Point_3_trajectory() {}

        /// @brief Creates a point trajectory for which the initial position corresponds to `current` and the final position corresponds to `next`.
        Point_3_trajectory(const Point& current, const Point& next)
        : Base(current, next)
        {}

        /// @brief Creates a point trajectory for which the initial position corresponds to the first element of the pair and the final position corresponds to the second.
        Point_3_trajectory(const Base& point_pair)
        : Base(point_pair)
        {}
        
        /// @}

        /// @brief Writes the current and next positions of the point trajectory to the standard output
        friend std::ostream& operator<<(std::ostream& os, Point_3_trajectory const& point_trajectory)
        {
            return (os << std::get<0>(point_trajectory) << " -> " << std::get<1>(point_trajectory));
        }

        /// \name Methods
        /// @{

        /// @brief Returns a point corresponding to the initial position of the trajectory.
        Point current() const {
            return std::get<0>(*this);
        }

        /// @brief Returns a point corresponding to the final position of the trajectory.
        Point next() const {
            return std::get<1>(*this);
        }

        /// @}
};

/// \brief The class `Segment_3_trajectory` serves as a container for the current and next positions of a `Segment_3` object
/// during the time of interest for a collision query.
template <class K>
class Segment_3_trajectory : public std::tuple<Point_3_trajectory<K>, Point_3_trajectory<K>> {

    using Segment            = typename K::Segment_3;
    using Point_trajectory   = Point_3_trajectory<K>;
    using Base               = std::tuple<Point_trajectory, Point_trajectory>;

    public:
        /// \name Creation
        /// @{

        /// @brief  Creates an empty trajectory.
        Segment_3_trajectory() {}

        /// @brief  Creates a segment trajectory characterized the trajectory of its two endpoints.
        Segment_3_trajectory(const Point_trajectory& source, const Point_trajectory& target)
        : Base(source, target)
        {}

        /// @brief  Creates a segment trajectory characterized the trajectory of its two endpoints.
        Segment_3_trajectory(Point_trajectory&& source, Point_trajectory&& target)
        : Base(std::move(source), std::move(target))
        {}

        /// @brief  Creates a segment trajectory characterized by the pair of trajectories corresponding to its endpoints.
        Segment_3_trajectory(const Base& point_trajectory_pair)
        : Base(point_trajectory_pair)
        {}

        /// @}

        /// @brief Writes the initial and final configurations of the segment trajectory to the standard output
        friend std::ostream& operator<<(std::ostream& os, Segment_3_trajectory const& segment_trajectory)
        {
            return (os << std::get<0>(segment_trajectory) << "\n" << std::get<1>(segment_trajectory));
        }

        /// \name Methods
        /// @{

        /// @brief Returns the segment corresponding to the configuration at the start of trajectory. 
        Segment current() const {
            return Segment(
                std::get<0>(*this).current(),
                std::get<1>(*this).current()
            );
        }

        /// @brief Returns the segment corresponding to the configuration at the end of trajectory. 
        Segment next() const {
            return Segment(
                std::get<0>(*this).next(),
                std::get<1>(*this).next()
            );
        }

        /// @}
};

/// @brief The class `Triangle_3_trajectory` serves as a container for the current and next positions of a `Triangle_3` object
/// during the time of interest for a collision query. 
template <class K>
class Triangle_3_trajectory : public std::tuple<Point_3_trajectory<K>, Point_3_trajectory<K>, Point_3_trajectory<K>> {

    using Triangle           = typename K::Triangle_3;
    using Point_trajectory   = Point_3_trajectory<K>;
    using Base               = std::tuple<Point_trajectory, Point_trajectory, Point_trajectory>;

    public:
        /// \name Creation
        /// @{

        /// @brief Creates an empty trajectory.
        Triangle_3_trajectory() {}

        /// @brief Creates a triangle trajectory characterized by the point trajectories of its three vertices.
        Triangle_3_trajectory(const Point_trajectory& v0, const Point_trajectory& v1, const Point_trajectory& v2)
        : Base(v0, v1, v2)
        {}

        /// @brief Creates a triangle trajectory characterized by the point trajectories of its three vertices.
        Triangle_3_trajectory(Point_trajectory&& v0, Point_trajectory&& v1, Point_trajectory&& v2)
        : Base(std::move(v0), std::move(v1), std::move(v2))
        {}

        /// @brief Creates a triangle trajectory characterized by a tuple of three point trajectories corresponding to its three vertices.
        Triangle_3_trajectory(const Base& point_trajectory_triplet)
        : Base(point_trajectory_triplet)
        {}

        /// @}

        /// @brief Writes the initial and final configurations of the triangle trajectory to the standard output
        friend std::ostream& operator<<(std::ostream& os, Triangle_3_trajectory const& triangle_trajectory)
        {
            return (os << std::get<0>(triangle_trajectory) << "\n" << std::get<1>(triangle_trajectory) << "\n" << std::get<2>(triangle_trajectory));
        }

        /// \name Methods
        /// @{

        /// @brief Returns a triangle corresponding to the configuration at the start of the trajectory. 
        Triangle current() const {
            return Triangle(
                std::get<0>(*this).current(),
                std::get<1>(*this).current(),
                std::get<2>(*this).current()
            );
        }

        /// @brief Returns a triangle corresponding to the configuration at the end of the trajectory.
        Triangle next() const {
            return Triangle(
                std::get<0>(*this).next(),
                std::get<1>(*this).next(),
                std::get<2>(*this).next()
            );
        }

        /// @}
};

/// @}

}

#endif