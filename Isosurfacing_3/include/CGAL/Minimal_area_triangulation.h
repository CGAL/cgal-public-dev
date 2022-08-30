// Copyright (c) 2022  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Ágoston Sipos

#ifndef CGAL_MINIMALAREATRIANGULATION_H
#define CGAL_MINIMALAREATRIANGULATION_H

#include <vector>

#include <CGAL/Simple_cartesian.h>

namespace CGAL {

namespace internal{

// area of 3D triangle
template <class Point_3>
typename Point_3::FT area (
    const Point_3& a,
    const Point_3& b,
    const Point_3& c)
{
    return sqrt(CGAL::cross_product(b - a, c - a).squared_length()) / 2;
}

}

template <class Point_3>
std::vector<std::vector<size_t>> minimal_area_triangulate (
    const std::vector<Point_3>& vertices,
    const std::vector<std::vector<size_t>>& polygon_indices)
{
    using FT = typename Point_3::FT;
    std::vector<std::vector<size_t>> triangle_indices;
    for (auto polygon : polygon_indices) {
        size_t n = polygon.size();
        std::vector<std::pair<FT, size_t>> values(n*n, std::make_pair(-1, -1));
        std::function<FT(size_t, size_t)> M = [&values, n, &M, &vertices, &polygon] (size_t i, size_t j) {
            if (values[i*n+j].second != -1) return values[i*n+j].first;
            else if (i == j || i == (j+1) % n) {
                values[i*n+j].first = 0.0; return 0.0;
            }
            else {
                FT min = std::numeric_limits<FT>::max();
                size_t ind;
                for (size_t k = j + 1; k < i; ++k) {
                    FT val = M(i,k) + M(k,j) +
                        internal::area(vertices[polygon[i]], vertices[polygon[j]], vertices[polygon[k]]);
                    if (val < min){ min = val; ind = k; }
                }
                values[i*n+j] = std::make_pair(min, ind); return min;
            }
        };
        for (size_t i = 0; i < n; ++i) for (size_t j = 0; j < i; ++j) M(i,j);
        std::function<void(size_t, size_t)> f = [&triangle_indices, &polygon, &values, n, &f] (size_t i, size_t j) {
            if (i == j || i == (j+1)%n) return;
            else {
                size_t k = values[i*n+j].second;
                triangle_indices.push_back(std::vector<size_t>{polygon[i],polygon[j],polygon[k]});
                f(i,k), f(k,j);
                return;
            }
        };
        f(n-1, 0);
    }
    return triangle_indices;
}

}

#endif