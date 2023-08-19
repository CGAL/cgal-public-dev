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

// #ifndef COLLISION_OBJECT_H
// #define COLLISION_OBJECT_H

// #include <cstddef>

// namespace CGAL {

// struct Collision_object {

//     typedef std::size_t size_type;

//     size_type idx_;

//     Mesh_index() {}
//     explicit Mesh_index(size_type idx) : idx_{ idx } {}
    
//     friend std::ostream& operator<<(std::ostream& os, Mesh_index const& mesh_index)
//     {
//         return (os << "m" << mesh_index.idx_ );
//     } 

//     bool operator==(Mesh_index const& other_mesh_index) const {
//         return this->idx_ == other_mesh_index.idx_;
//     }
    
//     bool operator<(Mesh_index const& other_mesh_index) const {
//         return this->idx_ < other_mesh_index.idx_;
//     }
    
//     operator size_type() const { return idx_; }

//     size_type idx() const {
//         return idx_;
//     }
    
//     size_type id() const {
//         return idx_;
//     }
// };

// }

// #endif