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

#ifndef COLLISION_CANDIDATE_3_H
#define COLLISION_CANDIDATE_3_H

namespace CGAL{



/*!
\ingroup PkgCollisions3Classes

The class `Collision_candidate` is a container for a pair of primitives that have been identified for potential collision.
*/
template <class Primitive>
class Collision_candidate : public std::pair<const Primitive*, const Primitive*> {

  using Base = std::pair<const Primitive*, const Primitive*>;

  public:

    /// \name Types
    /// @{

    /// @brief The kernel type
    using K     = typename Primitive::K;
    
    /// @brief Type of index used to identify the primitive in the mesh
    using Index = typename Primitive::Index;

    /// @}

    /// \name Creation
    /// @{

    /// @brief Creates a collision candidate from two pointers to primitive objects.
    Collision_candidate(const Primitive* first, const Primitive* second) : Base(first, second) {}

    /// @brief Creates a collision candidate from a pair of pointers to primitive objects.
    Collision_candidate(const Base& index_pair) : Base(index_pair) {}

    /// @}

    /// @brief Outputs the indices of the contained primitives to standard output.
    friend std::ostream& operator<<(std::ostream& os, Collision_candidate const& collision_candidate)
    {
        return (os << collision_candidate.first->index << " : " << collision_candidate.second->index );
    }

};



}

#endif