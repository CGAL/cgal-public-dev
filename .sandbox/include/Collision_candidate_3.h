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

#ifndef COLLISION_CANDIDATE_3_H
#define COLLISION_CANDIDATE_3_H

namespace CGAL{

  template <class Primitive>
  class Collision_candidate : public std::pair<const Primitive*, const Primitive*> {

    using Base = std::pair<const Primitive*, const Primitive*>;

    public:

      using K     = typename Primitive::K;
      using Index = typename Primitive::Index;

      Collision_candidate(const Primitive* first, const Primitive* second) : Base(first, second) {}
      Collision_candidate(const Base& index_pair) : Base(index_pair) {}

      friend std::ostream& operator<<(std::ostream& os, Collision_candidate const& collision_candidate)
      {
          return (os << collision_candidate.first->index << " : " << collision_candidate.second->index );
      } 
    
  };

}

#endif