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

  template <class Scene_face_index>
  class Collision_candidate : public std::pair<Scene_face_index, Scene_face_index> {

    using Base = std::pair<Scene_face_index, Scene_face_index>;

    public:

      Collision_candidate(const Scene_face_index& first, const Scene_face_index& second) : Base(first, second) {}
      Collision_candidate(const Base& index_pair) : Base(index_pair) {}

      friend std::ostream& operator<<(std::ostream& os, Collision_candidate const& collision_candidate)
      {
          return (os << collision_candidate.first << " : " << collision_candidate.second );
      } 
    
  };

}

#endif