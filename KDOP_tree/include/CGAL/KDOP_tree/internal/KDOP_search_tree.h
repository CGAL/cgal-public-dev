// Copyright (c) 2019  University of Cambridge (UK), GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s) : Xiao Xiao, Fehmi Cirak, Andreas Fabri

#ifndef CGAL_KDOP_TREE_INTERNAL_KDOP_SEARCH_TREE_H_
#define CGAL_KDOP_TREE_INTERNAL_KDOP_SEARCH_TREE_H_

#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>

namespace CGAL {
namespace KDOP_tree {

  template<typename Underlying, typename Id>
  class Add_decorated_point : public Underlying
  {

    //TODO class member functions

  }; // end Add_decorated_point class

  template<typename Traits>
  class KDOP_search_tree
  {
  public:
    typedef typename Traits::FT FT;
    typedef typename Traits::Point_3 Point;
    typedef typename Traits::Primitive Primitive;
    typedef typename Traits::Point_and_primitive_id Point_and_primitive_id;
    typedef typename CGAL::Search_traits_3< Add_decorated_point<Traits, typename Traits::Primitive::Id> > Tree_traits;
    typedef typename CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbour_search;
    typedef typename Neighbour_search::Tree Tree;

  private:
    Tree* m_p_tree;

    Point_and_primitive_id get_p_and_p(const Point_and_primitive_id& p) { return p; }

    Point_and_primitive_id get_p_and_p(const Point& p) { return Point_and_primitive_id(p, typename Primitive::Id()); }

  public:
    // constructor
    template<typename ConstPointIterator>
    KDOP_search_tree(ConstPointIterator begin, ConstPointIterator beyond)
      : m_p_tree(NULL)
    {
      typedef typename Add_decorated_point<Traits, typename Traits::Primitive::Id>::Point_3 Decorated_point;

      std::vector<Decorated_point> points;
      while(begin != beyond) {
        Point_and_primitive_id pp = get_p_and_p(*begin);
        points.push_back(Decorated_point(pp.first, pp.second));
        ++begin;
      }

      m_p_tree = new Tree(points.begin(), points.end());
      if (m_p_tree != NULL) {
        m_p_tree->build();
      }
      else {
        std::cerr << "unable to build the search tree!" << std::endl;
      }

    }

    ~KDOP_search_tree() { delete m_p_tree; }

    Point_and_primitive_id closest_point(const Point& query) const
    {
      Neighbour_search search(*m_p_tree, query, 1);
      return Point_and_primitive_id( static_cast<Point>(search.begin()->first), search.begin()->first.id() );
    }

  };

}
}

#endif // CGAL_KDOP_TREE_INTERNAL_KDOP_SEARCH_TREE_H_
