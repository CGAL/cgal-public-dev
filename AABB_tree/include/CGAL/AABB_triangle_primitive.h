// Copyright (c) 2012 INRIA Sophia-Antipolis (France).
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
//
//
// Author(s)     : Sebastien Loriot
//


#ifndef CGAL_AABB_TRIANGLE_PRIMITIVE_H_
#define CGAL_AABB_TRIANGLE_PRIMITIVE_H_

#include <CGAL/AABB_primitive.h>
#include <CGAL/result_of.h>
#include <iterator>

namespace CGAL {

namespace internal {

  template <class GeomTraits, class Iterator, int Dimension>
  struct Point_from_triangle_d_iterator_property_map
  {};

  template <class GeomTraits, class Iterator>
    struct Point_from_triangle_d_iterator_property_map<GeomTraits, Iterator, 2>{
      //classical typedefs
      typedef Iterator key_type;
      typedef typename GeomTraits::Point_2 value_type;
      typedef typename cpp11::result_of<
        typename GeomTraits::Construct_vertex_2(typename GeomTraits::Triangle_3,int)
      >::type reference;
      typedef boost::readable_property_map_tag category;

      inline friend
      typename Point_from_triangle_d_iterator_property_map<GeomTraits,Iterator,2>::reference
      get(Point_from_triangle_d_iterator_property_map<GeomTraits,Iterator,2>, Iterator it)
      {
        return typename GeomTraits::Construct_vertex_2()( *it, 0 );
      }
    };

  template <class GeomTraits, class Iterator>
  struct Point_from_triangle_d_iterator_property_map<GeomTraits, Iterator, 3>{
    //classical typedefs
    typedef Iterator key_type;
    typedef typename GeomTraits::Point_3 value_type;
    typedef typename cpp11::result_of<
      typename GeomTraits::Construct_vertex_3(typename GeomTraits::Triangle_3,int)
    >::type reference;
    typedef boost::readable_property_map_tag category;

    inline friend
    typename Point_from_triangle_d_iterator_property_map<GeomTraits,Iterator,3>::reference
    get(Point_from_triangle_d_iterator_property_map<GeomTraits,Iterator,3>, Iterator it)
    {
      return typename GeomTraits::Construct_vertex_3()( *it, 0 );
    }
  };
}//namespace internal


/*!
 * \ingroup PkgAABB_tree
 * Primitive type that uses as identifier an iterator with a triangle as `value_type`.
 * The iterator from which the primitive is built should not be invalided
 * while the AABB tree holding the primitive is in use.
 *
 * \cgalModels `AABBPrimitive`
 *
 * \tparam GeomTraits is a traits class providing the the necessary geometric data types (Ex. Point_2, Point_3, Triangle_2, Triangle_3).
 *         It also provides the functors (Ex. For 2D Construct_vertex_2) that has an operator taking a triangle type. (Ex. Triangle_2 for 2D)
 *         and an integer as parameters and returning a triangle point as a type convertible to `Point_d`.
 * \tparam Iterator is a model of `ForwardIterator` with its value type convertible to Ex. GeomTraits::Triangle_2 for 2D and accordingly.
 * \tparam cache_datum is either `CGAL::Tag_true` or `CGAL::Tag_false`. In the former case,
 *           the datum is stored in the primitive, while in the latter it is
 *           constructed on the fly to reduce the memory footprint.
 *           The default is `CGAL::Tag_false` (datum is not stored).
 *
 * \sa `AABBPrimitive`
 * \sa `AABB_primitive<Id,ObjectPropertyMap,PointPropertyMapPolyhedron,ExternalPropertyMaps,cache_datum>`
 * \sa `AABB_segment_primitive<Iterator,cache_datum>`
 * \sa `AABB_HalfedgeGraph_segment_primitive<HalfedgeGraph,OneHalfedgeGraphPerTree,cache_datum>`
 * \sa `AABB_FaceGraph_triangle_primitive<FaceGraph,OneFaceGraphPerTree,cache_datum>`
 */
template < class GeomTraits,
           class Iterator,
           class cache_datum=Tag_false>
class AABB_triangle_primitive
#ifndef DOXYGEN_RUNNING
  : public AABB_primitive<  Iterator,
                            Input_iterator_property_map<Iterator>,
                            internal::Point_from_triangle_d_iterator_property_map<GeomTraits, Iterator,Iterator::value_type::Ambient_dimension::value>,
                            Tag_false,
                            cache_datum >
#endif
{
  typedef AABB_primitive< Iterator,
                          Input_iterator_property_map<Iterator>,
                          internal::Point_from_triangle_d_iterator_property_map<GeomTraits, Iterator,Iterator::value_type::Ambient_dimension::value>,
                          Tag_false,
                          cache_datum > Base;
public:
  ///Constructor from an iterator
  AABB_triangle_primitive(Iterator it) : Base(it){}
};

}  // end namespace CGAL


#endif // CGAL_AABB_TRIANGLE_PRIMITIVE_H_

