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

#ifndef CGAL_KDOP_TREE_KDOP_TREE_H_
#define CGAL_KDOP_TREE_KDOP_TREE_H_

//#include <CGAL/license/KDOP_tree.h>

#include <CGAL/disable_warnings.h>

#include <vector>
#include <iterator>

#include <CGAL/internal/AABB_tree/Has_nested_type_Shared_data.h>

#include <CGAL/KDOP_tree/internal/KDOP_traversal_traits.h>
#include <CGAL/KDOP_tree/internal/KDOP_node.h>
//#include <CGAL/KDOP_tree/internal/KDOP_search_tree.h>
#include <CGAL/KDOP_tree/internal/Primitive_helper.h>

#include <boost/optional.hpp>
#include <boost/lambda/lambda.hpp>

#ifdef CGAL_HAS_THREADS
#include <CGAL/mutex.h>
#endif

/// \file KDOP_tree.h

//TODO std::forward cannot be resolved?

namespace CGAL {
namespace KDOP_tree {

/// \addtogroup PkgKDOPTree
/// @{

  /**
   * Class KDOP_tree is a static data structure for efficient
   * intersection in 3D. It builds implicitly a
   * hierarchy of k discrete oriented polytopes (an KDOP tree) from a set
   * of 3D geometric objects, and can receive intersection
   * queries, provided that the corresponding predicates are
   * implemented in the traits class KDOPTraits.
   * An instance of the class `KDOPTraits` is internally stored.
   *
   * \sa `KDOPTraits`
   * \sa `KDOPPrimitive`
   *
   */
  template <typename KDOPTraits>
  class KDOP_tree
  {
  private:
    // internal KD-tree used to accelerate the distance queries
    //typedef internal::KDOP_search_tree<KDOPTraits> Search_tree;

    // type of the primitives container
    typedef std::vector<typename KDOPTraits::Primitive> Primitives;

  public:
    typedef KDOPTraits KDOP_traits;

    /// \name Types
    ///@{

    /// Number type of the geometry kernel.
    typedef typename KDOPTraits::FT FT;


    /// Type of 3D point.
    typedef typename KDOPTraits::Point_3 Point;

    /// Type of input primitive.
    typedef typename KDOPTraits::Primitive Primitive;
    /// Identifier for a primitive in the tree.
    typedef typename Primitive::Id Primitive_id;
    /// Unsigned integer size type.
    typedef typename Primitives::size_type size_type;
    /// Type of kdop.
    typedef typename KDOPTraits::Kdop Kdop;
    /// 3D Point and Primitive Id type
    typedef typename KDOPTraits::Point_and_primitive_id Point_and_primitive_id;

    const static unsigned int num_directions = Kdop::num_directions;

    /*!
    An alias to `KDOPTraits::Intersection_and_primitive_id<Query>`
    @tparam Query should be the type of primitives.
    */
    #ifdef DOXYGEN_RUNNING
    template<typename Query>
    using Intersection_and_primitive_id = KDOPTraits::Intersection_and_primitive_id<Query>;
    #else
    template<typename Query>
    struct Intersection_and_primitive_id {
      typedef typename KDOPTraits::template Intersection_and_primitive_id<Query>::Type Type;
    };
    #endif


    ///@}

  public:
    /// \name Creation
    ///@{

    /// Construct an empty tree, and initializes the internally stored traits
    /// class using `traits`.
    KDOP_tree(const KDOPTraits& traits = KDOPTraits());

    /**
     * @brief Build the datastructure from a sequence of primitives.
     * @param first iterator over first primitive to insert
     * @param beyond past-the-end iterator
     *
     * It is equivalent to constructing an empty tree and calling `insert(first,last,t...)`.
     * The tree stays empty if the memory allocation is not successful.
     */
    template<typename InputIterator,typename ... T>
    KDOP_tree(InputIterator first, InputIterator beyond,T&& ...);

    /// After one or more calls to `insert()` the internal data
    /// structure of the tree must be reconstructed. This procedure
    /// has a complexity of \f$O(n log(n))\f$, where \f$n\f$ is the number of
    /// primitives of the tree.  This procedure is called implicitly
    /// at the first call to a query member function. You can call
    /// `build()` explicitly to ensure that the next call to
    /// query functions will not trigger the reconstruction of the
    /// data structure.
    /// A call to `KDOPTraits::set_shared_data(t...)`
    /// is made using the internally stored traits.
    template<typename ... T>
    void build(T&& ...);
#ifndef DOXYGEN_RUNNING
    void build();
#endif
    ///@}

    /// \name Operations
    ///@{

    /// Equivalent to calling `clear()` and then `insert(first,last,t...)`.
    template<typename ConstPrimitiveIterator,typename ... T>
    void rebuild(ConstPrimitiveIterator first, ConstPrimitiveIterator beyond,T&& ...);


    /// Add a sequence of primitives to the set of primitives of the KDOP tree.
    /// `%InputIterator` is any iterator and the parameter pack `T` are any types
    /// such that `Primitive` has a constructor with the following signature:
    /// `Primitive(%InputIterator, T...)`.
    template<typename InputIterator,typename ... T>
    void insert(InputIterator first, InputIterator beyond,T&& ...);

    /// Adds a primitive to the set of primitives of the tree.
    inline void insert(const Primitive& p);

    /// Clears and destroys the tree.
    ~KDOP_tree()
    {
      clear();
    }
    /// Returns a const reference to the internally stored traits class.
    const KDOPTraits& traits() const{
      return m_traits;
    }

    /// Clears the tree.
    void clear()
    {
      // clear KDOP tree
      clear_nodes();
      m_primitives.clear();
      //clear_search_tree();
      //m_default_search_tree_constructed = false;
    }

    // set parameters for k-dop tree
    void set_kdop_directions(std::vector< Point > directions) {
      m_directions = directions;
      m_num_directions = directions.size();
    }

    /// Returns the kdop of the whole tree.
    /// \pre '!empty()'
    const Kdop kdop() const {
      CGAL_precondition(!empty());
      if (size() > 1) {
        return root_node()->kdop();
      }
      else {
        return KDOP_traits().compute_kdop_object()(m_primitives[0], m_directions);
      }
    }

    /// Output the kdops
    void kdop_heights(std::vector< typename Kdop::Array_height >& heights) {
      CGAL_precondition(!empty());
      if (size() > 1) {
        typedef Compute_kdop_traits<KDOP_traits> Traversal_traits;
        Traversal_traits traversal_traits(m_traits);
        root_node()->template kdop_heights<Traversal_traits>(traversal_traits, m_primitives.size(), m_directions, heights);
      }
      else {
        Kdop kdop_primitive = KDOP_traits().compute_kdop_object()(m_primitives[0], m_directions);
        typename Kdop::Array_height heights_kdop = kdop_primitive.support_heights();
        heights.push_back(heights_kdop);
      }
    }

    /// Returns the number of primitives in the tree.
    size_type size() const { return m_primitives.size(); }

    /// Returns \c true, iff the tree contains no primitive.
    bool empty() const { return m_primitives.empty(); }
    ///@}

  private:
    template <typename ... T>
    void set_primitive_data_impl(CGAL::Boolean_tag<false>,T ... ){}
    template <typename ... T>
    void set_primitive_data_impl(CGAL::Boolean_tag<true>,T&& ... t)
    {m_traits.set_shared_data(std::forward<T>(t)...);}

    template <typename ... T>
    void set_shared_data(T&& ...t){
      set_primitive_data_impl(CGAL::Boolean_tag<CGAL::internal::Has_nested_type_Shared_data<Primitive>::value>(),std::forward<T>(t)...);
    }

    bool build_kd_tree() const;
    template<typename ConstPointIterator>
    bool build_kd_tree(ConstPointIterator first, ConstPointIterator beyond) const;
public:

    /// \name Intersection Tests
    ///@{

    /// Returns `true`, iff the query intersects at least one of
    /// the input primitives. \tparam Query must be a type for
    /// which `do_intersect` predicates are
    /// defined in the traits class `KDOPTraits`.
    template<typename Query>
    bool do_intersect(const Query& query) const;

    /// Returns the number of primitives intersected by the
    /// query. \tparam Query must be a type for which
    /// `do_intersect` predicates are defined
    /// in the traits class `KDOPTraits`.
    template<typename Query>
    size_type number_of_intersected_primitives(const Query& query) const;

    /// Outputs to the iterator the list of all intersected primitives
    /// ids. This function does not compute the intersection points
    /// and is hence faster than the function `all_intersections()`
    /// function below. \tparam Query must be a type for which
    /// `do_intersect` predicates are defined
    /// in the traits class `KDOPTraits`.
    template<typename Query, typename OutputIterator>
    OutputIterator all_intersected_primitives(const Query& query, OutputIterator out) const;


    /// Returns the intersected primitive id that is encountered first
    /// in the tree traversal, iff
    /// the query intersects at least one of the input primitives. No
    /// particular order is guaranteed over the tree traversal, such
    /// that, e.g, the primitive returned is not necessarily the
    /// closest from the source point of a ray query. \tparam Query
    /// must be a type for which
    /// `do_intersect` predicates are defined
    /// in the traits class `KDOPTraits`.
    template <typename Query>
    boost::optional<Primitive_id> any_intersected_primitive(const Query& query) const;
    ///@}

    /// \name Intersections
    ///@{

    /// Outputs the list of all intersections, as objects of
    /// `Intersection_and_primitive_id<Query>::%Type`,
    /// between the query and the input data to
    /// the iterator. `do_intersect()`
    /// predicates and intersections must be defined for `Query`
    /// in the `KDOPTraits` class.
    template<typename Query, typename OutputIterator>
    OutputIterator all_intersections(const Query& query, OutputIterator out) const;


    /// Returns the intersection that is encountered first
    /// in the tree traversal. No particular
    /// order is guaranteed over the tree traversal, e.g, the
    /// primitive returned is not necessarily the closest from the
    /// source point of a ray query. Type `Query` must be a type
    /// for which `do_intersect` predicates
    /// and intersections are defined in the traits class KDOPTraits.
    template <typename Query>
    boost::optional< typename Intersection_and_primitive_id<Query>::Type >
    any_intersection(const Query& query) const;



    /// Returns the intersection and  primitive id closest to the source point of the ray
    /// query.
    /// \tparam Ray must be the same as `KDOPTraits::Ray_3` and
    /// `do_intersect` predicates and intersections for it must be
    /// defined.
    /// \tparam Skip a functor with an operator
    /// `bool operator()(const Primitive_id& id) const`
    /// that returns `true` in order to skip the primitive.
    /// Defaults to a functor that always returns `false`.
    ///
    /// \note `skip` might be given some primitives that are not intersected by `query`
    ///       because the intersection test is done after the skip test. Also note that
    ///       the order the primitives are given to `skip` is not necessarily the
    ///       intersection order with `query`.
    ///
    ///
    /// `KDOPTraits` must be a model of `KDOPRayIntersectionTraits` to
    /// call this member function.
    template<typename Ray, typename SkipFunctor>
    boost::optional< typename Intersection_and_primitive_id<Ray>::Type >
    first_intersection(const Ray& query, const SkipFunctor& skip) const;

    /// \cond
    template<typename Ray>
    boost::optional< typename Intersection_and_primitive_id<Ray>::Type >
    first_intersection(const Ray& query) const
    {
      return first_intersection(query, boost::lambda::constant(false));
    }
    /// \endcond

    /// Returns the primitive id closest to the source point of the ray
    /// query.
    /// \tparam Ray must be the same as `KDOPTraits::Ray_3` and
    /// `do_intersect` predicates and intersections for it must be
    /// defined.
    /// \tparam Skip a functor with an operator
    /// `bool operator()(const Primitive_id& id) const`
    /// that returns `true` in order to skip the primitive.
    /// Defaults to a functor that always returns `false`.
    ///
    /// `KDOPTraits` must be a model of `KDOPRayIntersectionTraits` to
    /// call this member function.
    template<typename Ray, typename SkipFunctor>
    boost::optional<Primitive_id>
    first_intersected_primitive(const Ray& query, const SkipFunctor& skip) const;

    /// \cond
    template<typename Ray>
    boost::optional<Primitive_id>
    first_intersected_primitive(const Ray& query) const
    {
      return first_intersected_primitive(query, boost::lambda::constant(false));
    }
    /// \endcond
    ///@}

  private:
    template<typename KDOPTree, typename SkipFunctor>
    friend class KDOP_ray_intersection;

    // clear nodes
    void clear_nodes()
    {
      if( size() > 1 ) {
        delete [] m_p_root_node;
      }
      m_p_root_node = NULL;
    }

    /*
    // clears internal KD tree
    void clear_search_tree() const
    {
      if ( m_search_tree_constructed )
      {
        CGAL_assertion( m_p_search_tree!=NULL );
        delete m_p_search_tree;
        m_p_search_tree = NULL;
        m_search_tree_constructed = false;
                        }
    }
    */

  public:

    /// \internal
    template <class Traversal_traits>
    void kdop_traversal(Traversal_traits& traits)
    {
      switch(size())
      {
      case 0:
        break;
      case 1:
        traits.compute_kdop(m_primitives[0], m_directions);
        break;
      default:
        root_node()->template kdop_traversal<Traversal_traits>(traits, m_primitives.size(), m_directions);
      }
    }

    template <class QueryPair, class Traversal_traits>
    void traversal(const QueryPair& query_pair, Traversal_traits& traits) const
    {
      switch(size())
      {
      case 0:
        break;
      case 1:
        traits.intersection(query_pair.first, m_primitives[0]);
        break;
      default: // if(size() >= 2)
        root_node()->template traversal<Traversal_traits,QueryPair>(query_pair, traits, m_primitives.size());
      }
    }

  private:
    // parameters for k-dop computations
    int m_num_directions;
    std::vector< Point > m_directions;

    typedef internal::KDOP_node<KDOPTraits> Node;

  public:
    // returns a point which must be on one primitive
    Point_and_primitive_id any_reference_point_and_id() const
    {
      CGAL_assertion(!empty());
      return Point_and_primitive_id(
       internal::Primitive_helper<KDOP_traits>::get_reference_point(m_primitives[0],m_traits), m_primitives[0].id()
      );
    }

    /*
  public:
    Point_and_primitive_id best_hint(const Point& query) const
    {
      if(m_search_tree_constructed)
                        {
        return m_p_search_tree->closest_point(query);
                        }
      else
        return this->any_reference_point_and_id();
    }
    */

    //! Returns the datum (geometric object) represented `p`.
#ifndef DOXYGEN_RUNNING
    typename internal::Primitive_helper<KDOPTraits>::Datum_type
#else
    typename KDOPTraits::Primitive::Datum_reference
#endif
    datum(Primitive& p)const
    {
      return internal::Primitive_helper<KDOPTraits>::get_datum(p, this->traits());
    }

  private:
    //Traits class
    KDOPTraits m_traits;
    // set of input primitives
    Primitives m_primitives;
    // single root node
    Node* m_p_root_node;
    #ifdef CGAL_HAS_THREADS
    mutable CGAL_MUTEX internal_tree_mutex;//mutex used to protect const calls inducing build()
    mutable CGAL_MUTEX kd_tree_mutex;//mutex used to protect calls to accelerate_distance_queries
    #endif

    Node* root_node() const {
      CGAL_assertion(size() > 1);
      if(m_need_build){
        #ifdef CGAL_HAS_THREADS
        //this ensures that build() will be called once
        CGAL_SCOPED_LOCK(internal_tree_mutex);
        if(m_need_build)
        #endif
          const_cast< KDOP_tree<KDOPTraits>* >(this)->build();
      }
      return m_p_root_node;
    }

    const Primitive& singleton_data() const {
      CGAL_assertion(size() == 1);
      return *m_primitives.begin();
    }

    // search KD-tree
    //mutable const Search_tree* m_p_search_tree;
    //mutable bool m_search_tree_constructed;
    //mutable bool m_default_search_tree_constructed; // indicates whether the internal kd-tree should be built
    bool m_need_build;

  private:
    // Disabled copy constructor & assignment operator
    typedef KDOP_tree<KDOPTraits> Self;
    KDOP_tree(const Self& src);
    Self& operator=(const Self& src);

  };  // end class KDOP_tree

/// @}

  template<typename Tr>
  KDOP_tree<Tr>::KDOP_tree(const Tr& traits)
    : m_traits(traits)
    , m_primitives()
    , m_p_root_node(NULL)
    //, m_p_search_tree(NULL)
    //, m_search_tree_constructed(false)
    //, m_default_search_tree_constructed(false)
    , m_need_build(false)
    , m_num_directions(6) // default number of directions = 6
    , m_directions()
  {}

  template<typename Tr>
  template<typename ConstPrimitiveIterator, typename ... T>
  KDOP_tree<Tr>::KDOP_tree(ConstPrimitiveIterator first,
                           ConstPrimitiveIterator beyond,
                           T&& ... t)
    : m_traits()
    , m_primitives()
    , m_p_root_node(NULL)
    //, m_p_search_tree(NULL)
    //, m_search_tree_constructed(false)
    //, m_default_search_tree_constructed(false)
    , m_need_build(false)
    , m_num_directions(6) // default number of directions = 6
    , m_directions()
  {
    // Insert each primitive into tree
    insert(first, beyond,std::forward<T>(t)...);
  }

  template<typename Tr>
  template<typename ConstPrimitiveIterator, typename ... T>
  void KDOP_tree<Tr>::insert(ConstPrimitiveIterator first,
                             ConstPrimitiveIterator beyond,
                             T&& ... t)
  {
    set_shared_data(std::forward<T>(t)...);
    while(first != beyond)
    {
      m_primitives.push_back(Primitive(first,std::forward<T>(t)...));
      ++first;
    }
    m_need_build = true;
  }

  // Clears tree and insert a set of primitives
  template<typename Tr>
  template<typename ConstPrimitiveIterator, typename ... T>
  void KDOP_tree<Tr>::rebuild(ConstPrimitiveIterator first,
                              ConstPrimitiveIterator beyond,
                              T&& ... t)
  {
    // cleanup current tree and internal KD tree
    clear();

    // inserts primitives
    insert(first, beyond,std::forward<T>(t)...);

    build();
  }

        template<typename Tr>
        template<typename ... T>
        void KDOP_tree<Tr>::build(T&& ... t)
        {
          set_shared_data(std::forward<T>(t)...);
          build();
        }

  template<typename Tr>
  void KDOP_tree<Tr>::insert(const Primitive& p)
  {
    m_primitives.push_back(p);
    m_need_build = true;
  }

  // Build the data structure, after calls to insert(..)
  template<typename Tr>
  void KDOP_tree<Tr>::build()
  {
    clear_nodes();

    if(m_primitives.size() > 1) {

      // allocates tree nodes
      m_p_root_node = new Node[m_primitives.size()-1]();
      if(m_p_root_node == NULL)
      {
        std::cerr << "Unable to allocate memory for KDOP tree" << std::endl;
        CGAL_assertion(m_p_root_node != NULL);
        m_primitives.clear();
        clear();
      }

      // constructs the tree
      m_p_root_node->expand(m_primitives.begin(), m_primitives.end(),
                m_primitives.size(), m_traits);

      m_need_build = false;

      // compute k-dops of nodes after splitting
      Compute_kdop_traits<KDOP_traits> traversal_traits(m_traits);
      this->kdop_traversal(traversal_traits);

    }

    /*
    // In case the users has switched on the accelerated distance query
    // data structure with the default arguments, then it has to be
    // /built/rebuilt.
    //if(m_default_search_tree_constructed)
      //build_kd_tree();
    m_need_build = false;
    */
  }

  /*
  // constructs the search KD tree from given points
  // to accelerate the distance queries
  template<typename Tr>
  bool KDOP_tree<Tr>::build_kd_tree() const
  {
    // iterate over primitives to get reference points on them
    std::vector<Point_and_primitive_id> points;
    points.reserve(m_primitives.size());
    typename Primitives::const_iterator it;
    for(it = m_primitives.begin(); it != m_primitives.end(); ++it)
      points.push_back( Point_and_primitive_id(
        internal::Primitive_helper<KDOP_traits>::get_reference_point(
            *it,m_traits), it->id() ) );

    // clears current KD tree
    clear_search_tree();
    bool res = build_kd_tree(points.begin(), points.end());
    m_default_search_tree_constructed = true;
    return res;
  }

  // constructs the search KD tree from given points
  // to accelerate the distance queries
  template<typename Tr>
  template<typename ConstPointIterator>
  bool KDOP_tree<Tr>::build_kd_tree(ConstPointIterator first,
    ConstPointIterator beyond) const
  {
    m_p_search_tree = new Search_tree(first, beyond);
                m_default_search_tree_constructed = true;
    if(m_p_search_tree != NULL)
    {
      m_search_tree_constructed = true;
      return true;
    }
    else
    {
      std::cerr << "Unable to allocate memory for accelerating distance queries" << std::endl;
      return false;
    }
  }
  */

  template<typename Tr>
  template<typename Query>
  bool
    KDOP_tree<Tr>::do_intersect(const Query& query) const
  {
    // compute support heights of the query
    Kdop kdop_query;
    kdop_query.compute_support_heights_ray(m_directions, query);

    typedef typename KDOP_tree<Tr>::KDOP_traits KDOPTraits;
    typedef typename std::pair<Query, Kdop> QueryPair;

    QueryPair query_pair = std::make_pair(query, kdop_query);

    Do_intersect_traits<KDOPTraits, Query> traversal_traits(m_traits, m_directions);
    this->traversal(query_pair, traversal_traits);
    return traversal_traits.is_intersection_found();
  }
#ifndef DOXYGEN_RUNNING //To avoid doxygen to consider definition and declaration as 2 different functions (size_type causes problems)
  template<typename Tr>
  template<typename Query>
  typename KDOP_tree<Tr>::size_type
    KDOP_tree<Tr>::number_of_intersected_primitives(const Query& query) const
  {
    using CGAL::KDOP_tree::Counting_output_iterator;
    typedef typename KDOP_tree<Tr>::KDOP_traits KDOPTraits;
    typedef Counting_output_iterator<Primitive_id, size_type> Counting_iterator;

    size_type counter = 0;
    Counting_iterator out(&counter);

    Listing_primitive_traits<KDOPTraits,
      Query, Counting_iterator> traversal_traits(out,m_traits);
    this->traversal(query, traversal_traits);
    return counter;
  }
#endif
  template<typename Tr>
  template<typename Query, typename OutputIterator>
  OutputIterator
    KDOP_tree<Tr>::all_intersected_primitives(const Query& query,
    OutputIterator out) const
  {
    typedef typename KDOP_tree<Tr>::KDOP_traits KDOPTraits;
    Listing_primitive_traits<KDOPTraits,
      Query, OutputIterator> traversal_traits(out,m_traits);
    this->traversal(query, traversal_traits);
    return out;
  }

  template<typename Tr>
  template<typename Query, typename OutputIterator>
  OutputIterator
    KDOP_tree<Tr>::all_intersections(const Query& query,
    OutputIterator out) const
  {
    typedef typename KDOP_tree<Tr>::KDOP_traits KDOPTraits;
    Listing_intersection_traits<KDOPTraits,
      Query, OutputIterator> traversal_traits(out,m_traits);
    this->traversal(query, traversal_traits);
    return out;
  }

  template <typename Tr>
  template <typename Query>
  boost::optional< typename KDOP_tree<Tr>::template Intersection_and_primitive_id<Query>::Type >
    KDOP_tree<Tr>::any_intersection(const Query& query) const
  {
    typedef typename KDOP_tree<Tr>::KDOP_traits KDOPTraits;
    First_intersection_traits<KDOPTraits, Query> traversal_traits(m_traits);
    this->traversal(query, traversal_traits);
    return traversal_traits.result();
  }

  template <typename Tr>
  template <typename Query>
  boost::optional<typename KDOP_tree<Tr>::Primitive_id>
    KDOP_tree<Tr>::any_intersected_primitive(const Query& query) const
  {
    typedef typename KDOP_tree<Tr>::KDOP_traits KDOPTraits;
    First_primitive_traits<KDOPTraits, Query> traversal_traits(m_traits);
    this->traversal(query, traversal_traits);
    return traversal_traits.result();
  }

} // end namespace KDOP
} // end namespace CGAL

#include <CGAL/KDOP_tree/internal/KDOP_ray_intersection.h>

#include <CGAL/enable_warnings.h>

#endif // CGAL_KDOP_TREE_KDOP_TREE_H_