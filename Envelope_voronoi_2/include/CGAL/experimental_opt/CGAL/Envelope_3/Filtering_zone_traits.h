// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: $
// $Id: $
// 
//
// Author(s)     : Ophir Setter        <ophirset@post.cs.tau.ac.il>

/*! \file Filtering_zone_traits.h
*/

#ifndef CGAL_FILTERING_ZONE_TRAITS_H
#define CGAL_FILTERING_ZONE_TRAITS_H

#include <CGAL/Sweep_line_2/Arr_basic_insertion_traits_2.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/Arrangement_2/Arr_traits_adaptor_2.h>

#include <CGAL/Envelope_voronoi_2/envelope_voronoi_assertions.h>


namespace CGAL {

// class template for traits classes that have no Is_between_endpoints object.
// The answer to this question can be calculated with compare_y_at_x.
template <class EnvelopeTraits_>
class Basic_filtered_voronoi_traits 
: public Arr_traits_adaptor_2<EnvelopeTraits_>
{
public:
  typedef EnvelopeTraits_                                EnvelopeTraits;
  typedef typename 
    EnvelopeTraits_::X_monotone_curve_2                  X_monotone_curve_2;
  typedef typename
    EnvelopeTraits_::Point_2                             Point_2;
  typedef Arr_traits_adaptor_2<EnvelopeTraits>           Traits_adaptor_2;

  class Is_between_endpoints
  {
    const Traits_adaptor_2 *_traits_adaptor;
  public:
    Is_between_endpoints(const Traits_adaptor_2 *traits)
      : _traits_adaptor(traits)
    {}

    bool operator() (const X_monotone_curve_2 &c, const Point_2 &p)
    {
      // first check that the point is in the x-range of the curve.
      if (_traits_adaptor->is_in_x_range_2_object() (c, p) == false)
        return false;
      
      return (_traits_adaptor->compare_y_at_x_2_object()(p, c) == ZERO);
    }
  };

  Is_between_endpoints is_between_endpoints_object() const
  {
    return Is_between_endpoints(this);
  }
};

// class template for traits classes that that the answer for 
// Is_between_endpoints is always true, namely have only 1 x-monotone curve
template <class EnvelopeTraits_>
class Linear_filtered_voronoi_traits
: public Arr_traits_adaptor_2<EnvelopeTraits_>
{
public:
  typedef EnvelopeTraits_                                EnvelopeTraits;
  typedef typename 
    EnvelopeTraits_::X_monotone_curve_2                  X_monotone_curve_2;
  typedef typename
    EnvelopeTraits_::Point_2                             Point_2;
  typedef Arr_traits_adaptor_2<EnvelopeTraits>           Traits_adaptor_2;

  class Is_between_endpoints
  {
    const Traits_adaptor_2 *_traits_adaptor;
  public:
    Is_between_endpoints(const Traits_adaptor_2 *traits)
      : _traits_adaptor(traits)
    {}

    bool operator() (const X_monotone_curve_2 &c, const Point_2 &p)
    {
      CGAL_envelope_voronoi_assertion(_traits_adaptor->is_in_x_range_2_object() (c, p) == true);
      CGAL_envelope_voronoi_assertion(_traits_adaptor->compare_y_at_x_2_object()(p, c) == ZERO);
      return true;
    }
  };

  Is_between_endpoints is_between_endpoints_object() const
  {
    return Is_between_endpoints(this);
  }
};

// class that extends an envelope_3 traits class which extends 
// Arr_conic_traits_2 to support the (filtered) Voronoi traits.
template <class EnvelopeTraits_>
class Filtered_voronoi_conic_traits : public EnvelopeTraits_
{
public:
  typedef EnvelopeTraits_                                  Base;
  typedef typename Base::X_monotone_curve_2                X_monotone_curve_2;
  typedef typename Base::Point_2                           Point_2;
  
  class Is_between_endpoints
  {
  public:
    bool operator() (const X_monotone_curve_2 &c, const Point_2 &p)
    {
      return c.is_between_endpoints(p);
    }
  };
  
  Is_between_endpoints is_between_endpoints_object() const
  {
    return Is_between_endpoints();
  }
};

// class that extends an envelope_3 traits class which extends 
// Arr_geodesic_arc_on_sphere_traits_2 to support the (filtered) 
// Voronoi traits.
template <class EnvelopeTraits_>
class Filtered_voronoi_sphere_traits : public EnvelopeTraits_
{
public:
  typedef EnvelopeTraits_                                  Base;
  typedef typename Base::X_monotone_curve_2                X_monotone_curve_2;
  typedef typename Base::Point_2                           Point_2;
  
  class Is_between_endpoints
  {
  public:

    Is_between_endpoints(const Base *traits)
      : m_traits(traits) {}
    
    // this is a temporary solution until we can work without breaking
    // the bisector into 4.
    bool operator() (const X_monotone_curve_2 &xc, const Point_2 &p)
    {
      bool res = xc.is_in_x_range(p);
      if(res == false)
        return false;

      if (xc.is_vertical()) 
      {
        if (!xc.left().is_min_boundary()) 
        {
          Comparison_result cr = m_traits->compare_y(p, xc.left());
          if (cr != LARGER) 
            return false;
        }
        
        if (xc.right().is_max_boundary()) 
          return true;
        Comparison_result cr = m_traits->compare_y(p, xc.right());
        return (cr == LARGER) ? false : true;
      }

      return res;
    }
    
  private:
    const Base      *m_traits;
  };
  
  Is_between_endpoints is_between_endpoints_object() const
  {
    return Is_between_endpoints(this);
  }
};


// This class is used when using the algorithm to filter envelope_3 
// computations and voronoi diagram filtering computations. 
// It can be used when we are dealing with continous and GZIRIM surfaces. 
// The filtering mechanism takes the advantage that if we computed and 
// intersect 2 bisectors of 3 surfaces then the other bisector should 
// pass in the same point. Since in the overlay step of the envelope 
// algorithm created vertices of the arrangement are only intersection 
// between 2 bisectors of 4 sites, we need to use this tratis in merging 
// step of the envelope algorithm.
//
// To be used by this class the Traits_ class should model the FilteredVoronoi
// concept (see above) which allow us to differenciate between points on 
// different x-monotone curves of the same bisector.
//
template <class Traits_, class Arrangement_>
class Filtering_zone_traits : 
  public Arr_basic_insertion_traits_2<Traits_, Arrangement_>
{
public:
  typedef Traits_                                  Traits_2;
  typedef Filtering_zone_traits<Traits_, 
    Arrangement_>                                  Self;
  typedef Arr_basic_insertion_traits_2<Traits_,
    Arrangement_>                                  Base;


  typedef typename Traits_2::Split_2               Base_Split_2;
  typedef typename Traits_2::Intersect_2           Base_Intersect_2;

  typedef typename Base::Point_2                   Base_point_2;
  typedef typename Base::X_monotone_curve_2        Base_x_monotone_curve_2;

  typedef typename Traits_2::Compare_xy_2          Base_compare_xy_2;
  typedef typename Base::Equal_2                   Base_equal_2;
  typedef typename Base::Compare_x_2               Base_compare_x_2;
  typedef typename Base::Compare_x_near_boundary_2 
    Base_compare_x_near_boundary_2;
  typedef typename 
    Base::Construct_min_vertex_2                   Base_construct_min_vertex_2;
  typedef typename 
    Base::Construct_max_vertex_2                   Base_construct_max_vertex_2;
  typedef typename Base::Compare_y_at_x_2          Base_compare_y_at_x_2;
  typedef typename Base::Compare_y_near_boundary_2 
    Compare_y_near_boundary_2;
  typedef typename Traits_2::Compare_y_at_x_left_2 Base_compare_y_at_x_left_2;

  typedef typename Base::Compare_y_on_identification_2
    Base_compare_y_on_identification_2;

  typedef typename Traits_2::Curve_2               Curve_2;
  typedef typename Traits_2::Make_x_monotone_2     Make_x_monotone_2;
  typedef typename Traits_2::Multiplicity          Multiplicity;
  
  typedef typename Base::Boundary_category     Has_boundary_category;
  typedef typename Base::Has_left_category         Has_left_category;
  typedef typename Traits_2::Has_merge_category    Has_merge_category;

  typedef typename Base::Halfedge_handle           Halfedge_handle;
  typedef typename Base::Vertex_handle             Vertex_handle;

  // The traits should be used with minimization diagram, and extended dcel
  // for filtering zone.
  typedef typename Traits_2::Xy_monotone_surface_3 Xy_monotone_surface_3;
  typedef typename Traits_2::Surface_3             Surface_3;
  typedef typename Traits_2::Surface_id            Surface_id;
  typedef typename Traits_2::Is_between_endpoints  Is_between_endpoints;

  // We need to know that a curve is the intersecting curve of the current
  // two surfaces, so we add a flag to the X_monotone class and point class.

  class X_monotone_curve_2 : public Base_x_monotone_curve_2
  {
  public:
    
    X_monotone_curve_2()
    {}
    
  X_monotone_curve_2(const typename Base::Base_x_monotone_curve_2& base, 
                     Halfedge_handle he = Halfedge_handle(),
                     bool current = false)
    : Base_x_monotone_curve_2(base, he), is_current_curve(current)
    {}
    
    // The flag indicate if this is a curve that it is still not 
    // connected to an edge because it is being inserted (was created 
    // by the two current surfaces).
    bool is_current_curve;    
  };

  friend std::ostream& operator<< (std::ostream &out, 
                                   const X_monotone_curve_2& curve)
  {
    return out << curve.base();
  }

  class Point_2 : public Base_point_2
  {
  public:
    
    Point_2()
    {}

    Point_2(const Base_point_2 &p)
      : Base_point_2(p), is_current_point(false)
    {}
    
  Point_2(const typename Base::Base_point_2& base, 
          Vertex_handle vh = Vertex_handle(),
          bool current = false)
    : Base_point_2(base, vh), is_current_point(current)
    {}
    
    // The flag indicate if this is a point that it is still not 
    // connected to a vertex because it is being inserted (was created 
    // by the two current surfaces).
    bool is_current_point;
  };

  friend std::ostream& operator<< (std::ostream &out, 
                                   const Point_2& point)
  {
    return out << point.base();
  }


protected:
  // we keep the current surfaces that we insert their bisector.
  const Xy_monotone_surface_3 &_surf1;
  const Xy_monotone_surface_3 &_surf2;

public:

  
  Filtering_zone_traits (Traits_2& tr, const Xy_monotone_surface_3 &surf1, 
                         const Xy_monotone_surface_3 &surf2)
    : Base (tr), _surf1(surf1), _surf2(surf2)
  { }

    /*! \class
   * The Comapre_y_at_x_left_2 functor.
   */
  class Compare_y_at_x_left_2
  {
  private:
    Base_compare_y_at_x_left_2 m_base_cmp_y_at_x_left;

  public:
    Compare_y_at_x_left_2(const Base_compare_y_at_x_left_2& base):
        m_base_cmp_y_at_x_left(base)
    {}

    Comparison_result operator() (const X_monotone_curve_2& cv1,
                                  const X_monotone_curve_2& cv2,
                                  const Point_2& p) const
    {
      return (m_base_cmp_y_at_x_left(cv1.base(),
                                      cv2.base(),
                                      p.base()));
    }
  };

  Compare_y_at_x_left_2 compare_y_at_x_left_2_object () const
  {
    return (Compare_y_at_x_left_2
	    (this->m_base_traits->compare_y_at_x_left_2_object()));
  }
  
  /*! \class
   * The Compare_y_on_identification_2 functor.
   */
  class Compare_y_on_identification_2
  {
  private:
    Base_compare_y_on_identification_2 m_base_cmp_y_on_idnt;

  public:
    Compare_y_on_identification_2(const Base_compare_y_on_identification_2& base):
        m_base_cmp_y_on_idnt(base)
    {}

    Comparison_result operator() (const Point_2& p1,
                                  const Point_2& p2) const
    {
      return (m_base_cmp_y_on_idnt(p1.base(), p2.base()));
    }
  };

  Compare_y_on_identification_2 compare_y_on_identification_2_object () const
  {
    return (Compare_y_on_identification_2
	    (this->Base::compare_y_on_identification_2_object()));
  }




  /*! \class
   * The Intersect_2 functor.
   */
  class Intersect_2
  {
  private:

    const Self*     _traits;
    Base_Intersect_2 m_base_intersect;
    
  public:
   
    /*! Constructor. */
    Intersect_2 (const Self *traits, Base_Intersect_2 base_inter)
      : m_base_intersect (base_inter)
    {
      _traits = traits;
    }

    template<class OutputIterator>
    OutputIterator operator() (const X_monotone_curve_2& cv1,
                               const X_monotone_curve_2& cv2,
                               OutputIterator oi)
    {
//      std::clog << "intersecting " << cv1 << " and " << cv2 << std::endl;

      if(cv1.halfedge_handle() != Halfedge_handle() &&
         cv2.halfedge_handle() != Halfedge_handle() )
        return (oi); // the curves are disjoint-interior because they
                     // are already at the Arrangement

      if (cv1.halfedge_handle() != Halfedge_handle() ||
          cv2.halfedge_handle() != Halfedge_handle())
      {
        // if one of the curves have an halfedge, we can use the filtering
        // mechanism.
        Halfedge_handle hh;
        const X_monotone_curve_2 *cur;
        if (cv1.halfedge_handle() != Halfedge_handle())
        {
          hh = cv1.halfedge_handle();
          cur = &cv2;
        }
        else
        {
          hh = cv2.halfedge_handle();
          cur = &cv1;
        }
        
        bool s_at_inf = hh->source()->is_at_infinity();
        bool t_at_inf = hh->target()->is_at_infinity();

        Is_between_endpoints is_bet = 
          _traits->m_base_traits->is_between_endpoints_object();
        bool source_on_surfaces = false, target_on_surfaces = false;
        
        if (s_at_inf == false && t_at_inf == false)
        {
          if (hh->source()->is_on_surfaces(
                _traits->_surf1.id, _traits->_surf2.id) &&
              is_bet(*cur, hh->source()->point()))
          {
            source_on_surfaces = true;
          }
          if (hh->target()->is_on_surfaces(
                _traits->_surf1.id, _traits->_surf2.id) &&
              is_bet(*cur, hh->target()->point()))
          {
            target_on_surfaces = true;
          }
        }

        // if both vertices are on the curve, there could be overlaps, and if
        // both vertices are not on the curve there is nothing to do. This is
        // why we want only one vertex to be on the curve.
        if (s_at_inf == false && t_at_inf == false &&
            source_on_surfaces != target_on_surfaces)
        {
          if (source_on_surfaces)
            *oi++ = make_object (std::make_pair<Point_2, Multiplicity>
                                 (Point_2(hh->source()->point(), hh->source()), 1));
          else if (target_on_surfaces)
            *oi++ = make_object (std::make_pair<Point_2, Multiplicity> 
                                 (Point_2(hh->target()->point(), hh->target()), 1));
          
          return oi;
        }
      }

      typedef std::pair<typename Traits_2::Point_2, Multiplicity> Intersection_point;
      typedef typename Traits_2::X_monotone_curve_2    Intersection_curve;
      typedef std::list<CGAL::Object>  Object_list;
      Object_list obj_list;
      m_base_intersect(cv1.base(), cv2.base(), std::back_inserter(obj_list));
      const Intersection_curve                *base_overlap_cv;
      const Intersection_point  *intersect_p;

      // convert objects that are associated with Base_x_monotone_curve_2 to
      // X_monotone_curve_2 
      // OutputIterator concept has no '!=' operator, so we fixed this function.
      Object_list::iterator it;
      for(it = obj_list.begin(); it != obj_list.end(); ++it)
      {
        base_overlap_cv = object_cast<Intersection_curve> (&(*it));
        if (base_overlap_cv != NULL)
        {
          // Add halfedge handles to the resulting curve.
          Halfedge_handle  he;

          if (cv1.halfedge_handle() != Halfedge_handle())
            he = cv1.halfedge_handle();
          else if (cv2.halfedge_handle() != Halfedge_handle())
            he = cv2.halfedge_handle();

          X_monotone_curve_2    overlap_cv (*base_overlap_cv, he);
          overlap_cv.set_overlapping();
          *oi++ = make_object (overlap_cv);
        }
        else
        {
          intersect_p = 
            object_cast<Intersection_point> (&(*it));

          CGAL_envelope_voronoi_assertion (intersect_p != NULL);

          *oi++ = make_object (std::make_pair (Point_2(intersect_p->first),
                                               intersect_p->second));
        }
      }

      // Return a past-the-end iterator.
      return oi;
    }
  };

  Intersect_2 intersect_2_object () const
  {
    return (Intersect_2 (this, this->m_base_traits->intersect_2_object()));
  }


   /*! \class
   * The Split_2 functor.
   */
  class Split_2
  {
  private:
    Base_Split_2    m_base_split;

  public:

    /*! Constructor. */
    Split_2 (const Base_Split_2& base) :
        m_base_split (base)
    {}

    void operator() (const X_monotone_curve_2& cv, const Point_2 & p,
                     X_monotone_curve_2& c1, X_monotone_curve_2& c2)
    {
      m_base_split(cv.base(),
                   p.base(),
                   c1.base(),
                   c2.base());
      c1.is_current_curve = cv.is_current_curve;
      c2.is_current_curve = cv.is_current_curve;
    }
  };

  Split_2 split_2_object () const
  {
    return (Split_2 (this->m_base_traits->split_2_object()));
  }

  

  /*! \class
   * The Equal_2 functor.
   */
  class Equal_2 : public Base_equal_2
  {
  public:

    Equal_2(const Base_equal_2& base) :
        Base_equal_2(base)
    {}

    using Base_equal_2::operator();

    /*! Check if the two points are the same. */
    bool operator() (const Point_2& p1, const Point_2& p2) const
    {
      if(p1.vertex_handle() == p2.vertex_handle() &&
         p1.vertex_handle() != Vertex_handle())
      {
        CGAL_envelope_voronoi_assertion(m_base_eq(p1.base(), p2.base()) == true);
        return true;
      }

      return (m_base_eq(p1.base(), p2.base()));
    }
  };

  Equal_2 equal_2_object () const
  {
    return (Equal_2 (this->Base::equal_2_object()));
  }

  /*! \class
   * The Compare_x_2 functor.
   */
  class Compare_x_2
  {
  private:
    const Self*     _traits;

  public:

  Compare_x_2(const Self* traits)
    : _traits(traits)
    {}

    Comparison_result operator() (const Point_2& p1, 
                                  const Point_2& p2) const
    {
      if(p1.vertex_handle() == p2.vertex_handle() &&
         p1.vertex_handle() != Vertex_handle())
      {
        CGAL_envelope_voronoi_assertion(_traits->Base::compare_x_2_object() (p1, p2) == EQUAL);
        return EQUAL;
      }

      return _traits->Base::compare_x_2_object() (p1, p2);
    }
  };

  Compare_x_2 compare_x_2_object () const
  {
    return (Compare_x_2 (this));
  }

  class Compare_x_near_boundary_2
  {
  private:
    const Self*     _traits;

  public:

  Compare_x_near_boundary_2(const Self* traits)
    : _traits(traits)
    {}
    Comparison_result operator() (const Point_2& p, 
                                  const X_monotone_curve_2& cv, 
                                  Arr_curve_end ce) const
    {
      Vertex_handle vh = p.vertex_handle();
      Halfedge_handle hh = cv.halfedge_handle();
      if (vh != Vertex_handle())
      {
        if (cv.is_current_curve) 
        {
          if (vh->is_on_surfaces(_traits->_surf1.id, _traits->_surf2.id))
          {
            CGAL_envelope_voronoi_assertion(_traits->Base::compare_x_near_boundary_2_object () \
                           (p, cv, ce) == EQUAL);
            return EQUAL;
          }
        }
        else
        {
          CGAL_envelope_voronoi_assertion (hh != Halfedge_handle());
          if (vh->is_on_surfaces(hh->surfaces().begin(), hh->surfaces().end()))
          {
            CGAL_envelope_voronoi_assertion(_traits->Base::compare_x_near_boundary_2_object () \
                           (p, cv, ce) == EQUAL);
            return EQUAL;
          }
        }
      }
      else if (p.is_current_point)
      {
        if (cv.is_current_curve)
        {
          CGAL_envelope_voronoi_assertion(_traits->Base::compare_x_near_boundary_2_object () \
                         (p, cv, ce) == EQUAL);
          return EQUAL;
        }
        else
        {
          CGAL_envelope_voronoi_assertion (hh != Halfedge_handle());
          if ((hh->surfaces().find(_traits->_surf1.id) != hh->surfaces().end()) &&
              (hh->surfaces().find(_traits->_surf2.id) != hh->surfaces().end()))
          {
            CGAL_envelope_voronoi_assertion(_traits->Base::compare_x_near_boundary_2_object () \
                           (p, cv, ce) == EQUAL);
            return EQUAL;
          }
        }
      }
      
      return _traits->Base::compare_x_near_boundary_2_object () (p, cv, ce);
    }
    
    Comparison_result operator() (const X_monotone_curve_2& cv1,
                                  Arr_curve_end ind1,
                                  const X_monotone_curve_2& cv2,
                                  Arr_curve_end ind2) const
    {
      if (cv1.is_current_curve && cv2.is_current_curve)
      {
        CGAL_envelope_voronoi_assertion(_traits->Base::compare_x_near_boundary_2_object () \
                       (cv1, ind1, cv2, ind2) == EQUAL);
        return EQUAL;
      }

      if(cv1.halfedge_handle() == cv2.halfedge_handle() &&
         cv1.halfedge_handle() != Halfedge_handle())
      {
        CGAL_envelope_voronoi_assertion(_traits->Base::compare_x_near_boundary_2_object () \
                       (cv1, ind1, cv2, ind2) == EQUAL);
        return EQUAL;
      }
      
      return _traits->Base::compare_x_near_boundary_2_object () 
        (cv1, ind1, cv2, ind2);
    }
  };


  /*! \class
   * The Construct_min_vertex_2 functor.
   */
  class Construct_min_vertex_2
  {
  private:
    Base_construct_min_vertex_2 m_base_min_v;

  public:

    Construct_min_vertex_2 (const Base_construct_min_vertex_2& base_min_v)
      : m_base_min_v (base_min_v)
    {}
    
    Point_2 operator() (const X_monotone_curve_2 & cv) 
    {
      return Point_2 (m_base_min_v(cv));
    }
  };

  Construct_min_vertex_2 construct_min_vertex_2_object () const
  {
    return Construct_min_vertex_2 (this->Base::construct_min_vertex_2_object());
  }

  /*! \class
   * The Construct_max_vertex_2 functor.
   */
  class Construct_max_vertex_2
  {
  private:
    Base_construct_max_vertex_2 m_base_max_v;

  public:

    Construct_max_vertex_2 (const Base_construct_max_vertex_2& base_max_v)
      : m_base_max_v (base_max_v)
    {}
    
    Point_2 operator() (const X_monotone_curve_2 & cv) 
    {
      return Point_2 (m_base_max_v(cv));
    }
  };

  Construct_max_vertex_2 construct_max_vertex_2_object () const
  {
    return Construct_max_vertex_2 (this->Base::construct_max_vertex_2_object());
  }

  class Compare_y_at_x_2
  {
  private:
    const Self*     _traits;
  public:
  Compare_y_at_x_2(const Self* traits)
    : _traits(traits)
    {}
    
    Comparison_result operator() (const Point_2& p,
                                  const X_monotone_curve_2& cv) const
    {
      Is_between_endpoints is_bet = 
        _traits->m_base_traits->is_between_endpoints_object();

      Vertex_handle vh = p.vertex_handle();
      Halfedge_handle hh = cv.halfedge_handle();
      if (vh != Vertex_handle())
      {
        if (cv.is_current_curve) 
        {
          if (vh->is_on_surfaces(_traits->_surf1.id, _traits->_surf2.id) && 
              is_bet(cv, p))
          {
            CGAL_envelope_voronoi_assertion(_traits->Base::compare_y_at_x_2_object() (p, cv) == EQUAL);
            return EQUAL;
          }
        }
        else
        {
          CGAL_envelope_voronoi_assertion (hh != Halfedge_handle());
          if (vh->is_on_surfaces(hh->surfaces().begin(), hh->surfaces().end()) && 
              is_bet(cv, p))
          {
            CGAL_envelope_voronoi_assertion(_traits->Base::compare_y_at_x_2_object() (p, cv) == EQUAL);
            return EQUAL;
          }
        }
      }
      else if (p.is_current_point)
      {
        if (cv.is_current_curve)
        {
          CGAL_envelope_voronoi_assertion(_traits->Base::compare_y_at_x_2_object() (p, cv) == EQUAL);
          return EQUAL;
        }
        else
        {
          CGAL_envelope_voronoi_assertion (hh != Halfedge_handle());
          if ((hh->surfaces().find(_traits->_surf1.id) != hh->surfaces().end()) &&
              (hh->surfaces().find(_traits->_surf2.id) != hh->surfaces().end()) && 
              is_bet(cv, p))
          {
            CGAL_envelope_voronoi_assertion(_traits->compare_y_at_x_2_object() (p, cv) == EQUAL);
            return EQUAL;
          }
        }
      }
      
      return _traits->Base::compare_y_at_x_2_object () (p, cv);
    }
    
  };

  Compare_y_at_x_2 compare_y_at_x_2_object () const
  {
    return Compare_y_at_x_2(this);
  }

};

} //namespace CGAL

#endif
