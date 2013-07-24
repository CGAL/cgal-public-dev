// Copyright (c) 2011 Andreas von Dziegielewski and Michael Hemmer (Germany).
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
// $URL:$
// $Id:$
//
// Author(s)     : Michael Hemmer (mhsaar@googlemail.com)
//
// ================================================================================

#ifndef SV_MESH_DOMAIN_3_H
#define SV_MESH_DOMAIN_3_H

#include <CGAL/basic.h>
#include <CGAL/Implicit_mesh_domain_3.h>
#include <CGAL/Surface_mesher/Implicit_surface_oracle_3.h>

// #include <Swept_volume_3_type.h>

namespace SV{

// For now I derive from Implicit_function 
template< class SweptVolume_3 >
class Mesh_domain_3
  :public CGAL::Implicit_mesh_domain_3<SweptVolume_3,typename SweptVolume_3::Kernel>{

public:
  typedef SweptVolume_3 Swept_volume_3; 
  typedef Mesh_domain_3<Swept_volume_3> Self;
  typedef typename Swept_volume_3::Kernel Kernel; 
  typedef typename Kernel::Point_3 Point_3; 
  typedef typename Kernel::FT FT; 
  typedef CGAL::Implicit_mesh_domain_3<Swept_volume_3,Kernel> Base; 
  typedef typename Base::Intersection Intersection;
  typedef typename Base::Segment_3 Segment_3;
  typedef typename Base::Surface_index Surface_index;
  typedef typename Base::Index Index;

public: 
  // additional typedef for Surface_mesh_traits; 
  typedef Swept_volume_3 Surface_3;
  typedef Point_3 Intersection_point; 
  
public:
  Mesh_domain_3(Swept_volume_3& sv3)
    :m_swept_volume_3(sv3),Base(sv3,sv3.bounding_sphere(),FT(1e-3)){
  }
      
  Swept_volume_3& m_swept_volume_3; 
  Swept_volume_3& swept_volume_3(){return m_swept_volume_3;}
  const Swept_volume_3& swept_volume_3() const {return m_swept_volume_3;}
  
  std::vector<Point_3> m_initial_points; 
  std::vector<Point_3>& initial_points(){return m_initial_points;}
  const std::vector<Point_3>& initial_points() const {return m_initial_points;}
  

#if SV_HAS_INTERSECT 
  struct Construct_intersection
  {
  private:
    const Mesh_domain_3<Swept_volume_3>& m_domain; 
  public:
    Construct_intersection(const Mesh_domain_3<Swept_volume_3>& domain)
      : m_domain(domain) {}
    
    // Query may be an CGAL Line/Ray/Segment 
    template <class Query> 
    Intersection operator()(const Query& query) const
    {
      // clip Query to segment, if empty return no intersection 
      boost::optional<Segment_3> op_segment = clip_to_segment(query);
      if(!op_segment) return Intersection();
      
      // get intersection index if empty return no intersection 
      boost::optional<Surface_index>  op_index = m_domain.do_intersect_surface_object()(*op_segment);
      if(!op_index) return Intersection();  

      Point_3 p = m_domain.swept_volume_3().intersect(op_segment->source(),op_segment->target());
      return Intersection(p,*op_index,2);

      return  m_domain.Base::construct_intersection_object()(query);

    }

    /// Clips \c query to a segment \c s or nothing 
    template<typename Query>
    boost::optional<Segment_3> clip_to_segment(const Query& query) const
    {
      const CGAL::Object clipped = CGAL::intersection(query, m_domain.bounding_box());
      if ( const Segment_3* s = CGAL::object_cast<Segment_3>(&clipped) )
        return  boost::optional<Segment_3>(*s);
      else
        return boost::optional<Segment_3>();
    }
  };

  /// Returns Construct_intersection object
  Construct_intersection construct_intersection_object() const
  {
    return Construct_intersection(*this);
  }
#endif // SV_HAS_INTERSECT


  struct Construct_initial_points {
  private:
    const Mesh_domain_3<Swept_volume_3>& m_domain; 
  public:
    Construct_initial_points(const Mesh_domain_3<Swept_volume_3>& domain)
      : m_domain(domain) {}
    
    template <class OutputIterator> 
    OutputIterator operator()(OutputIterator oit, int n = 20) const {
      Index index(Surface_index(-1,1));
      
      BOOST_FOREACH(Point_3 point, m_domain.initial_points()){
        *oit++ = std::make_pair(point,index);
      }
      
      if(n>m_domain.initial_points().size()){
        m_domain.Base::construct_initial_points_object()(
            oit,n-m_domain.initial_points().size());
      }
          
      return oit; 
    }
    
    // dummy operator to make it a surface_mesh_traits_3
    // 
    template <class Surface_3, class OutputIterator> 
    OutputIterator operator()(const Surface_3&, OutputIterator oit, int n = 20) const{
      // OutputIterator takes just points not point index pairs
      typedef std::pair<Point_3,Index> PIX;
      std::vector<PIX> pixs; 
      this->operator()(std::back_inserter(pixs),n); 
      BOOST_FOREACH(const PIX& pix, pixs){
        *oit++ = pix.first; 
      }
      return oit;
    }
  };
  /// Returns Construct_initial_points object
  Construct_initial_points construct_initial_points_object() const
  {
    return Construct_initial_points(*this);
  }
  
  // Intersect_3 for Surface_mesh_3 
  
  class Intersect_3{
  private:
    const Mesh_domain_3<Swept_volume_3>& m_domain; 
  public:
    Intersect_3(const Mesh_domain_3<Swept_volume_3>& domain)
      : m_domain(domain) {}
 
  public: 
    template <class Query> 
    CGAL::Object operator()(const Surface_3&, const Query& query){
      if(!m_domain.do_intersect_surface_object()(query)){
        return CGAL::Object();
      }
      // otherwise intersect: 
      Intersection intersection = 
        m_domain.construct_intersection_object()(query);

      // take the point 
      return CGAL::make_object(CGAL::cpp0x::get<0>(intersection));
    }
  };
  
  Intersect_3 intersect_3_object() const
  {
    return Intersect_3(*this);
  }



private:
  // Disabled copy constructor & assignment operator
  Mesh_domain_3();
  Mesh_domain_3(const Self& src);
  Self& operator=(const Self& src);  
};

} // namespace SV 

// make it a Surface_mesh_traits;
namespace SV{
template <class Kernel> class Swept_volume_with_vhull_3; 
}
namespace CGAL{
template <class Kernel> class Surface_mesh_traits_generator_3; 
template <class Kernel> 
class Surface_mesh_traits_generator_3<
SV::Swept_volume_with_vhull_3<Kernel> >{
public:
  typedef SV::Swept_volume_with_vhull_3<Kernel> Surface_3; 
  typedef SV::Mesh_domain_3<Surface_3> type; 
};
}

#endif // SV_MESH_DOMAIN_3_H
