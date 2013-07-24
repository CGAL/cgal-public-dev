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

#ifndef SV_FACET_CRITERION_3_H
#define SV_FACET_CRITERION_3_H

#include <CGAL/basic.h>
#include <CGAL/Mesh_3/mesh_standard_criteria.h> 

namespace SV {

// checks whether a facet is valid with respect to a swept volume 
template <typename Tr_, typename Visitor_, typename SweptVolume_3>
class Facet_criterion_3 :
  public CGAL::Mesh_3::Abstract_criterion<Tr_, Visitor_>
{
private:
  typedef SweptVolume_3 Swept_volume_3; 
  typedef Tr_ Tr;
  typedef Visitor_ Visitor;
  
  typedef typename Tr::Facet Facet;
  typedef typename Tr::Geom_traits::FT FT;

  typedef CGAL::Mesh_3::Abstract_criterion<Tr,Visitor_> Base;
  typedef typename Base::Quality Quality;
  typedef typename Base::Badness Badness;

  typedef Facet_criterion_3<Tr,Visitor,Swept_volume_3> Self;
  
  
public: 
  Facet_criterion_3(const Swept_volume_3* sv3):
    p_swept_volume_3(sv3){}

protected:
  virtual void do_accept(Visitor_& v) const
  {
    v.visit(*this);
  }

  virtual Self* do_clone() const
  {
    // Call copy ctor on this
    return new Self(*this);
  }
  
  virtual Badness do_is_bad (const Facet& f) const
  {
    

    typedef typename Tr::Point Point_3;
    const Point_3& p1 = f.first->vertex((f.second+1)&3)->point();
    const Point_3& p2 = f.first->vertex((f.second+2)&3)->point();
    const Point_3& p3 = f.first->vertex((f.second+3)&3)->point();
    
    Badness badness = p_swept_volume_3->compute_facet_badness(p1,p2,p3);
    
    if ( badness )
      return Badness(Quality(*badness));
    else
      return Badness();
  }
  
private:
  const Swept_volume_3* p_swept_volume_3; 

};  // end Facet_criterion_3


} // namespace SV 


#endif // SV_FACET_CRITERION_3_H
