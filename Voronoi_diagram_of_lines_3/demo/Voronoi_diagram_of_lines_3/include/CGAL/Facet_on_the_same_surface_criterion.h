//
// C++ Interface: Facet_on_the_same_surface_criterion
//
// Description: 
//
//
// Author:  <>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef FACET_ON_THE_SAME_SURFACE_CRITERION_H
#define FACET_ON_THE_SAME_SURFACE_CRITERION_H

#include <CGAL/Mesh_3/mesh_standard_criteria.h>

namespace CGAL {

namespace Mesh_3 {


template <typename Tr, typename Visitor_>
class Facet_on_the_same_surface_criterion :
  public Mesh_3::Abstract_criterion<Tr, Visitor_>
{
private:
  typedef typename Tr::Facet Facet;

  typedef Mesh_3::Abstract_criterion<Tr,Visitor_> Base;
  typedef typename Base::Quality Quality;
  typedef typename Base::Badness Badness;

  typedef Facet_on_the_same_surface_criterion<Tr,Visitor_> Self;

public:
  /// Constructor
  Facet_on_the_same_surface_criterion() {};
  /// Destructor
  ~Facet_on_the_same_surface_criterion() {};

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
    typedef typename Tr::Vertex_handle Vertex_handle;
    typedef typename Tr::Cell_handle Cell_handle;
	typedef typename Tr::Vertex::Index Index;

    const Cell_handle& ch = f.first;
    const int i = f.second;
    Vertex_handle v1 = ch->vertex((i+1)&3);
    Vertex_handle v2 = ch->vertex((i+2)&3);
    Vertex_handle v3 = ch->vertex((i+3)&3);
	

	if ( (v1->in_dimension() == 3) ||
         (v2->in_dimension() == 3) ||
         (v3->in_dimension() == 3) )
    {
		return Badness(Quality(1));
    }
	else
	{
		std::vector<Vertex_handle> points;

 		if(v1->in_dimension() == 2) points.push_back(v1);
 		if(v2->in_dimension() == 2) points.push_back(v2);
 		if(v3->in_dimension() == 2) points.push_back(v3);

 		if ( !points.empty() )
		{
			typename std::vector<Vertex_handle>::iterator vit = points.begin();
			Index index = (*vit)->index();

			while (++vit !=  points.end())
			{
				if(!(index == (*vit)->index()))
				{
					return Badness(Quality(1));
				}
			}
		}
		
		return  Badness();			
	}	
  }
}; // end class Facet_on_the_same_surface_criterion

} // end namespace Mesh_3
} //end namespace CGAL

#endif