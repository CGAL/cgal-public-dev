//
// C++ Interface: facet_criteria_visitor_with_balls
//
// Description: 
//
//
// Author:  <>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
//--------------DOB----------------------------------------------

#include <CGAL/Mesh_3/mesh_standard_facet_criteria.h>
#include <CGAL/Facet_on_the_same_surface_criterion.h>

namespace CGAL {
namespace Mesh_3 {

template <typename Tr>
  class Facet_criterion_visitor_with_balls
  : public Mesh_3::Criterion_visitor<Tr, typename Tr::Facet>
  {
    typedef Mesh_3::Criterion_visitor<Tr, typename Tr::Facet> Base;
    typedef Facet_criterion_visitor_with_balls<Tr> Self;
   
  public:
    typedef Mesh_3::Abstract_criterion<Tr, Self> Criterion;
	typedef Mesh_3::Curvature_size_criterion<Tr, Self> Curvature_size_criterion;
	typedef Mesh_3::Aspect_ratio_criterion<Tr, Self> Aspect_ratio_criterion;
	typedef Mesh_3::Facet_on_surface_criterion<Tr, Self> Facet_on_surface_criterion;
	typedef Mesh_3::Uniform_size_criterion<Tr, Self> Uniform_size_criterion;
	typedef Mesh_3::Facet_on_the_same_surface_criterion<Tr, Self> Facet_on_the_same_surface_criterion;
	

    typedef typename Base::Quality Facet_quality;
    typedef typename Base::Badness Facet_badness;
    typedef typename Base::Handle Handle;
    typedef Handle Facet;

	typedef typename Tr::Point Point_3;
	typedef typename Tr::Geom_traits::Compute_squared_radius_smallest_orthogonal_sphere_3 
								Squared_radius_orthogonal_sphere;
	int nb_weighted_points;
	std::vector<Point_3> points;
	double radius_ortho_shpere;
	
	//typedef typename Tr::Cell::Surface_index Surface_index;
 	//typedef typename Tr::Vertex_handle Vertex_handle;

   
    // Constructor
    Facet_criterion_visitor_with_balls(const Facet& fh)
    : Base(fh)
    {
		const Point_3 p1 = fh.first->vertex ((fh.second+1)&3)->point();
		const Point_3 p2 = fh.first->vertex ((fh.second+2)&3)->point();
		const Point_3 p3 = fh.first->vertex ((fh.second+3)&3)->point();

		if(p1.weight() > 0) points.push_back(p1);
		if(p2.weight() > 0) points.push_back(p2);
		if(p3.weight() > 0) points.push_back(p3);

		nb_weighted_points = points.size();
		
		if(nb_weighted_points ==3)
			radius_ortho_shpere = Squared_radius_orthogonal_sphere()( points[0], points[1], points[2]);
		else if(nb_weighted_points ==2)
			radius_ortho_shpere = Squared_radius_orthogonal_sphere()( points[0], points[1]);

// 		std::cout<<p1<<std::endl;
// 		std::cout<<p2<<std::endl;
// 		std::cout<<p3<<std::endl;
// 		std::cout<< " nb_weighted_points = "<<nb_weighted_points<<" rad = "<<radius_ortho_shpere<<std::endl;
    }
   
    // Destructor
    ~Facet_criterion_visitor_with_balls() { };
   
//     void visit(const Criterion& criterion)
//     {
//       Base::do_visit(criterion);
//     }
   
    void visit(const Curvature_size_criterion& criterion)
    {
		//Base::increment_counter();
		if ( nb_weighted_points >= 2 && radius_ortho_shpere <= 0.)
			Base::increment_counter();
		else if ( nb_weighted_points == 1)
			Base::increment_counter();
		else
			Base::do_visit(criterion);
	
    }
   
    void visit(const Aspect_ratio_criterion& criterion)
    {
		if ( nb_weighted_points >=2  && radius_ortho_shpere <= 0.)
			Base::increment_counter();
		else if ( nb_weighted_points == 1)
			Base::increment_counter();
		else
			Base::do_visit(criterion);

    }

    void visit(const Facet_on_surface_criterion& criterion)
    {
		if ( nb_weighted_points == 3 && radius_ortho_shpere <= 0.)
			Base::increment_counter();
		else
			Base::do_visit(criterion);
    }

    void visit(const Uniform_size_criterion& criterion)
    {
        if ( nb_weighted_points == 3 && radius_ortho_shpere <= 0.)
			Base::increment_counter();
		else
			Base::do_visit(criterion);
    }

	void visit(const Facet_on_the_same_surface_criterion& criterion)
    {
		if ( nb_weighted_points == 3 && radius_ortho_shpere <= 0.)
			Base::increment_counter();
		else
			Base::do_visit(criterion);
    }

  };  // end class Facet_criterion_visitor


}
}
