//
// C++ Interface: Mesh_criteria_3_with_balls
//
// Description: 
//
//
// Author:  <>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef MESH_CRITERIA_3_DB_H
#define MESH_CRITERIA_3_DB_H

#include <CGAL/Mesh_facet_criteria_3.h>
#include <CGAL/Mesh_cell_criteria_3.h>
#include <CGAL/facet_criteria_visitor_with_balls.h>
#include <CGAL/Cell_criteria_visitor_with_balls.h>

namespace CGAL {

// Class Mesh_criteria_3
// Provides default meshing criteria to drive Mesh_3 process
template <class Tr>
class Mesh_criteria_3_with_balls
{
public:
  typedef Mesh_facet_criteria_3<Tr,Mesh_3::Facet_criterion_visitor_with_balls<Tr>  >     Facet_criteria;
  typedef Mesh_cell_criteria_3<Tr,Mesh_3::Cell_criteria_visitor_with_balls<Tr>  >      Cell_criteria;

  // Constructor
  Mesh_criteria_3_with_balls(const Facet_criteria& facet_criteria,
                  const Cell_criteria& cell_criteria)
    : facet_criteria_(facet_criteria)
    , cell_criteria_(cell_criteria) { };

  // Destructor
  ~Mesh_criteria_3_with_balls() { };

  const Facet_criteria& facet_criteria() const { return facet_criteria_; };
  const Cell_criteria& cell_criteria() const { return cell_criteria_; };

private:
  Facet_criteria facet_criteria_;
  Cell_criteria cell_criteria_;

};  // end class Mesh_criteria_3


}  // end namespace CGAL


#endif // MESH_CRITERIA_3_H
