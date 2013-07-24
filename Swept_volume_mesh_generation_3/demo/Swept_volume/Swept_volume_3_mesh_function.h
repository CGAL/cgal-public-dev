// Copyright (c) 2010 INRIA Sophia-Antipolis (France).
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
// $URL: svn+ssh://hemmer@scm.gforge.inria.fr/svn/cgal/trunk/Mesh_3/demo/Mesh_3/Mesh_function.h $
// $Id: Mesh_function.h 56883 2010-06-18 16:16:50Z stayeb $
//
//
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description : 
//******************************************************************************

#ifndef SWEEPT_VOLUME_3_MESH_FUNCTION_H
#define SWEEPT_VOLUME_3_MESH_FUNCTION_H

#define CGAL_MESH_3_MESHER_STATUS_ACTIVATED 0

#include <QStringList>
#include <QString>

#include <CGAL/Mesh_3/Mesher_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <Meshing_thread.h>
#include <Mesh_parameters.h>
#include <SV/Facet_criterion_3.h>


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_3/Robust_weighted_circumcenter_filtered_traits_3.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>

#include <SV/Swept_volume_with_vhull_3.h>
#include <SV/Mesh_domain_3.h>

template < typename Domain_ >
class Swept_volume_3_mesh_function
  : public Mesh_function_interface
{
public:
  typedef Domain_ Domain;
  typedef typename CGAL::Mesh_triangulation_3<Domain>::type      Triangulation; 
  typedef CGAL::Mesh_complex_3_in_triangulation_3<Triangulation> C3t3;
  
public:
  Swept_volume_3_mesh_function(C3t3& c3t3, Domain* domain, const Mesh_parameters& p);
  
  ~Swept_volume_3_mesh_function();
  
  // Launch
  virtual void launch();
  void launch(int);
  
  // Stop
  virtual void stop();
  
  // Logs
  virtual QStringList parameters_log() const;
  virtual QString status(double time_period) const;

private:
  typedef typename Domain::Point_3                  Point_3;
  typedef typename Domain::Index                    Index;
  typedef std::vector<std::pair<Point_3, Index> >   Initial_points_vector;
  typedef typename Initial_points_vector::iterator  Ipv_iterator;
  typedef typename C3t3::Vertex_handle                       Vertex_handle;
  
  typedef typename C3t3::Triangulation              Tr;
  typedef CGAL::Mesh_criteria_3<Tr>                 Mesh_criteria;
  typedef typename Mesh_criteria::Facet_criteria    Facet_criteria;
  typedef typename Mesh_criteria::Cell_criteria     Cell_criteria;
  
  typedef CGAL::Mesh_3::Mesher_3<C3t3, Mesh_criteria, Domain>   Mesher;
  
private:
  C3t3& c3t3_;
  Domain* domain_;
  Mesh_parameters p_;
  bool continue_;
  Mesher* mesher_;
  mutable typename Mesher::Mesher_status last_report_;
};



// -----------------------------------
// Class Mesh_parameters
// -----------------------------------
inline
QStringList
Mesh_parameters::
log() const
{
  return QStringList()
  << QString("facet min angle: %1").arg(facet_angle)
  << QString("facet max size: %1").arg(facet_sizing)
  << QString("facet approx error: %1").arg(facet_approx)
  << QString("tet shape (radius-edge): %1").arg(tet_shape)
  << QString("tet max size: %1").arg(tet_sizing);
}


// -----------------------------------
// Class Swept_volume_3_mesh_function
// -----------------------------------
template < typename D_ >
Swept_volume_3_mesh_function<D_>::
Swept_volume_3_mesh_function(C3t3& c3t3, Domain* domain, const Mesh_parameters& p)
: c3t3_(c3t3)
, domain_(domain)
, p_(p)
, continue_(true)
, mesher_(NULL)
, last_report_(0,0,0)
{
  std::cerr << "Swept_volume_3_mesh_function()" << std::endl;  
}


template < typename D_ >
Swept_volume_3_mesh_function<D_>::
~Swept_volume_3_mesh_function()
{
  delete domain_;
  delete mesher_;
}


template < typename D_ >
void
Swept_volume_3_mesh_function<D_>::
launch(){this->launch(-1);}

template < typename D_ >
void
Swept_volume_3_mesh_function<D_>::
launch( int max_vertecies)
{
  std::cerr <<" Swept_volume_3_mesh_function::launch()" <<std::endl; 
 
  
  // Mesh initialization : get some points and add them to the mesh
  Initial_points_vector initial_points;
  domain_->construct_initial_points_object()(std::back_inserter(initial_points),20);
  
  std::cerr <<" Swept_volume_3_mesh_function::launch() initial_points: " << initial_points.size() << std::endl;
  
  // Insert points and set their index and dimension
  for ( Ipv_iterator it = initial_points.begin() ;
       it != initial_points.end() ;
       ++it )
  {
    Vertex_handle v = c3t3_.triangulation().insert(it->first);
    c3t3_.set_dimension(v,2); // by construction, points are on surface
    c3t3_.set_index(v,it->second);
  }
  
  std::cerr <<" Swept_volume_3_mesh_function::launch() vertices : " 
            << c3t3_.triangulation().number_of_vertices ()  << std::endl;

  // Create mesh criteria
  Facet_criteria facet_criteria(p_.facet_angle,p_.facet_sizing,p_.facet_approx);
  Cell_criteria  cell_criteria(p_.tet_shape,p_.tet_sizing);
  
  typedef typename SV::Facet_criterion_3<Tr, typename Facet_criteria::Visitor,typename Domain::Swept_volume_3> SV3_facet_criterion; 
  facet_criteria.add(new SV3_facet_criterion(&(domain_->swept_volume_3())));
  
  Mesh_criteria criteria(facet_criteria, cell_criteria);
  
  // Build mesher and launch refinement process
  mesher_ = new Mesher(c3t3_, *domain_, criteria);
  
  mesher_->initialize(); 
  while ( ! mesher_->is_algorithm_done() && continue_)
  {
    mesher_->one_step(); 
    //if(max_vertecies != -1 && c3t3_.number_of_vertices() >= max_vertecies){
    if(max_vertecies != -1 && c3t3_.number_of_facets () >= max_vertecies){
      continue_ = false;
    }
  }

  // Ensure c3t3 is ok (usefull if process has been stop by the user)
  mesher_->fix_c3t3(); 
}


template < typename D_ >
void
Swept_volume_3_mesh_function<D_>::
stop()
{
  continue_ = false;
}


template < typename D_ >
QStringList
Swept_volume_3_mesh_function<D_>::
parameters_log() const
{
  return p_.log();
}


template < typename D_ >
QString
Swept_volume_3_mesh_function<D_>::
status(double time_period) const
{

  typename Mesher::Mesher_status s = 
    (mesher_ == NULL)? last_report_ : mesher_->status();
  
  QString result = QString("Vertices: %1 \n"
                           "Vertices inserted last %2s: %3 \n\n"
                           "Bad facets: %4 \n"
                           "Bad cells: %5")
    .arg(s.vertices)
    .arg(time_period)
    .arg(s.vertices - last_report_.vertices)
    .arg(s.facet_queue)
    .arg(s.cells_queue);
  
  last_report_ = s;
  
  return result;
}

#endif // SWEEPT_VOLUME_3_MESH_FUNCTION_H
