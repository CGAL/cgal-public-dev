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
// $URL: 
// $Id: 
//
// Author(s)     : Ophir Setter        <ophir.setter@cs.tau.ac.il
//


#ifndef CGAL_MWA_ARRANGEMENT_2_H
#define CGAL_MWA_ARRANGEMENT_2_H

#include <CGAL/Envelope_3/Envelope_pm_dcel.h>

namespace CGAL {

namespace Min_width_annulus_2
{
  template <class Point_2, class Data>
    class MWA_arr_vertex : public CGAL::Arr_vertex_base<Point_2>,
    public std::pair< Envelope_3::Dcel_info<Data>, Envelope_3::Dcel_info<Data> >
  {
  public:
    typedef CGAL::Arr_vertex_base<Point_2>       Base_vertex;
    typedef std::pair< Envelope_3::Dcel_info<Data>, 
      Envelope_3::Dcel_info<Data> >                          Base_info;
    typedef MWA_arr_vertex<Point_2, Data>        Self;

    /*! Assign from another vertex.
     * \param v the other vertex.
     */
    virtual void assign(const Base_vertex & v)
    {
      Base_vertex::assign(v);
      
      const Self & ex_v = static_cast<const Self&>(v);
      this->Base_info::operator=(ex_v);
    }
  };
  
  template <class X_monotone_curve_2, class Data>
    class MWA_arr_halfedge 
    : public CGAL::Arr_halfedge_base<X_monotone_curve_2>,
    public std::pair< Envelope_3::Dcel_info<Data>, Envelope_3::Dcel_info<Data> >
  {
  public:
    typedef CGAL::Arr_halfedge_base<X_monotone_curve_2>     Base_halfedge;
    typedef std::pair< Envelope_3::Dcel_info<Data>, 
      Envelope_3::Dcel_info<Data> >                                     Base_info;
    typedef MWA_arr_halfedge<X_monotone_curve_2, Data>      Self;
    
    /*! Assign from another halfedge.
     * \param h the other halfedge.
     */
    virtual void assign(const Base_halfedge & h)
    {
      Base_halfedge::assign(h);
      
      const Self & ex_h = static_cast<const Self&>(h);
      this->Base_info::operator=(ex_h);
    }
  };

  template <class Data>
    class MWA_arr_face : public CGAL::Arr_face_base,
    public std::pair< Envelope_3::Dcel_info<Data>, Envelope_3::Dcel_info<Data> >
  {
  public:
    typedef CGAL::Arr_face_base                             Base_face;
    typedef std::pair< Envelope_3::Dcel_info<Data>, 
      Envelope_3::Dcel_info<Data> >                                     Base_info;
    typedef MWA_arr_face<Data>                              Self;

    typedef typename Envelope_3::Dcel_info<Data>::Data_container        Data_container;
    typedef typename Data_container::const_iterator         Data_const_iterator;

    /*! Assign from another face.
     * \param f the other face.
     */
    virtual void assign (const Base_face & f)
    {
      Base_face::assign(f);
      
      const Self & ex_f = static_cast<const Self&>(f);
      this->Base_info::operator=(ex_f);
    }
  };

  template <class Traits, class Data>
    class MWA_arr_dcel 
    : public CGAL::Arr_dcel_base<MWA_arr_vertex<typename Traits::Point_2, Data>,
    MWA_arr_halfedge<typename Traits::X_monotone_curve_2, Data>,
    MWA_arr_face<Data> >
    {
      typedef Data                                    Face_data;
      typedef typename MWA_arr_face<Data>::Data_const_iterator
        Face_data_const_iterator;
    public: 
      typedef Face_data_const_iterator                Dcel_data_const_iterator;
    };

}

template <class Traits_> 
class MWA_arrangement_2
: public Arrangement_2< Traits_,
  Min_width_annulus_2::MWA_arr_dcel< Traits_,
  typename Traits_::Xy_monotone_surface_3> >
{  
public:
  typedef MWA_arrangement_2< Traits_ >                  Self;
  typedef Arrangement_2< Traits_,
    Min_width_annulus_2::MWA_arr_dcel< Traits_,
    typename Traits_::Xy_monotone_surface_3> >          Base;
  typedef typename Base::Dcel                           Base_dcel;

  typedef typename Base_dcel::Dcel_data_const_iterator  Surface_const_iterator;
};

} //namespace CGAL

#endif // CGAL_MWA_ARRANGEMENT_2_H
