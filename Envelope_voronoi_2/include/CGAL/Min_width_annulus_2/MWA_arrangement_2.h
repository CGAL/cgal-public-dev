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

#ifndef CGAL_MWA_ARRANGEMENT_2_H
#define CGAL_MWA_ARRANGEMENT_2_H

#include <CGAL/Envelope_3/Envelope_pm_dcel.h>

namespace CGAL {

namespace Min_width_annulus_2 {
  template <typename Point_2, class Data>
    class MWA_arr_vertex : public CGAL::Arr_vertex_base<Point_2>,
                           public std::pair<Envelope_3::Dcel_info<Data>,
                                            Envelope_3::Dcel_info<Data>> {
  public:
    using Base_vertex = CGAL::Arr_vertex_base<Point_2>;
    using Base_info =
      std::pair< Envelope_3::Dcel_info<Data>, Envelope_3::Dcel_info<Data>>;
    using Self = MWA_arr_vertex<Point_2, Data>;

    /*! Assign from another vertex.
     * \param v the other vertex.
     */
    virtual void assign(const Base_vertex& v) {
      Base_vertex::assign(v);
      const Self & ex_v = static_cast<const Self&>(v);
      this->Base_info::operator=(ex_v);
    }
  };

  template <typename X_monotone_curve_2, class Data>
    class MWA_arr_halfedge : public CGAL::Arr_halfedge_base<X_monotone_curve_2>,
    public std::pair< Envelope_3::Dcel_info<Data>, Envelope_3::Dcel_info<Data>>
  {
  public:
    using Base_halfedge = CGAL::Arr_halfedge_base<X_monotone_curve_2>;
    using Base_info = std::pair<Envelope_3::Dcel_info<Data>,
                                Envelope_3::Dcel_info<Data>>;
    using Self = MWA_arr_halfedge<X_monotone_curve_2, Data>;

    /*! Assign from another halfedge.
     * \param h the other halfedge.
     */
    virtual void assign(const Base_halfedge& h) {
      Base_halfedge::assign(h);

      const Self & ex_h = static_cast<const Self&>(h);
      this->Base_info::operator=(ex_h);
    }
  };

  template <typename Data>
    class MWA_arr_face : public CGAL::Arr_face_base,
                         public std::pair< Envelope_3::Dcel_info<Data>,
                                           Envelope_3::Dcel_info<Data>> {
  public:
    using Base_face = CGAL::Arr_face_base;
    using Base_info = std::pair< Envelope_3::Dcel_info<Data>,
                                 Envelope_3::Dcel_info<Data>>;
    using Self = MWA_arr_face<Data>;

    using Data_container = typename Envelope_3::Dcel_info<Data>::Data_container;
    using Data_const_iterator = typename Data_container::const_iterator;

    /*! Assign from another face.
     * \param f the other face.
     */
    virtual void assign (const Base_face& f) {
      Base_face::assign(f);

      const Self & ex_f = static_cast<const Self&>(f);
      this->Base_info::operator=(ex_f);
    }
  };

  template <typename Traits, typename Data>
  class MWA_arr_dcel :
    public CGAL::Arr_dcel_base<MWA_arr_vertex<typename Traits::Point_2, Data>,
                               MWA_arr_halfedge<typename Traits::X_monotone_curve_2, Data>,
                               MWA_arr_face<Data>> {
    using Face_data = Data;
    using Face_data_const_iterator =
      typename MWA_arr_face<Data>::Data_const_iterator;

  public:
    using Dcel_data_const_iterator = Face_data_const_iterator;
  };
}

template <typename Traits_>
class MWA_arrangement_2 : public Arrangement_2< Traits_,
  Min_width_annulus_2::MWA_arr_dcel<Traits_,
                                    typename Traits_::Xy_monotone_surface_3>> {
public:
  using Self = MWA_arrangement_2< Traits_ >;
  using Base = Arrangement_2<Traits_,
                             Min_width_annulus_2::MWA_arr_dcel
                             <Traits_, typename Traits_::Xy_monotone_surface_3>>;
  using Base_dcel = typename Base::Dcel;

  using Surface_const_iterator = typename Base_dcel::Dcel_data_const_iterator;
};

} //namespace CGAL

#endif // CGAL_MWA_ARRANGEMENT_2_H
