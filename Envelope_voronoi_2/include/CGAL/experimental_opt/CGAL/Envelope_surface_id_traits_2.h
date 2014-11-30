// Copyright (c) 2008  Tel-Aviv University (Israel).
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
// Author(s)     : Ophir Setter       <ophir.setter@cs.tau.ac.il>


#ifndef CGAL_ENV_SURFACE_ID_TRAITS_2_H
#define CGAL_ENV_SURFACE_ID_TRAITS_2_H

#include <CGAL/basic.h>
#include <CGAL/tags.h>
#include <CGAL/functions_on_enums.h>
#include <CGAL/number_utils.h> 
#include <CGAL/Envelope_3/Envelope_base.h>

#include <CGAL/utility.h>

namespace CGAL {

/*!
 * \class The traits class
 */
template <class Traits_>
class Envelope_surface_id_traits_2 : public Traits_
{
 public:
  typedef Traits_                                               Traits;
  typedef Envelope_surface_id_traits_2<Traits>                          Self;
  typedef Traits                                                Base;
  typedef typename Base::Xy_monotone_surface_3                  Base_xy_monotone_surface_3;
  
  typedef typename Base::Surface_3                              Surface_3;

  typedef unsigned int                                          Surface_id;

  class Xy_monotone_surface_3 : public Base_xy_monotone_surface_3
  {
  public:
    Xy_monotone_surface_3(const Base_xy_monotone_surface_3& base, Surface_id in_id)
      : Base_xy_monotone_surface_3(base), id(in_id) 
      {}

    Surface_id id;
  };

protected:
  Surface_id                                                   _current_id;

public:

  Envelope_surface_id_traits_2<Base>() : _current_id(1)
  {
  }
  
  Surface_id get_id()
  {
    return _current_id++;
  }

  class Make_xy_monotone_3
  {
  private:
    Self &_traits;

  public:
    Make_xy_monotone_3(Self &traits) : _traits(traits) {}
    
    template <class OutputIterator>
      OutputIterator operator()(const Surface_3& s,
                                bool is_lower,
                                OutputIterator o) 
    {
      Surface_id k = _traits.get_id();
      
      typedef std::list<Base_xy_monotone_surface_3> Surface_list;
      Surface_list l;
      _traits.Base::make_xy_monotone_3_object() (s, is_lower, std::back_inserter(l));
      
      for (typename Surface_list::iterator i = l.begin(); i != l.end(); ++i)
      {
        *o++ = Xy_monotone_surface_3(*i, k);
      }
      return o;
    }
  };
  
  Make_xy_monotone_3 make_xy_monotone_3_object()
  {
    return Make_xy_monotone_3(*this);
  }
};

} //namespace CGAL

#endif // CGAL_ENV_NUMBER_XY_TRAITS_H
