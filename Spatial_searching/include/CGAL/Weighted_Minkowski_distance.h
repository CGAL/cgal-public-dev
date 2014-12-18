// Copyright (c) 2002,2011 Utrecht University (The Netherlands).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Hans Tangelder (<hanst@cs.uu.nl>)

// Note: Use p=0 to denote the weighted Linf-distance 
// For 0<p<1 Lp is not a metric

#ifndef CGAL_WEIGHTED_MINKOWSKI_DISTANCE_H
#define CGAL_WEIGHTED_MINKOWSKI_DISTANCE_H

#include <cmath>
#include <vector>

#include <CGAL/number_utils.h>
#include <CGAL/Kd_tree_rectangle.h>
#include <CGAL/Dimension.h>

namespace CGAL {

  template <class SearchTraits, class Query_type = typename SearchTraits::Point_d>
  class Weighted_Minkowski_distance;

  namespace internal{
	    #ifndef HAS_DIMENSION_TAG
		#define HAS_DIMENSION_TAG
		BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(has_dimension,Dimension,false);
	    #endif

	  template <class SearchTraits, bool has_dim = has_dimension<SearchTraits>::value>
	  struct Weighted_Minkowski_distance_base;

	  template <class SearchTraits>
	  struct Weighted_Minkowski_distance_base<SearchTraits,true>{
		  typedef typename SearchTraits::Dimension Dimension;
	  };

	  template <class SearchTraits>
	  struct Weighted_Minkowski_distance_base<SearchTraits,false>{
		  typedef Dynamic_dimension_tag Dimension;
	  };
   }

  template <class SearchTraits, class Query_type>
  class Weighted_Minkowski_distance {
    SearchTraits traits;
    public:

    typedef typename SearchTraits::Point_d Point_d;
    typedef Point_d                        Query_item;
    typedef typename SearchTraits::FT      FT;
    typedef std::vector<FT>                Weight_vector;
    typedef typename internal::Euclidean_distance_base<SearchTraits>::Dimension Dimension;

    private:

    typedef typename SearchTraits::Cartesian_const_iterator_d Coord_iterator;
    FT power; 

    Weight_vector the_weights;

    public:


    // default constructor
    Weighted_Minkowski_distance(const SearchTraits& traits_=SearchTraits())
      : traits(traits_),power(2) 
    {}

    Weighted_Minkowski_distance(const int d,const SearchTraits& traits_=SearchTraits()) 
      : traits(traits_),power(2), the_weights(d)
    {
      for (int i = 0; i < d; ++i) the_weights[i]=FT(1);
    }

    //default copy constructor and destructor
    

    Weighted_Minkowski_distance (FT pow, int dim,
				 const Weight_vector& weights,
                                 const SearchTraits& traits_=SearchTraits()) 
      : traits(traits_),power(pow)
    {
      CGAL_assertion(power >= FT(0));
      CGAL_assertion(dim==weights.size());
      for (unsigned int i = 0; i < weights.size(); ++i)
	CGAL_assertion(weights[i]>=FT(0));
      the_weights.resize(weights.size());
      the_weights = weights;
    }

    template <class InputIterator>
    Weighted_Minkowski_distance (FT pow, int dim,
				 InputIterator begin, InputIterator end,
                                 const SearchTraits& traits_=SearchTraits()) 
      : traits(traits_),power(pow)
    {
      CGAL_assertion(power >= FT(0));
      the_weights.resize(dim);
      std::copy(begin, end, the_weights.begin());
      for (int i = 0; i < dim; ++i){
	the_weights[i] = *begin;
	++begin;
	CGAL_assertion(the_weights[i]>=FT(0));
      }
      CGAL_assertion(begin == end);
    }


    inline FT transformed_distance(const Query_item& q, const Point_d& p) const {
        return transformed_distance(q,p, Dimension());
    }

    //Dynamic version for runtime dimension
    inline 
    FT 
    transformed_distance(const Query_item& q, const Point_d& p, Dynamic_dimension_tag dt) const
    {
      FT distance = FT(0);
      typename SearchTraits::Construct_cartesian_const_iterator_d construct_it=
        traits.construct_cartesian_const_iterator_d_object();
      Coord_iterator qit = construct_it(q),
	             qe = construct_it(q,1), 
	             pit = construct_it(p);
      if (power == FT(0)) {
	for (unsigned int i = 0; qit != qe; ++qit, ++i)
	  if (the_weights[i] * CGAL::abs((*qit) - (*pit)) > distance)
	    distance = the_weights[i] * CGAL::abs((*qit)-(*pit));
      }
      else
	for (unsigned int i = 0; qit != qe; ++qit, ++i)
	  distance += 
	    the_weights[i] * std::pow(CGAL::abs((*qit)-(*pit)),power);
      return distance;
    }

    //Generic version for DIM > 3
    template <int DIM>
    inline FT 
    transformed_distance(const Query_item& q, const Point_d& p, Dimension_tag<DIM> dt) const
    {
      FT distance = FT(0);
      typename SearchTraits::Construct_cartesian_const_iterator_d construct_it=
        traits.construct_cartesian_const_iterator_d_object();
      Coord_iterator qit = construct_it(q),
	             qe = construct_it(q,1), 
	             pit = construct_it(p);
      if (power == FT(0)) {
	for (unsigned int i = 0; qit != qe; ++qit, ++i)
	  if (the_weights[i] * CGAL::abs((*qit) - (*pit)) > distance)
	    distance = the_weights[i] * CGAL::abs((*qit)-(*pit));
      }
      else
	for (unsigned int i = 0; qit != qe; ++qit, ++i)
	  distance += 
	    the_weights[i] * std::pow(CGAL::abs((*qit)-(*pit)),power);
      return distance;
    }

    //DIM = 2 loop unrolled
    inline FT 
    transformed_distance(const Query_item& q, const Point_d& p, Dimension_tag<2> dt) const
    {
      FT distance = FT(0);
      typename SearchTraits::Construct_cartesian_const_iterator_d construct_it=
        traits.construct_cartesian_const_iterator_d_object();
      Coord_iterator qit = construct_it(q),
	             pit = construct_it(p);
      if (power == FT(0)) {
	  if (the_weights[0] * CGAL::abs((*qit) - (*pit)) > distance)
	    distance = the_weights[0] * CGAL::abs((*qit)-(*pit));
          qit++;pit++;
          if (the_weights[1] * CGAL::abs((*qit) - (*pit)) > distance)
	    distance = the_weights[1] * CGAL::abs((*qit)-(*pit));
      }
      else{
	  distance += 
	    the_weights[0] * std::pow(CGAL::abs((*qit)-(*pit)),power);
          qit++;pit++;
          distance += 
	    the_weights[1] * std::pow(CGAL::abs((*qit)-(*pit)),power);
      }
      return distance;
    }

    //DIM = 3 loop unrolled
    inline FT 
    transformed_distance(const Query_item& q, const Point_d& p, Dimension_tag<3> dt) const
    {
      FT distance = FT(0);
      typename SearchTraits::Construct_cartesian_const_iterator_d construct_it=
        traits.construct_cartesian_const_iterator_d_object();
      Coord_iterator qit = construct_it(q),
	             pit = construct_it(p);
      if (power == FT(0)) {
	  if (the_weights[0] * CGAL::abs((*qit) - (*pit)) > distance)
	    distance = the_weights[0] * CGAL::abs((*qit)-(*pit));
          qit++;pit++;
          if (the_weights[1] * CGAL::abs((*qit) - (*pit)) > distance)
	    distance = the_weights[1] * CGAL::abs((*qit)-(*pit));
          qit++;pit++;
          if (the_weights[2] * CGAL::abs((*qit) - (*pit)) > distance)
	    distance = the_weights[2] * CGAL::abs((*qit)-(*pit));
      }
      else{
	  distance += 
	    the_weights[0] * std::pow(CGAL::abs((*qit)-(*pit)),power);
          qit++;pit++;
          distance += 
	    the_weights[1] * std::pow(CGAL::abs((*qit)-(*pit)),power);
          qit++;pit++;
          distance += 
	    the_weights[2] * std::pow(CGAL::abs((*qit)-(*pit)),power);
      }
      return distance;
    }

    inline 
    FT 
    min_distance_to_rectangle(const Query_item& q,
			      const Kd_tree_rectangle<FT,Dimension>& r) const 
    {
      FT distance = FT(0);
      typename SearchTraits::Construct_cartesian_const_iterator_d construct_it=
        traits.construct_cartesian_const_iterator_d_object();
      Coord_iterator qit = construct_it(q), qe = construct_it(q,1);
      if (power == FT(0))
	{
	  for (unsigned int i = 0; qit != qe; ++qit, ++i) {
	    if (the_weights[i]*(r.min_coord(i) - 
				(*qit)) > distance)
	      distance = the_weights[i] * (r.min_coord(i)-
					   (*qit));
	    if (the_weights[i] * ((*qit) - r.max_coord(i)) > 
		distance)
	      distance = the_weights[i] * 
		((*qit)-r.max_coord(i));
	  }
	}
      else
	{
	  for (unsigned int i = 0; qit != qe; ++qit, ++i) {
	    if ((*qit) < r.min_coord(i))
	      distance += the_weights[i] * 
		std::pow(r.min_coord(i)-(*qit),power);
	    if ((*qit) > r.max_coord(i))
	      distance += the_weights[i] * 
		std::pow((*qit)-r.max_coord(i),power);
	  }
	};
      return distance;
    }

    inline 
    FT 
    min_distance_to_rectangle(const Query_item& q,
			      const Kd_tree_rectangle<FT,Dimension>& r,std::vector<FT>& dists) {
      FT distance = FT(0);
      typename SearchTraits::Construct_cartesian_const_iterator_d construct_it=
        traits.construct_cartesian_const_iterator_d_object();
      Coord_iterator qit = construct_it(q), qe = construct_it(q,1);
      if (power == FT(0))
	{
	  for (unsigned int i = 0; qit != qe; ++qit, ++i) {
	    if (the_weights[i]*(r.min_coord(i) - 
				(*qit)) > distance){
              dists[i] = (r.min_coord(i)-
		(*qit));
	      distance = the_weights[i] * dists[i];
            }
	    if (the_weights[i] * ((*qit) - r.max_coord(i)) > 
		distance){
                  dists[i] = 
		((*qit)-r.max_coord(i));
	      distance = the_weights[i] * dists[i];
            }
	  }
	}
      else
	{
	  for (unsigned int i = 0; qit != qe; ++qit, ++i) {
	    if ((*qit) < r.min_coord(i)){
              dists[i] = r.min_coord(i)-(*qit);
	      distance += the_weights[i] * 
		std::pow(dists[i],power);
            }
	    if ((*qit) > r.max_coord(i)){
              dists[i] = (*qit)-r.max_coord(i);
	      distance += the_weights[i] * 
		std::pow(dists[i],power);
            }
	  }
	};
      return distance;
    }

    inline 
    FT
    max_distance_to_rectangle(const Query_item& q,
			      const Kd_tree_rectangle<FT,Dimension>& r) const {
      FT distance=FT(0);
      typename SearchTraits::Construct_cartesian_const_iterator_d construct_it=
        traits.construct_cartesian_const_iterator_d_object();
      Coord_iterator qit = construct_it(q), qe = construct_it(q,1);
      if (power == FT(0))
	{
	  for (unsigned int i = 0; qit != qe; ++qit, ++i) {
	    if ((*qit) >= (r.min_coord(i) + 
			 r.max_coord(i))/FT(2.0)) {
	      if (the_weights[i] * ((*qit) - 
				    r.min_coord(i)) > distance)
		distance = the_weights[i] * 
		  ((*qit)-r.min_coord(i));
	      else
		if (the_weights[i] * 
		    (r.max_coord(i) - (*qit)) > distance)
		  distance = the_weights[i] * 
		    ( r.max_coord(i)-(*qit));
            }
	  }
	}
      else
	{
	  for (unsigned int i = 0; qit != qe; ++qit, ++i) {
	    if ((*qit) <= (r.min_coord(i)+r.max_coord(i))/FT(2.0))
	      distance += the_weights[i] * std::pow(r.max_coord(i)-(*qit),power);
	    else
	      distance += the_weights[i] * std::pow((*qit)-r.min_coord(i),power);
	  }
	};
      return distance;
    }

     inline 
    FT
    max_distance_to_rectangle(const Query_item& q,
			      const Kd_tree_rectangle<FT,Dimension>& r,std::vector<FT>& dists) {
      FT distance=FT(0);
      typename SearchTraits::Construct_cartesian_const_iterator_d construct_it=
        traits.construct_cartesian_const_iterator_d_object();
      Coord_iterator qit = construct_it(q), qe = construct_it(q,1);
      if (power == FT(0))
	{
	  for (unsigned int i = 0; qit != qe; ++qit, ++i) {
	    if ((*qit) >= (r.min_coord(i) + 
			 r.max_coord(i))/FT(2.0)) {
	      if (the_weights[i] * ((*qit) - 
				    r.min_coord(i)) > distance){
                dists[i] = (*qit)-r.min_coord(i);
		distance = the_weights[i] * 
		  (dists[i]);
              }
	      else
		if (the_weights[i] * 
		    (r.max_coord(i) - (*qit)) > distance){
                      dists[i] =  r.max_coord(i)-(*qit);
		  distance = the_weights[i] * 
		    (dists[i]);
                }
            }
	  }
	}
      else
	{
	  for (unsigned int i = 0; qit != qe; ++qit, ++i) {
	    if ((*qit) <= (r.min_coord(i)+r.max_coord(i))/FT(2.0)){
              dists[i] = r.max_coord(i)-(*qit);
	      distance += the_weights[i] * std::pow(dists[i],power);
            }
	    else{
              dists[i] = (*qit)-r.min_coord(i);
	      distance += the_weights[i] * std::pow(dists[i],power);
            }
	  }
	};
      return distance;
    }
    
    inline 
    FT 
    new_distance(FT dist, FT old_off, FT new_off,
		 int cutting_dimension)  const 
    {
      FT new_dist;
      if (power == FT(0))
	{
	  if (the_weights[cutting_dimension]*CGAL::abs(new_off) 
	      > dist) 
	    new_dist= 
	      the_weights[cutting_dimension]*CGAL::abs(new_off);
	  else new_dist=dist;
	}
      else
	{
	  new_dist = dist + the_weights[cutting_dimension] * 
	    (std::pow(CGAL::abs(new_off),power)-std::pow(CGAL::abs(old_off),power));
	}
      return new_dist;
    }
    
    inline 
    FT 
    transformed_distance(FT d) const 
    {
      if (power <= FT(0)) return d;
      else return std::pow(d,power);
      
    }
    
    inline 
    FT 
    inverse_of_transformed_distance(FT d) const 
    {
      if (power <= FT(0)) return d;
      else return std::pow(d,1/power);
      
    }

  }; // class Weighted_Minkowski_distance

  //Partial specialization for rectangles
  template <class SearchTraits>
  class Weighted_Minkowski_distance<SearchTraits, typename SearchTraits::Iso_box_d> {
    SearchTraits traits;
    public:

    typedef typename SearchTraits::Point_d Point_d;
    typedef typename SearchTraits::Iso_box_d Query_item;
    typedef typename SearchTraits::FT      FT;
    typedef std::vector<FT>                Weight_vector;
    typedef typename internal::Euclidean_distance_base<SearchTraits>::Dimension Dimension;

    private:

    typedef typename SearchTraits::Cartesian_const_iterator_d Coord_iterator;
    FT power; 

    Weight_vector the_weights;

    public:


    // default constructor
    Weighted_Minkowski_distance(const SearchTraits& traits_=SearchTraits())
      : traits(traits_),power(2) 
    {}

    Weighted_Minkowski_distance(const int d,const SearchTraits& traits_=SearchTraits()) 
      : traits(traits_),power(2), the_weights(d)
    {
      for (int i = 0; i < d; ++i) the_weights[i]=FT(1);
    }

    //default copy constructor and destructor
    

    Weighted_Minkowski_distance (FT pow, int dim,
				 const Weight_vector& weights,
                                 const SearchTraits& traits_=SearchTraits()) 
      : traits(traits_),power(pow)
    {
      CGAL_assertion(power >= FT(0));
      CGAL_assertion(dim==weights.size());
      for (unsigned int i = 0; i < weights.size(); ++i)
	CGAL_assertion(weights[i]>=FT(0));
      the_weights.resize(weights.size());
      the_weights = weights;
    }

    template <class InputIterator>
    Weighted_Minkowski_distance (FT pow, int dim,
				 InputIterator begin, InputIterator end,
                                 const SearchTraits& traits_=SearchTraits()) 
      : traits(traits_),power(pow)
    {
      CGAL_assertion(power >= FT(0));
      the_weights.resize(dim);
      std::copy(begin, end, the_weights.begin());
      for (int i = 0; i < dim; ++i){
	the_weights[i] = *begin;
	++begin;
	CGAL_assertion(the_weights[i]>=FT(0));
      }
      CGAL_assertion(begin == end);
    }


    inline FT transformed_distance(const Query_item& q, const Point_d& p) const {
        return transformed_distance(q,p, Dimension());
    }

    //Dynamic version for runtime dimension
    inline 
    FT 
    transformed_distance(const Query_item& q, const Point_d& p, Dynamic_dimension_tag dt) const
    {
      FT distance = FT(0);
      typename SearchTraits::Construct_cartesian_const_iterator_d construct_it=
          traits.construct_cartesian_const_iterator_d_object();
      typename SearchTraits::Construct_min_vertex_d construct_min_vertex;
      typename SearchTraits::Construct_max_vertex_d construct_max_vertex;
      typename SearchTraits::Cartesian_const_iterator_d qmaxit = construct_it(construct_max_vertex(q)),
	  qe = construct_it(construct_max_vertex(q),1), qminit = construct_it(construct_min_vertex(q)),
	  pit = construct_it(p);
      if (power == FT(0)) {
	for (unsigned int i = 0; qmaxit != qe; ++pit,++qmaxit,++qminit)
          FT dist = 0;
          if ((*pit)>(*qmaxit)) dist = ((*pit)-(*qmaxit)); 
	  else if ((*pit)<(*qminit)) dist = ((*qminit)-(*pit));

	  if (the_weights[i] * dist > distance)
	    distance = the_weights[i] * dist;
      }
      else{
	for (unsigned int i = 0; qmaxit != qe; ++pit,++qmaxit,++qminit)
          if ((*pit)>(*qmaxit)) distance += 
			the_weights[i] * std::pow((*pit)-(*qmaxit),power);
	  else if ((*pit)<(*qminit)) distance += 
			the_weights[i] * std::pow((*qminit)-(*pit),power);
      }
      return distance;
    }

    //Generic version for DIM > 3
    template <int DIM>
    inline FT 
    transformed_distance(const Query_item& q, const Point_d& p, Dimension_tag<DIM> dt) const
    {
       FT distance = FT(0);
      typename SearchTraits::Construct_cartesian_const_iterator_d construct_it=
          traits.construct_cartesian_const_iterator_d_object();
      typename SearchTraits::Construct_min_vertex_d construct_min_vertex;
      typename SearchTraits::Construct_max_vertex_d construct_max_vertex;
      typename SearchTraits::Cartesian_const_iterator_d qmaxit = construct_it(construct_max_vertex(q)),
	  qe = construct_it(construct_max_vertex(q),1), qminit = construct_it(construct_min_vertex(q)),
	  pit = construct_it(p);
      if (power == FT(0)) {
	for (unsigned int i = 0; qmaxit != qe; ++pit,++qmaxit,++qminit)
          FT dist = 0;
          if ((*pit)>(*qmaxit)) dist = ((*pit)-(*qmaxit)); 
	  else if ((*pit)<(*qminit)) dist = ((*qminit)-(*pit));

	  if (the_weights[i] * dist > distance)
	    distance = the_weights[i] * dist;
      }
      else{
	for (unsigned int i = 0; qmaxit != qe; ++pit,++qmaxit,++qminit)
          if ((*pit)>(*qmaxit)) distance += 
			the_weights[i] * std::pow((*pit)-(*qmaxit),power);
	  else if ((*pit)<(*qminit)) distance += 
			the_weights[i] * std::pow((*qminit)-(*pit),power);
      }
      return distance;
    }

    //DIM = 2 loop unrolled
    inline FT 
    transformed_distance(const Query_item& q, const Point_d& p, Dimension_tag<2> dt) const
    {
       FT distance = FT(0);
     typename SearchTraits::Construct_cartesian_const_iterator_d construct_it=
          traits.construct_cartesian_const_iterator_d_object();
      typename SearchTraits::Construct_min_vertex_d construct_min_vertex;
      typename SearchTraits::Construct_max_vertex_d construct_max_vertex;
      typename SearchTraits::Cartesian_const_iterator_d qmaxit = construct_it(construct_max_vertex(q)),
	  qminit = construct_it(construct_min_vertex(q)), pit = construct_it(p);
      if (power == FT(0)) {
	  FT dist = 0;
          if ((*pit)>(*qmaxit)) dist = ((*pit)-(*qmaxit)); 
	  else if ((*pit)<(*qminit)) dist = ((*qminit)-(*pit));
          if (the_weights[0] * dist > distance)
	    distance = the_weights[0] * dist;
          qmaxit++;qminit++;pit++;
          if ((*pit)>(*qmaxit)) dist = ((*pit)-(*qmaxit)); 
	  else if ((*pit)<(*qminit)) dist = ((*qminit)-(*pit));
          if (the_weights[1] * dist > distance)
	    distance = the_weights[1] * dist;
      }
      else{
	  if ((*pit)>(*qmaxit)) distance += 
			the_weights[0] * std::pow((*pit)-(*qmaxit),power);
	  else if ((*pit)<(*qminit)) distance += 
			the_weights[0] * std::pow((*qminit)-(*pit),power);
          qmaxit++;qminit++;pit++;
          if ((*pit)>(*qmaxit)) distance += 
			the_weights[1] * std::pow((*pit)-(*qmaxit),power);
	  else if ((*pit)<(*qminit)) distance += 
			the_weights[1] * std::pow((*qminit)-(*pit),power);
      }
      return distance;
    }

    //DIM = 3 loop unrolled
    inline FT 
    transformed_distance(const Query_item& q, const Point_d& p, Dimension_tag<3> dt) const
    {
      FT distance = FT(0);
      typename SearchTraits::Construct_cartesian_const_iterator_d construct_it=
          traits.construct_cartesian_const_iterator_d_object();
      typename SearchTraits::Construct_min_vertex_d construct_min_vertex;
      typename SearchTraits::Construct_max_vertex_d construct_max_vertex;
      typename SearchTraits::Cartesian_const_iterator_d qmaxit = construct_it(construct_max_vertex(q)),
	  qminit = construct_it(construct_min_vertex(q)), pit = construct_it(p);
      if (power == FT(0)) {
	  FT dist = 0;
          if ((*pit)>(*qmaxit)) dist = ((*pit)-(*qmaxit)); 
	  else if ((*pit)<(*qminit)) dist = ((*qminit)-(*pit));
          if (the_weights[0] * dist > distance)
	    distance = the_weights[0] * dist;
          qmaxit++;qminit++;pit++;
          if ((*pit)>(*qmaxit)) dist = ((*pit)-(*qmaxit)); 
	  else if ((*pit)<(*qminit)) dist = ((*qminit)-(*pit));
          if (the_weights[1] * dist > distance)
	    distance = the_weights[1] * dist;
          qmaxit++;qminit++;pit++;
          if ((*pit)>(*qmaxit)) dist = ((*pit)-(*qmaxit)); 
	  else if ((*pit)<(*qminit)) dist = ((*qminit)-(*pit));
          if (the_weights[2] * dist > distance)
	    distance = the_weights[2] * dist;
      }
      else{
	   if ((*pit)>(*qmaxit)) distance += 
			the_weights[0] * std::pow((*pit)-(*qmaxit),power);
	  else if ((*pit)<(*qminit)) distance += 
			the_weights[0] * std::pow((*qminit)-(*pit),power);
          qmaxit++;qminit++;pit++;
          if ((*pit)>(*qmaxit)) distance += 
			the_weights[1] * std::pow((*pit)-(*qmaxit),power);
	  else if ((*pit)<(*qminit)) distance += 
			the_weights[1] * std::pow((*qminit)-(*pit),power);
          qmaxit++;qminit++;pit++;
          if ((*pit)>(*qmaxit)) distance += 
			the_weights[2] * std::pow((*pit)-(*qmaxit),power);
	  else if ((*pit)<(*qminit)) distance += 
			the_weights[2] * std::pow((*qminit)-(*pit),power);
      }
      return distance;
    }

    inline 
    FT 
    min_distance_to_rectangle(const Query_item& q,
			      const Kd_tree_rectangle<FT,Dimension>& r) const 
    {
      FT distance = FT(0);
		typename SearchTraits::Construct_cartesian_const_iterator_d construct_it=
                  traits.construct_cartesian_const_iterator_d_object();
		typename SearchTraits::Construct_min_vertex_d construct_min_vertex;
		typename SearchTraits::Construct_max_vertex_d construct_max_vertex;
                typename SearchTraits::Cartesian_const_iterator_d qmaxit = construct_it(construct_max_vertex(q)),
		  qe = construct_it(construct_max_vertex(q),1), qminit = construct_it(construct_min_vertex(q));
      if (power == FT(0))
	{
	  for (unsigned int i = 0; qmaxit != qe; ++ qmaxit, ++qminit, ++i) {
	    FT dist = 0;
            if (r.min_coord(i)>(*qmaxit)) 
	        dist =(r.min_coord(i)-(*qmaxit)); 
	    if (r.max_coord(i)<(*qminit)) 
		dist = ((*qminit)-r.max_coord(i));
            
            if ((the_weights[i]*dist) > distance)
	      distance = the_weights[i] * dist;
	  }
	}
      else
	{
	  for (unsigned int i = 0; qmaxit != qe; ++ qmaxit, ++qminit, ++i) {
	    if (r.min_coord(i)>(*qmaxit)) 
	        distance += the_weights[i] * std::pow(r.min_coord(i)-(*qmaxit),power);
	    if (r.max_coord(i)<(*qminit)) 
		distance += the_weights[i] * std::pow((*qminit)-r.max_coord(i),power);
           }
	}
      return distance;
    }

    inline 
    FT 
    min_distance_to_rectangle(const Query_item& q,
			      const Kd_tree_rectangle<FT,Dimension>& r,std::vector<FT>& dists) {
     FT distance = FT(0);
     unsigned int j;
		typename SearchTraits::Construct_cartesian_const_iterator_d construct_it=
                  traits.construct_cartesian_const_iterator_d_object();
		typename SearchTraits::Construct_min_vertex_d construct_min_vertex;
		typename SearchTraits::Construct_max_vertex_d construct_max_vertex;
                typename SearchTraits::Cartesian_const_iterator_d qmaxit = construct_it(construct_max_vertex(q)),
		  qe = construct_it(construct_max_vertex(q),1), qminit = construct_it(construct_min_vertex(q));
      if (power == FT(0))
	{
          FT temp;
	  for (unsigned int i = 0; qmaxit != qe; ++ qmaxit, ++qminit, ++i) {
	    FT dist = 0;
            if (r.min_coord(i)>(*qmaxit)) 
	        dist =(r.min_coord(i)-(*qmaxit)); 
	    if (r.max_coord(i)<(*qminit)) 
		dist = ((*qminit)-r.max_coord(i));
            
            if ((the_weights[i]*dist) > distance){
              j = i;
              temp = dist;
	      distance = the_weights[i] * dist;
            }
	  }
          dists[j] = temp;
	}
      else
	{
	  for (unsigned int i = 0; qmaxit != qe; ++ qmaxit, ++qminit, ++i) {
	    if (r.min_coord(i)>(*qmaxit)){
                dists[i] = r.min_coord(i)-(*qmaxit);
	        distance += the_weights[i] * std::pow(dists[i],power);
            }
	    if (r.max_coord(i)<(*qminit)) {
                dists[i] = (*qminit)-r.max_coord(i);
		distance += the_weights[i] * std::pow(dists[i],power);
            }
           }
	}
      return distance;
    }

    inline 
    FT
    max_distance_to_rectangle(const Query_item& q,
			      const Kd_tree_rectangle<FT,Dimension>& r) const {
      FT distance=FT(0);
     typename SearchTraits::Construct_cartesian_const_iterator_d construct_it=
                  traits.construct_cartesian_const_iterator_d_object();
		typename SearchTraits::Construct_min_vertex_d construct_min_vertex;
		typename SearchTraits::Construct_max_vertex_d construct_max_vertex;
                typename SearchTraits::Cartesian_const_iterator_d qmaxit = construct_it(construct_max_vertex(q)),
		  qe = construct_it(construct_max_vertex(q),1), qminit = construct_it(construct_min_vertex(q));
      if (power == FT(0))
	{
	  for (unsigned int i = 0; qmaxit != qe; ++ qmaxit, ++qminit, ++i) {
	    FT dist = 0;
            if (r.min_coord(i)>(*qmaxit)) 
	        dist =(r.max_coord(i)-(*qmaxit)); 
	    if (r.max_coord(i)<(*qminit)) 
		dist = ((*qminit)-r.min_coord(i));
            
	      if (the_weights[i] * dist > distance)
		distance = the_weights[i] * dist;
            }
	  }
      else
	{
	  for (unsigned int i = 0; qmaxit != qe; ++ qmaxit, ++qminit, ++i) {
            if (r.min_coord(i)>(*qmaxit)) 
	        distance += the_weights[i] * std::pow(r.max_coord(i)-(*qmaxit),power);
	    if (r.max_coord(i)<(*qminit)) 
		distance += the_weights[i] * std::pow((*qminit)-r.min_coord(i),power);
	  }
	}
      return distance;
    }

     inline 
    FT
    max_distance_to_rectangle(const Query_item& q,
			      const Kd_tree_rectangle<FT,Dimension>& r,std::vector<FT>& dists) {
      FT distance = FT(0);
     unsigned int j;
		typename SearchTraits::Construct_cartesian_const_iterator_d construct_it=
                  traits.construct_cartesian_const_iterator_d_object();
		typename SearchTraits::Construct_min_vertex_d construct_min_vertex;
		typename SearchTraits::Construct_max_vertex_d construct_max_vertex;
                typename SearchTraits::Cartesian_const_iterator_d qmaxit = construct_it(construct_max_vertex(q)),
		  qe = construct_it(construct_max_vertex(q),1), qminit = construct_it(construct_min_vertex(q));
      if (power == FT(0))
	{
          FT temp;
	  for (unsigned int i = 0; qmaxit != qe; ++ qmaxit, ++qminit, ++i) {
	    FT dist = 0;
            if (r.min_coord(i)>(*qmaxit)) 
	        dist =(r.max_coord(i)-(*qmaxit)); 
	    if (r.max_coord(i)<(*qminit)) 
		dist = ((*qminit)-r.min_coord(i));
            
	      if (the_weights[i] * dist > distance){
                j=i;
                temp = dist;
		distance = the_weights[i] * dist;
                }
            }
          dists[j] = temp;
	  }
      else
	{
	 for (unsigned int i = 0; qmaxit != qe; ++ qmaxit, ++qminit, ++i) {
            if (r.min_coord(i)>(*qmaxit)){
                dists[i] = r.max_coord(i)-(*qmaxit);
	        distance += the_weights[i] * std::pow(dists[i],power);
            }
	    if (r.max_coord(i)<(*qminit)) {
              dists[i] = (*qminit)-r.min_coord(i);
              distance += the_weights[i] * std::pow(dists[i],power);
            }
	  }
	}
      return distance;
    }
    
    inline 
    FT 
    new_distance(FT dist, FT old_off, FT new_off,
		 int cutting_dimension)  const 
    {
      FT new_dist;
      if (power == FT(0))
	{
	  if (the_weights[cutting_dimension]*CGAL::abs(new_off) 
	      > dist) 
	    new_dist= 
	      the_weights[cutting_dimension]*CGAL::abs(new_off);
	  else new_dist=dist;
	}
      else
	{
	  new_dist = dist + the_weights[cutting_dimension] * 
	    (std::pow(CGAL::abs(new_off),power)-std::pow(CGAL::abs(old_off),power));
	}
      return new_dist;
    }
    
    inline 
    FT 
    transformed_distance(FT d) const 
    {
      if (power <= FT(0)) return d;
      else return std::pow(d,power);
      
    }
    
    inline 
    FT 
    inverse_of_transformed_distance(FT d) const 
    {
      if (power <= FT(0)) return d;
      else return std::pow(d,1/power);
      
    }

  }; // Specialization for rectangles

} // namespace CGAL

#endif // CGAL_WEIGHTED_MINKOWSKI_DISTANCE_H
