// Copyright (c) 2010  Tel-Aviv University (Israel).
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
// Author(s)     : Asaf Porat          <asafpor1@post.tau.ac.il>

#ifndef LINE_THROUGH_POLYTOPES_CON_COMPONENT_H
#define LINE_THROUGH_POLYTOPES_CON_COMPONENT_H

#include <map>
#include <CGAL/Lines_through_segments_bounded_seg.h>

namespace CGAL {

template <typename Vertex, typename Edge, typename Traits_2_adapt>
class Lines_through_polytopes_con_component 
{
public:
  typedef enum {
    CGAL_LTP_COMP_UNINITIALIZED,
    CGAL_LTP_COMP_SINGLE,
    CGAL_LTP_COMP_LEFT,
    CGAL_LTP_COMP_RIGHT,
    CGAL_LTP_COMP_TOP,
    CGAL_LTP_COMP_BOTTOM,
    CGAL_LTP_COMP_TOP_LEFT,
    CGAL_LTP_COMP_TOP_RIGHT,
    CGAL_LTP_COMP_BOTTOM_LEFT,
    CGAL_LTP_COMP_BOTTOM_RIGHT
  } cgal_ltp_comp_type;

  class Lines_through_polytopes_component
  {
    /* Bitwise or holds the compenent type. */
    
    std::list<Edge> m_component;
    int m_index;
    Lines_through_polytopes_component* m_father;
    Traits_2_adapt m_traits_2_adapt;
    typedef typename Traits_2_adapt::Rational Rational;
    typedef typename Traits_2_adapt::Algebraic_real_1 Algebraic_real_1;
               
    cgal_ltp_comp_type m_comp_type;
    typedef typename Traits_2_adapt::Algebraic Algebraic;
    typedef Lines_through_segments_rbound_unbound_union<Algebraic> Rbound_alg_core;
    typedef Lines_through_segments_rbound_unbound_union<Algebraic_real_1> Rbound_alg_real;

    /* Holds the boundaries of the component. */
    Rbound_alg_real m_comp_max_x;
    Rbound_alg_real m_comp_min_x;
    Rbound_alg_core m_comp_max_y;
    Rbound_alg_core m_comp_min_y;
     
    
  public:
    typedef Lines_through_segments_rbound_unbound_union<Rational>  LTS_rbound;  
    typedef typename std::list<Edge>::iterator iterator;
    typedef typename std::list<Edge>::const_iterator const_iterator;
               
    /* iterator support */
    iterator begin() { return m_component.begin(); }
    const_iterator begin() const { return m_component.begin(); }
    iterator end() { return m_component.end(); }
    const_iterator end() const { return m_component.end(); }


    Lines_through_polytopes_component(int index)
    {
      m_index = index;
      m_father = NULL;
      m_comp_type = CGAL_LTP_COMP_UNINITIALIZED;
      m_comp_max_x = Rbound_alg_real(LTS_rbound::LTS_BS_UNBOUNDED_MINUS_INFINITY);
      m_comp_min_x = Rbound_alg_real(LTS_rbound::LTS_BS_UNBOUNDED_PLUS_INFINITY);
      m_comp_max_y = Rbound_alg_core(LTS_rbound::LTS_BS_UNBOUNDED_MINUS_INFINITY);
      m_comp_min_y = Rbound_alg_core(LTS_rbound::LTS_BS_UNBOUNDED_PLUS_INFINITY);
    }
               
    void add_element(const Edge& edge)
    {
      if (m_father == NULL)
        m_component.push_back(edge);
      else
        m_father->add_element(edge);
    }
               
    int size()
    {
      if (m_father == NULL)
        return m_component.size();
      else
        return m_father->size();
    }
               
    void splice(Lines_through_polytopes_component& comp)
    {
      Lines_through_polytopes_component* comp_ptr = &comp;
      while (comp_ptr->father() != NULL)
      {
        comp_ptr = comp_ptr->father();
      }
                  
      if (m_father == NULL)
      {
        /* Add the elements of comp to my list */
        m_component.splice(m_component.begin(),comp_ptr->m_component);
        comp_ptr->m_father = this;
        comp_ptr->m_index = this->index();
        CGAL_assertion(comp_ptr->m_component.size() == 0);
      }
      else
      {
        m_father->splice(*comp_ptr);
        comp_ptr->m_father = this->father();
        comp_ptr->m_index = this->index();
      }                  
    }

    int index()
    {
      if (m_father == NULL)
        return m_index;
      else
        return m_father->index();
    }

    Lines_through_polytopes_component* father()
    {
      return m_father;
    }

    void set_comp_type(cgal_ltp_comp_type type)
    {
      m_comp_type = type;
    }

    cgal_ltp_comp_type comp_type()
    {
      return m_comp_type;
    }

    Rbound_alg_core max_y()
    {
      return m_comp_max_y;
    }
    Rbound_alg_core min_y()
    {
      return m_comp_min_y;
    }
    Rbound_alg_real max_x()
    {
      return m_comp_max_x;
    }
    Rbound_alg_real min_x()
    {
      return m_comp_min_x;
    }

    /* The following function get over all of the vertical asymptote and
       adds a roof/floor above/below all of them.
    */

    void update_bounded_box()
    {
      for (iterator it = this->begin();
           it != this->end();
           ++it)
      {
        /* Update the bounded box of comp */
        set_min_max_x(m_comp_min_x, m_comp_max_x, *it);
        set_min_max_y(m_comp_min_y, m_comp_max_y, *it);
      }
    }
                  
    void eliminate_asymptotes(const Rational& global_max_y,
                              const Rational& global_min_y,
                              const Rational& global_max_x,
                              const Rational& global_min_x,
                              std::list<Edge>& vertical_edges_list,
                              std::list<Edge>& horizontal_edges_list)
    {

      bool first = true;
                  
      /* The asymptote type must be the same for all of the edges 
         at the connected component. */
      int asym_type = -1;
      /* The left vertex of the horizontal edge */
      Algebraic_real_1 left_vertex;

      /* The left vertex of the horizontal edge */
      Algebraic_real_1 right_vertex;
                  
      bool add_horizon_seg = false;


      /* Iterate over all of the edges and find the left most and
         right most edges. 
         For each edge cut the part of the edge that goes to infinity.
      */
      iterator it = this->begin();
      while (it != this->end())
      {
        Algebraic_real_1 new_x;
        /* If right or left vertices are at infinity bound them to
           the global min and max y. */
        Algebraic_real_1 right_x;
        Algebraic_real_1 left_x;
        bound_left_and_right_vertices(*it, right_x, left_x, global_max_x,
                                      global_min_x);
                                        
        bool update_bounderies = false;
                     
        get_left_and_right_bounds(                     
                                  global_max_y, global_min_y,
                                  global_max_x, global_min_x,
                                  first,
                                  right_x, left_x, update_bounderies,
                                  asym_type, new_x, it);
                     
        if (update_bounderies)
        {
          add_horizon_seg = true;
          set_left_and_right_vertices(
                                      left_vertex,
                                      right_vertex,
                                      new_x,
                                      first);
        }
      }
#if LTS_POLY_CON_COMP_DEBUG
      std::cout << "COMP MAX_Y = " << m_comp_max_y << std::endl;
      std::cout << "COMP MIN_Y = " << m_comp_min_y << std::endl;
      std::cout << "COMP MAX_X = " << m_comp_max_x << std::endl;
      std::cout << "COMP MIN_X = " << m_comp_min_x << std::endl;

      std::cout << "GLOBAL MAX_Y = " << global_max_y << std::endl;
      std::cout << "GLOBAL MIN_Y = " << global_min_y << std::endl;
      std::cout << "GLOBAL MAX_X = " << global_max_x << std::endl;
      std::cout << "GLOBAL MIN_X = " << global_min_x << std::endl;
#endif

                  
      if (add_horizon_seg)
      {
        /* One for each intersection of S_i (i = 1..2) with the polytope.
           - If size == 2, S_i intersect the polytope at two different edges
           and both of these lines are used, the two lines must be at the same
           component.

           - if size == 1, S_i intersect the polytope at only one point and
           the vertical line that was created is either an entry 
           point to the component or exit point or none.
           In each one of these three cases this line is used only in one of the 
           component.
        */
        for (iterator it = vertical_edges_list.begin();
             it != vertical_edges_list.end();
             ++it)
        {
          Algebraic_real_1 new_x_real = it->source_x();
                        
          /* Check if the vertical line as at the left or right components. */
          if (
              (m_comp_type != CGAL_LTP_COMP_RIGHT && 
               m_comp_type != CGAL_LTP_COMP_LEFT) ||
              (m_comp_min_x <= new_x_real &&
               m_comp_max_x >= new_x_real))
          {
            if (m_comp_min_x > new_x_real)
              m_comp_min_x = new_x_real;

            if (m_comp_max_x < new_x_real)
              m_comp_max_x = new_x_real;
                           
            set_left_and_right_vertices(left_vertex,
                                        right_vertex,
                                        it->source_x(),
                                        first);
          }
        }
      }
                  
      if (add_horizon_seg)
      {
        Edge horizontal_edge;
        set_bounded_horizontal_edge(horizontal_edge,
                                    global_min_y, global_max_y, left_vertex, right_vertex, asym_type);

#if LTS_POLY_CON_COMP_DEBUG
        std::cout << change_color(CGAL_BLUE,"NEW HORIZONTAL EDGE = ",
                                  horizontal_edge) << std::endl;
#endif                     
        m_component.push_back(horizontal_edge);
      }
                  
      for (iterator it = horizontal_edges_list.begin();
           it != horizontal_edges_list.end();
           ++it)
      {
        Algebraic y_coord_alg;
        m_traits_2_adapt.get_horizontal_asymptote_y_val(*it,y_coord_alg);

        Rbound_alg_core y_coord(y_coord_alg);
        /* One for each intersection of S_i (i = 1..2) with the polytope.
           - If size == 2, S_i intersect the polytope at two different edges
           and both of these lines are used.

           - if size == 1, S_i intersect the polytope at only one point and
           the horizontal line that was created is either an entry 
           point to the compenet or exit point or none.
           In each one of these three cases this line is used only in one of the 
           componenet.
        */

        if (
            (m_comp_type != CGAL_LTP_COMP_TOP && 
             m_comp_type != CGAL_LTP_COMP_BOTTOM) ||
            (m_comp_min_y <= y_coord &&
             m_comp_max_y >= y_coord))
        {
          Edge bounded_edge;
                        
          Algebraic_real_1 min_x_temp(global_min_x);
                        
          if (m_comp_min_x.is_bound())
          {
            Algebraic comp_min_x_core = 
              m_traits_2_adapt.convert_real_to_algebraic(m_comp_min_x.bound());


            if (global_min_x < comp_min_x_core)
              min_x_temp = m_comp_min_x.bound();
          }


          Algebraic_real_1 max_x_temp(global_max_x);
                        
          if (m_comp_max_x.is_bound())
          {
            Algebraic comp_max_x_core = 
              m_traits_2_adapt.convert_real_to_algebraic(m_comp_max_x.bound());
                           
            if (global_max_x > comp_max_x_core)
              max_x_temp = m_comp_max_x.bound();
          }
                        
          m_traits_2_adapt.create_curve_on_plane_arr(
                                                     bounded_edge,
                                                     min_x_temp,
                                                     max_x_temp,
                                                     *it);
                        
          m_component.push_back(bounded_edge);
        }
      }
    }
               
    void get_left_and_right_bounds(const Rational& global_max_y,
                                   const Rational& global_min_y,
                                   const Rational& global_max_x,
                                   const Rational& global_min_x,
                                   bool& first,
                                   const Algebraic_real_1& right_x,
                                   const Algebraic_real_1& left_x,
                                   bool& update_bounderies,
                                   int& asym_type,
                                   Algebraic_real_1& new_x,
                                   iterator& it)
    {
      Edge bounded_edge;

      if (m_traits_2_adapt.is_vertical(*it))
      {
        /* Remove the vertical segments.
           They do not affect the top and bottom envelope.*/
        new_x = it->source_x();
        iterator temp_it = it;
                     
        if (
            it->left_infinite_in_y() == ARR_TOP_BOUNDARY ||
            it->right_infinite_in_y() == ARR_TOP_BOUNDARY)
        {
          asym_type = ARR_TOP_BOUNDARY;
          update_bounderies = true;
        }
        else if (
                 it->left_infinite_in_y() == ARR_BOTTOM_BOUNDARY ||
                 it->right_infinite_in_y() == ARR_BOTTOM_BOUNDARY)
        {
          asym_type = ARR_BOTTOM_BOUNDARY;
          update_bounderies = true;
        }
                     
        it = m_component.erase(temp_it);
      }
      else if (it->left_infinite_in_y() == ARR_BOTTOM_BOUNDARY)
      {
        update_bounderies = true;
        iterator temp_it = it;
        Rational new_x_coord = m_traits_2_adapt.get_x_val(*it, global_min_y);
        new_x = Algebraic_real_1(new_x_coord);
        m_traits_2_adapt.create_curve_on_plane_arr(bounded_edge,
                                                   new_x_coord,
                                                   right_x,
                                                   *it);
                     
        it = m_component.erase(temp_it);
                     
        m_component.insert(it,bounded_edge);
        asym_type = ARR_BOTTOM_BOUNDARY;
      }
      else if (it->left_infinite_in_y() == ARR_TOP_BOUNDARY)
      {
        update_bounderies = true;
        iterator temp_it = it;
        Rational new_x_coord = m_traits_2_adapt.get_x_val(*it, global_max_y);
        new_x = Algebraic_real_1(new_x_coord);

        m_traits_2_adapt.create_curve_on_plane_arr(bounded_edge,
                                                   new_x_coord,
                                                   right_x,
                                                   *it);
                     
        it = m_component.erase(temp_it);
        m_component.insert(it,bounded_edge);
        asym_type = ARR_TOP_BOUNDARY;
      }
      else if (it->right_infinite_in_y() == ARR_BOTTOM_BOUNDARY)
      {
        update_bounderies = true;
        iterator temp_it = it;
        Rational new_x_coord = m_traits_2_adapt.get_x_val(*it, global_min_y);
        new_x = Algebraic_real_1(new_x_coord);
                     
        m_traits_2_adapt.create_curve_on_plane_arr(bounded_edge,
                                                   left_x,
                                                   new_x_coord,
                                                   *it);
                     
        it = m_component.erase(temp_it);
        m_component.insert(it,bounded_edge);
        asym_type = ARR_BOTTOM_BOUNDARY;
      }
      else if (it->right_infinite_in_y() == ARR_TOP_BOUNDARY)
      {
        update_bounderies = true;
        iterator temp_it = it;
        Rational new_x_coord = m_traits_2_adapt.get_x_val(*it, global_max_y);
        new_x = Algebraic_real_1(new_x_coord);
                     
        m_traits_2_adapt.create_curve_on_plane_arr(bounded_edge,
                                                   left_x,
                                                   new_x_coord,
                                                   *it);
                     
        it = m_component.erase(temp_it);
        m_component.insert(it,bounded_edge);
        asym_type = ARR_TOP_BOUNDARY;
      }
      else
      {
        iterator temp_it = it;
        m_traits_2_adapt.create_curve_on_plane_arr(bounded_edge,
                                                   left_x,
                                                   right_x,
                                                   *it);
                     
        it = m_component.erase(temp_it);
        m_component.insert(it,bounded_edge);
      }
    }
               

    ~Lines_through_polytopes_component()
    {
    }
               

    std::string to_string() const
    {
      std::ostringstream o;
                  
      for (const_iterator it = this->begin();
           it != this->end();
           it++)
      {
                     
        o <<  *it << std::endl;
      }
                  
      return o.str();
    }
               
  private:
               
    template <typename Algebraic_real_1>
    void bound_left_and_right_vertices(const Edge& edge,
                                       Algebraic_real_1& right_x,
                                       Algebraic_real_1& left_x,
                                       const Rational& global_max_x, 
                                       const Rational& global_min_x)
    {
      if (edge.right_infinite_in_x() == ARR_RIGHT_BOUNDARY)
        right_x = Algebraic_real_1(global_max_x);
      else
        right_x = edge.right_x();
                  
      if (edge.left_infinite_in_x() == ARR_LEFT_BOUNDARY)
        left_x = Algebraic_real_1(global_min_x);
      else
        left_x = edge.left_x();

    }
               
    template <typename NT>
    void set_min_max_y(NT& min_y, NT& max_y, const Edge& edge)
    {
      if (edge.left_infinite_in_y() == ARR_BOTTOM_BOUNDARY)
      {
        min_y = LTS_rbound::LTS_BS_UNBOUNDED_MINUS_INFINITY;
      }
      else if (edge.left_infinite_in_y() == ARR_TOP_BOUNDARY)
      {
        max_y = LTS_rbound::LTS_BS_UNBOUNDED_PLUS_INFINITY;
      }
      else/* ARR_INTERIOR */
      {
        if (m_traits_2_adapt.is_vertical(edge))
        {
          /* Get the bottom y val */
          Vertex left(edge.left());
          if (min_y > left.y())
            min_y = left.y();
          if (max_y < left.y())
            max_y = left.y();
        }
        else
        {
          if (edge.left_infinite_in_x() == ARR_INTERIOR)
          {
            NT y_val(m_traits_2_adapt.get_y_val(edge,edge.left_x()));
            if (min_y > y_val)
              min_y = y_val;
            if (max_y < y_val)
              max_y = y_val;
          }
          else
          {
            Algebraic y_coord_alg;
            m_traits_2_adapt.get_horizontal_asymptote_y_val(edge,y_coord_alg);
            NT temp(y_coord_alg);

            if (min_y > temp)
              min_y = temp;
            if (max_y < temp)
              max_y = temp;
          }
        }                 
      }
                  
      if (edge.right_infinite_in_y() == ARR_BOTTOM_BOUNDARY)
      {
        min_y = LTS_rbound::LTS_BS_UNBOUNDED_MINUS_INFINITY;
      }
      else if (edge.right_infinite_in_y() == ARR_TOP_BOUNDARY)
      {
        max_y = LTS_rbound::LTS_BS_UNBOUNDED_PLUS_INFINITY;
      }
      else/* ARR_INTERIOR */
      {
        if (m_traits_2_adapt.is_vertical(edge))
        {
          Vertex right(edge.right());
          if (min_y > right.y())
            min_y = right.y();
          if (max_y < right.y())
            max_y = right.y();
        }
        else
        {
          if (edge.right_infinite_in_x() == ARR_INTERIOR)
          {
            NT y_val(m_traits_2_adapt.get_y_val(edge,edge.right_x()));
            if (min_y > y_val)
              min_y = y_val;
            if (max_y < y_val)
              max_y = y_val;
          }
          else
          {
            Algebraic y_hor_asym;
            m_traits_2_adapt.get_horizontal_asymptote_y_val(edge,y_hor_asym);
            NT temp(y_hor_asym);
            if (min_y > temp)
              min_y = temp;
            if (max_y < temp)
              max_y = temp;
          }
        }      
      }
    }

    template <typename NT>
    void set_min_max_x(NT& min_x, NT& max_x, const Edge& edge)
    {
      if (edge.left_infinite_in_x() == ARR_LEFT_BOUNDARY)
      {
        min_x = LTS_rbound::LTS_BS_UNBOUNDED_MINUS_INFINITY;
      }
      else/* ARR_INTERIOR */
      {
        NT x_val = edge.left_x();
        if (min_x > x_val)
          min_x = x_val;
        if (max_x < x_val)
          max_x = x_val;
      }                 
                  

      if (edge.right_infinite_in_x() == ARR_RIGHT_BOUNDARY)
      {
        max_x = LTS_rbound::LTS_BS_UNBOUNDED_PLUS_INFINITY;
      }
      else/* ARR_INTERIOR */
      {
        NT x_val = edge.right_x();
        if (min_x > x_val)
          min_x = x_val;
        if (max_x < x_val)
          max_x = x_val;
      }
    }
               
    template <typename NT>
    void set_bounded_horizontal_edge(
                                     Edge& horizontal_edge,
                                     const Rational& min_y,
                                     const Rational& max_y,
                                     const NT& left_most_ver,
                                     const NT& right_most_ver,
                                     int asym_type)
    {

      switch(asym_type)
      {
       case ARR_BOTTOM_BOUNDARY:
        /* Create horizontal curve to bound from the bottom. */
        m_traits_2_adapt.create_horizontal_curve_on_plane_arr(
                                                              horizontal_edge,
                                                              min_y,
                                                              left_most_ver,
                                                              right_most_ver);
        break;
                        
       case ARR_TOP_BOUNDARY:
        /* Create horizontal curve to bound from the bottom. */
        m_traits_2_adapt.create_horizontal_curve_on_plane_arr(
                                                              horizontal_edge,
                                                              max_y,
                                                              left_most_ver,
                                                              right_most_ver);
        break;

       default:
        CGAL_error_msg("Unexpected boundary type");
        break;

      }
    }
               
    template <typename NT>
    void set_left_and_right_vertices(
                                     NT& left_most_ver,
                                     NT& right_most_ver,
                                     const NT& new_x,
                                     bool& first)
    {
      if (first)
      {
        first = false;
        left_most_ver = new_x;
        right_most_ver = new_x;
      }
      else
      {
        if (new_x < left_most_ver)
        {
          left_most_ver = new_x;
        }
        if (new_x > right_most_ver)
        {
          right_most_ver = new_x;
        }
      }
                  
    }
               
    bool set_asym_type(bool first,int& asym_type,int expected)
    {
      if (first)
      {
        asym_type = expected;
        return true;
      }
      else
      {
        return (asym_type == expected);
      }
    }
               
  };

private:            
  class Lines_through_polytopes_map
  {
    int m_comp_index;
    int m_num_of_comp;
    std::map<const Vertex, Lines_through_polytopes_component*> m_map;
    std::vector<Lines_through_polytopes_component*> m_con_comps;
               

  public:
    typedef typename std::map<const Vertex, Lines_through_polytopes_component*>::iterator 
    map_iterator;

    typedef typename std::map<const Vertex, Lines_through_polytopes_component*>::const_iterator 
    map_const_iterator;

    typedef typename std::vector<Lines_through_polytopes_component*>::iterator 
    comp_iterator;

    typedef typename std::vector<Lines_through_polytopes_component*>::const_iterator 
    comp_const_iterator;
               
    /* iterator support */
    map_iterator map_begin() { return m_map.begin(); }
    map_iterator map_end() { return m_map.end(); }

    map_const_iterator map_begin() const { return m_map.begin(); }
    map_const_iterator map_end() const { return m_map.end(); }

    comp_iterator comp_begin() { return m_con_comps.begin(); }
    comp_iterator comp_end() { return m_con_comps.end(); }

    comp_const_iterator comp_begin() const { return m_con_comps.begin(); }
    comp_const_iterator comp_end() const { return m_con_comps.end(); }

    Lines_through_polytopes_map()
    {
      m_comp_index = 0;
      m_num_of_comp = 0;



    }
               
    ~Lines_through_polytopes_map()
    {
      comp_iterator it;
      for (it = comp_begin();
           it != comp_end();
           it++)
      {
        delete(*it);
      }
    }
               
    void add_edge(const Edge& edge, 
                  Vertex& left,
                  Vertex& right)
    {
      if (edge.left_infinite_in_x() == CGAL::ARR_INTERIOR &&
          edge.left_infinite_in_y() == CGAL::ARR_INTERIOR)
      {
        left = edge.left();
      }
      else if (edge.left_infinite_in_x() == CGAL::ARR_INTERIOR)
      {
        /* y at infinity */
        if (edge.left_infinite_in_y() == ARR_BOTTOM_BOUNDARY)
           left = Vertex(edge.left_x(),Vertex::LTS_VERTEX_Y_MINUS_INFINITY);
        else
          left = Vertex(edge.left_x(),Vertex::LTS_VERTEX_Y_PLUS_INFINITY);
      }
      else
      {
        Traits_2_adapt traits_2_adapt;
        /* x at infinity */
        typename Traits_2_adapt::Algebraic y_coord_alg;
        traits_2_adapt.get_horizontal_asymptote_y_val(edge,y_coord_alg);
        left = Vertex(y_coord_alg, Vertex::LTS_VERTEX_X_MINUS_INFINITY);
      }
                  
      if (
          edge.right_infinite_in_x() == CGAL::ARR_INTERIOR &&
          edge.right_infinite_in_y() == CGAL::ARR_INTERIOR)
      {
        right = edge.right();
      }
      else if (edge.right_infinite_in_x() == CGAL::ARR_INTERIOR)
      {
        /* y at infinity */
        if (edge.right_infinite_in_y() == ARR_BOTTOM_BOUNDARY)
          right = Vertex(edge.right_x(),Vertex::LTS_VERTEX_Y_MINUS_INFINITY);
        else
          right = Vertex(edge.right_x(),Vertex::LTS_VERTEX_Y_PLUS_INFINITY);
      }
      else
      {
        /* x at infinity */
        Traits_2_adapt traits_2_adapt;
        typename Traits_2_adapt::Algebraic asym_y_val;
        traits_2_adapt.get_horizontal_asymptote_y_val(edge,asym_y_val);
        right = Vertex(asym_y_val,Vertex::LTS_VERTEX_X_PLUS_INFINITY);
      }

      map_iterator ver_map_it_left = m_map.find(left);
      map_iterator ver_map_it_right = m_map.find(right);
                  
      if (ver_map_it_left != m_map.end())
      {
        /* The left element is alreay at the map. */
        ver_map_it_left->second->add_element(edge);
                     
        if (ver_map_it_right == m_map.end())
          m_map.insert(std::pair<const Vertex, Lines_through_polytopes_component*>(
                                                                                   right, ver_map_it_left->second));
        else
        {
          /* Connect two connected components. */
          if (ver_map_it_left->second->index() != ver_map_it_right->second->index())
          {
            ver_map_it_right->second->splice(*ver_map_it_left->second);
            m_num_of_comp--;
          }
        }
      }
      else if (ver_map_it_right != m_map.end())
      {
        /* right element is alreay at the map. */
        ver_map_it_right->second->add_element(edge);
        m_map.insert( std::pair<const Vertex,Lines_through_polytopes_component*>(
                                                                                 left, ver_map_it_right->second));
      }
      else
      {
        Lines_through_polytopes_component* new_comp = 
          new Lines_through_polytopes_component(m_comp_index);
        m_con_comps.push_back(new_comp);
        new_comp->add_element(edge);
                     
        m_map.insert( std::pair<const Vertex,Lines_through_polytopes_component*>(
                                                                                 right, new_comp));
        m_map.insert( std::pair<const Vertex,Lines_through_polytopes_component*>(
                                                                                 left, new_comp));

        m_comp_index++;
        m_num_of_comp++;
      }
    }

    comp_iterator erase(comp_iterator to_erase)
    {
      return m_con_comps.erase(to_erase);
    }
               
    int num_of_comp() const
    {
      return m_num_of_comp;
    }
               
    void dec_num_of_comp()
    {
      m_num_of_comp--;
    }

    map_iterator find(const Vertex& ver)
    {
      return m_map.find(ver);
    }
               
    int size() const
    {
      return m_map.size();
    }
               
  };
         
  /* Member variables. */
  Lines_through_polytopes_map m_vertices_map;
  Traits_2_adapt m_traits_2_adapt;
  /* These edges were created from intersection of S1 with one of the segments. */
  std::list<Edge> m_vertical_unbounded_edges_list;
  /* These edges were created from intersection of S2 with one of the segments. */
  std::list<Edge> m_horizontal_unbounded_edges_list;
  /* Since the envelope does not support y = infinite, we 
     bound the vertical asymptote. at the maximum y that we
     use.
  */
  typename Traits_2_adapt::Algebraic m_max_y;
  typename Traits_2_adapt::Algebraic m_min_y;
  static const int LTS_CON_COMP_MAX = 10;
  typename Traits_2_adapt::Rational m_max_y_rat;
  typename Traits_2_adapt::Rational m_min_y_rat;

  typename Traits_2_adapt::Algebraic m_max_x;
  typename Traits_2_adapt::Algebraic m_min_x;
  typename Traits_2_adapt::Rational m_max_x_rat;
  typename Traits_2_adapt::Rational m_min_x_rat;
  
public:      
         
  typedef Lines_through_polytopes_component Lines_through_polytopes_component;
  typedef typename Lines_through_polytopes_component::LTS_rbound LTS_rbound;
  typedef typename Lines_through_polytopes_map::comp_iterator iterator;

  typedef typename std::list<Edge>::iterator vertical_edges_iterator;
  typedef typename std::list<Edge>::iterator horizontal_edges_iterator;
         
  /* iterator support */
  iterator begin() { return m_vertices_map.comp_begin(); }
  iterator end() { return m_vertices_map.comp_end(); }

  vertical_edges_iterator vertical_edges_begin() 
  {
    return m_vertical_unbounded_edges_list.begin();
  }
  vertical_edges_iterator vertical_edges_end() 
  { 
    return m_vertical_unbounded_edges_list.end();
  }

  horizontal_edges_iterator horizon_edges_begin() 
  {
    return m_horizontal_unbounded_edges_list.begin();
  }
  horizontal_edges_iterator horizon_edges_end() 
  { 
    return m_horizontal_unbounded_edges_list.end(); 
  }

  template <typename Input_edge_iterator>
  Lines_through_polytopes_con_component(Input_edge_iterator edges_begin,
                                        Input_edge_iterator edges_end)
  {
    Vertex left,right;
    typename Lines_through_polytopes_map::map_iterator ver_map_it_left, ver_map_it_right;
            
    Input_edge_iterator edges_it = edges_begin;

    m_max_y = LTS_CON_COMP_MAX;
    m_min_y = -LTS_CON_COMP_MAX;
    m_max_y_rat = 0;
    m_min_y_rat = 0;

    m_max_x = LTS_CON_COMP_MAX;
    m_min_x = -LTS_CON_COMP_MAX;
    m_max_x_rat = 0;
    m_min_x_rat = 0;

    /* Find all of the vertices and add them to m_vertices_map. */
    m_vertices_map = Lines_through_polytopes_map();
    while (edges_it != edges_end)
    {
      /* Ignore vertical and horizontal edges that are
         from minus infinity to plus infinity. These edges were 
         created from intersection of edge with S1 or S2.
         Later on these edges will be added to all of the connected components. */

      if (m_traits_2_adapt.is_vertical(*edges_it))
      {
        /* If unbounded */
        if (edges_it->left_infinite_in_y() != CGAL::ARR_INTERIOR &&
            edges_it->right_infinite_in_y() != CGAL::ARR_INTERIOR)
        {
          m_vertical_unbounded_edges_list.push_back(*edges_it);  
        }
        else
        {
          m_vertices_map.add_edge(*edges_it,left,right);
          update_min_max_x_y(left,right);
        }                 
      }
      else if (m_traits_2_adapt.is_horizontal(*edges_it) && 
               edges_it->right_infinite_in_x() != CGAL::ARR_INTERIOR &&
               edges_it->left_infinite_in_x() != CGAL::ARR_INTERIOR)
      {
        m_horizontal_unbounded_edges_list.push_back(*edges_it);
      }
      else
      {
        m_vertices_map.add_edge(*edges_it,left,right);
        update_min_max_x_y(left,right);
      }
      edges_it++;
    }

    /* Update the bounding box of all componenets and remove empty components.*/
    typename Lines_through_polytopes_map::comp_iterator comp_it;
    comp_it = m_vertices_map.comp_begin();
    while (comp_it != m_vertices_map.comp_end())
    {
      if ((*comp_it)->size() == 0)
      {
        comp_it = m_vertices_map.erase(comp_it);
      }
      else
      {
        (*comp_it)->update_bounded_box();
        comp_it++;
      }
    }
    set_componenets_type();

    /* Connect componenets that are connected with infinite vertical line.
       These components are connected at the start/end point of S1 and
       the edges are degenerate hyperoblas.
       This scenario takes place when the imagnary plane that sweep through
       S2 meet the polytope at one of the end points of S1.
    */
    comp_it = m_vertices_map.comp_begin();
            
    for (comp_it = m_vertices_map.comp_begin();
         comp_it != m_vertices_map.comp_end();
         comp_it++)
    {
      typename Lines_through_polytopes_map::comp_iterator comp_it2;
      comp_it2 = m_vertices_map.comp_begin();
      while(comp_it2 != m_vertices_map.comp_end())
      {
        if (
            (**comp_it2).index() != (**comp_it).index() &&
            common_comps(comp_it2,comp_it) && 
            !opposite_comps(comp_it2,comp_it) &&
            (contain_the_same_vertical_line(comp_it2,comp_it) ||
             contain_the_same_horizontal_line(comp_it2,comp_it)))
        {
          (*comp_it)->splice(**comp_it2);
          comp_it2 = m_vertices_map.erase(comp_it2);
          m_vertices_map.dec_num_of_comp();
        }
        else
        {
          comp_it2++;
        }
      }
    }

    /* Connect components that are connected with a vertex at another component.
       This scenario occurs when edges are split into two components by vertical
       asymptote.
       These components are of the same type.
    */
    comp_it = m_vertices_map.comp_begin();
            
    for (comp_it = m_vertices_map.comp_begin();
         comp_it != m_vertices_map.comp_end();
         comp_it++)
    {
      typename Lines_through_polytopes_map::comp_iterator comp_it2;
      comp_it2 = m_vertices_map.comp_begin();
      while(comp_it2 != m_vertices_map.comp_end())
      {
        if (
            (**comp_it2).index() != (**comp_it).index() &&
            (**comp_it2).comp_type() == (**comp_it).comp_type())
        {
          (*comp_it)->splice(**comp_it2);
          comp_it2 = m_vertices_map.erase(comp_it2);
          m_vertices_map.dec_num_of_comp();
        }
        else
        {
          comp_it2++;
        }
      }
    }

            
    /* Update the bounding box of all componenets */
    for (comp_it = m_vertices_map.comp_begin();
         comp_it != m_vertices_map.comp_end();
         comp_it++)
    {
      (*comp_it)->update_bounded_box();
    }
            
    CGAL_assertion(m_vertices_map.num_of_comp() <= 4);
  }

  bool contain_the_same_vertical_line(typename Lines_through_polytopes_map::comp_iterator comp_it1,
                                      typename Lines_through_polytopes_map::comp_iterator comp_it2)
  {
    for (vertical_edges_iterator vei = this->vertical_edges_begin();
         vei != this->vertical_edges_end();
         vei++)
    {
      if (
          (**comp_it1).max_x() >= vei->source_x() &&
          (**comp_it1).min_x() <= vei->source_x() &&
          (**comp_it2).max_x() >= vei->source_x() &&
          (**comp_it2).min_x() <= vei->source_x())
        return true;
    }

    return false;
  }

  bool contain_the_same_horizontal_line(typename Lines_through_polytopes_map::comp_iterator comp_it1,
                                        typename Lines_through_polytopes_map::comp_iterator comp_it2)
  {
    for (horizontal_edges_iterator hei = this->horizon_edges_begin();
         hei != this->horizon_edges_end();
         hei++)
    {
      typename Traits_2_adapt::Algebraic y_val;
      m_traits_2_adapt.get_horizontal_asymptote_y_val(*hei,y_val);

      if (
          (**comp_it1).max_y() >= y_val &&
          (**comp_it1).min_y() <= y_val &&
          (**comp_it2).max_y() >= y_val &&
          (**comp_it2).min_y() <= y_val)
        return true;
    }
    return false;
  }
         
  bool common_comps(
                    typename Lines_through_polytopes_map::comp_iterator comp_it1,
                    typename Lines_through_polytopes_map::comp_iterator comp_it2)
  {
    if (
        (**comp_it1).max_y() == LTS_rbound::LTS_BS_UNBOUNDED_PLUS_INFINITY &&
        (**comp_it2).max_y() == LTS_rbound::LTS_BS_UNBOUNDED_PLUS_INFINITY)
      return true;

    if (
        (**comp_it2).min_y() == LTS_rbound::LTS_BS_UNBOUNDED_MINUS_INFINITY &&
        (**comp_it1).min_y() == LTS_rbound::LTS_BS_UNBOUNDED_MINUS_INFINITY)
      return true;

    if (
        (**comp_it2).max_x() == LTS_rbound::LTS_BS_UNBOUNDED_PLUS_INFINITY &&
        (**comp_it1).max_x() == LTS_rbound::LTS_BS_UNBOUNDED_PLUS_INFINITY)
      return true;

    if (
        (**comp_it1).min_x() == LTS_rbound::LTS_BS_UNBOUNDED_MINUS_INFINITY &&
        (**comp_it2).min_x() == LTS_rbound::LTS_BS_UNBOUNDED_MINUS_INFINITY)
      return true;
            
    return false;
  }

  bool opposite_comps(typename Lines_through_polytopes_map::comp_iterator comp_it1,
                      typename Lines_through_polytopes_map::comp_iterator comp_it2)
  {
    if ((**comp_it1).max_y() == LTS_rbound::LTS_BS_UNBOUNDED_PLUS_INFINITY &&
        (**comp_it2).min_y() == LTS_rbound::LTS_BS_UNBOUNDED_MINUS_INFINITY)
      return true;

    if ((**comp_it2).max_y() == LTS_rbound::LTS_BS_UNBOUNDED_PLUS_INFINITY &&
        (**comp_it1).min_y() == LTS_rbound::LTS_BS_UNBOUNDED_MINUS_INFINITY)
      return true;

    if ((**comp_it2).max_x() == LTS_rbound::LTS_BS_UNBOUNDED_PLUS_INFINITY &&
        (**comp_it1).min_x() == LTS_rbound::LTS_BS_UNBOUNDED_MINUS_INFINITY)
      return true;

    if ((**comp_it1).max_x() == LTS_rbound::LTS_BS_UNBOUNDED_PLUS_INFINITY &&
        (**comp_it2).min_x() == LTS_rbound::LTS_BS_UNBOUNDED_MINUS_INFINITY)
      return true;
            
    return false;
  }
         
  void set_componenets_type()
  {
    /* Set the type of componenet */
    typename Lines_through_polytopes_map::comp_iterator it = m_vertices_map.comp_begin();
            
    for (it = m_vertices_map.comp_begin();
         it != m_vertices_map.comp_end();
         it++)
    {
      if (
          (*it)->max_y() == LTS_rbound::LTS_BS_UNBOUNDED_PLUS_INFINITY &&
          (*it)->max_x() == LTS_rbound::LTS_BS_UNBOUNDED_PLUS_INFINITY)
      {
        (*it)->set_comp_type(CGAL_LTP_COMP_TOP_RIGHT);
      }
      else if (
               (*it)->max_y() == LTS_rbound::LTS_BS_UNBOUNDED_PLUS_INFINITY &&
               (*it)->min_x() == LTS_rbound::LTS_BS_UNBOUNDED_MINUS_INFINITY)
      {
        (*it)->set_comp_type(CGAL_LTP_COMP_TOP_LEFT);
      }
      else if (
               (*it)->min_y() == LTS_rbound::LTS_BS_UNBOUNDED_MINUS_INFINITY &&
               (*it)->max_x() == LTS_rbound::LTS_BS_UNBOUNDED_PLUS_INFINITY)
      {
        (*it)->set_comp_type(CGAL_LTP_COMP_BOTTOM_RIGHT);
      }
      else if (
               (*it)->min_y() == LTS_rbound::LTS_BS_UNBOUNDED_MINUS_INFINITY &&
               (*it)->min_x() == LTS_rbound::LTS_BS_UNBOUNDED_MINUS_INFINITY)
      {
        (*it)->set_comp_type(CGAL_LTP_COMP_BOTTOM_LEFT);
      }
      else if ((*it)->max_y() == LTS_rbound::LTS_BS_UNBOUNDED_PLUS_INFINITY)
      {
        (*it)->set_comp_type(CGAL_LTP_COMP_TOP);
      }
      else if ((*it)->min_y() == LTS_rbound::LTS_BS_UNBOUNDED_MINUS_INFINITY)
      {
        (*it)->set_comp_type(CGAL_LTP_COMP_BOTTOM);
      }
      else if ((*it)->max_x() == LTS_rbound::LTS_BS_UNBOUNDED_PLUS_INFINITY)
      {
        (*it)->set_comp_type(CGAL_LTP_COMP_RIGHT);
      }
      else if ((*it)->min_x() == LTS_rbound::LTS_BS_UNBOUNDED_MINUS_INFINITY)
      {
        (*it)->set_comp_type(CGAL_LTP_COMP_LEFT);
      }
      else
      {
        (*it)->set_comp_type(CGAL_LTP_COMP_SINGLE);
      }
    }
  }
         
  void update_min_max_x_y(const Vertex& left,
                          const Vertex& right)
  {
    if (left.type() != Vertex::LTS_VERTEX_Y_PLUS_INFINITY &&
        left.type() != Vertex::LTS_VERTEX_Y_MINUS_INFINITY &&
        left.y() > m_max_y)
    {
      m_max_y = left.y();
    }

    if (left.type() != Vertex::LTS_VERTEX_Y_PLUS_INFINITY &&
        left.type() != Vertex::LTS_VERTEX_Y_MINUS_INFINITY &&
        left.y() < m_min_y)
    {
      m_min_y = left.y();
    }

    if (right.type() != Vertex::LTS_VERTEX_Y_PLUS_INFINITY &&
        right.type() != Vertex::LTS_VERTEX_Y_MINUS_INFINITY &&
        right.y() > m_max_y)
    {
      m_max_y = right.y();
    }

    if (right.type() != Vertex::LTS_VERTEX_Y_PLUS_INFINITY &&
        right.type() != Vertex::LTS_VERTEX_Y_MINUS_INFINITY &&
        right.y() < m_min_y)
    {
      m_min_y = right.y();
    }               
    /* Update min max x */
    if (left.type() != Vertex::LTS_VERTEX_X_PLUS_INFINITY &&
        left.type() != Vertex::LTS_VERTEX_X_MINUS_INFINITY &&
        left.x() > m_max_x)
    {
      m_max_x = left.x();
    }

    if (left.type() != Vertex::LTS_VERTEX_X_PLUS_INFINITY &&
        left.type() != Vertex::LTS_VERTEX_X_MINUS_INFINITY &&
        left.x() < m_min_x)
    {
      m_min_x = left.x();
    }

    if (right.type() != Vertex::LTS_VERTEX_X_PLUS_INFINITY &&
        right.type() != Vertex::LTS_VERTEX_X_MINUS_INFINITY &&
        right.x() > m_max_x)
    {
      m_max_x = right.x();
    }

    if (right.type() != Vertex::LTS_VERTEX_X_PLUS_INFINITY &&
        right.type() != Vertex::LTS_VERTEX_X_MINUS_INFINITY &&
        right.x() < m_min_x)
    {
      m_min_x = right.x();
    }               
  }


  void eliminate_asymptotes(Lines_through_polytopes_component& comp)
  {
    comp.eliminate_asymptotes(this->max_y(),this->min_y(),
                              this->max_x(),this->min_x(),
                              this->m_vertical_unbounded_edges_list,
                              this->m_horizontal_unbounded_edges_list);
            
  }
         
  typename Traits_2_adapt::Rational max_y()
  {
    if (m_max_y_rat != 0)
      return m_max_y_rat;
            
    typename Traits_2_adapt::Algebraic temp(LTS_CON_COMP_MAX);
    m_max_y_rat = LTS_CON_COMP_MAX;
    while (temp < m_max_y)
    {
      temp *= 2;
      m_max_y_rat *= 2;
    }

    return m_max_y_rat;
  }

  typename Traits_2_adapt::Rational min_y()
  {
    if (m_min_y_rat != 0)
      return m_min_y_rat;
            
    typename Traits_2_adapt::Algebraic temp(-LTS_CON_COMP_MAX);
    m_min_y_rat = -LTS_CON_COMP_MAX;
    while (temp > m_min_y)
    {
      temp *= 2;
      m_min_y_rat *= 2;
    }

    return m_min_y_rat;
  }

  typename Traits_2_adapt::Rational max_x()
  {
    if (m_max_x_rat != 0)
      return m_max_x_rat;
            
    typename Traits_2_adapt::Algebraic temp(LTS_CON_COMP_MAX);
    m_max_x_rat = LTS_CON_COMP_MAX;
    while (temp < m_max_x)
    {
      temp *= 2;
      m_max_x_rat *= 2;
    }

    return m_max_x_rat;
  }

  typename Traits_2_adapt::Rational min_x()
  {
    if (m_min_x_rat != 0)
      return m_min_x_rat;
            
    typename Traits_2_adapt::Algebraic temp(-LTS_CON_COMP_MAX);
    m_min_x_rat = -LTS_CON_COMP_MAX;
    while (temp > m_min_x)
    {
      temp *= 2;
      m_min_x_rat *= 2;
    }

    return m_min_x_rat;
  }
         
public:
  std::string to_string() const
  {
    std::ostringstream o;
    o << "num of vertices = " << m_vertices_map.size() << std::endl;
            
    for (typename Lines_through_polytopes_map::map_const_iterator it = m_vertices_map.map_begin();
         it != m_vertices_map.map_end();
         it++)
    {
      o << (it->first) << " Index = " << it->second->index() << std::endl;
    }
            
    o << "m_num_of_comp = " << m_vertices_map.num_of_comp() << std::endl;
            
    for (typename Lines_through_polytopes_map::comp_const_iterator it = m_vertices_map.comp_begin();
         it != m_vertices_map.comp_end();
         it++)
    {
      o << "New comp index = " << (*it)->index() << std::endl;
               
      o << (*it)->to_string() << std::endl;
    }
            
    return o.str();
  }
};

template <typename Vertex, typename Edge, typename Traits_2_adapt>
inline std::ostream&
operator<<(std::ostream & out,
           const Lines_through_polytopes_con_component<Vertex, Edge,
                                                       Traits_2_adapt>& to_print)
{
  out << to_print.to_string();
  return out;
}

   
} //namespace CGAL


#endif /*LINE_THROUGH_POLYTOPES_CON_COMPONENT_H*/
