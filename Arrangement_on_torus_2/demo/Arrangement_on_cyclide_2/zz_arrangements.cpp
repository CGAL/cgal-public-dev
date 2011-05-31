// ============================================================================
//
// Copyright (c) 2001-2010 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of EXACUS (http://www.mpi-inf.mpg.de/projects/EXACUS/);
// you may redistribute it under the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with EXACUS.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRcomANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// ----------------------------------------------------------------------------
//
// Library       : QdX
// File          : demos/xsurface/zz_arrangements.C
// QdX_release   : $Name:  $
// Revision      : $Revision$
// Revision_date : $Date$
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

#define NDEBUG 1
#include "include/arrangements.h"
#include "include/includes_ckva.h"

#include <CGAL/Algebraic_kernel_3/IO/Algebraic_surface_3_iostream.h>

#include <CGAL/Arrangement_on_torus_2.h>
#include <CGAL/Arrangement_2/Arrangement_on_surface_2_global.h>

#include <CGAL/Arr_overlay_2.h>
#include "include/Arr_red_blue_overlay_traits.h"
#include <CGAL/Arr_extended_dcel.h>

#include <CGAL/Algebraic_kernel_d/Polynomial_parser_d.h>

#include <CGAL/Timer.h>

typedef CGAL::Arr_extended_dcel< Geo_traits, 
    std::pair< bool, bool >, 
    std::pair< bool, bool >, 
    std::pair< bool, bool > > Extended_dcel;

typedef CGAL::Arrangement_on_torus_2< Geo_traits > Arr_on_surface_2;  

typedef CGAL::Arrangement_on_torus_2< Geo_traits, Extended_dcel > 
   Arr_on_surface_ext_2;  

typedef CGAL::Arr_red_blue_overlay_traits_base< Arr_on_surface_2,
      Arr_on_surface_2, Arr_on_surface_ext_2 > Arr_overlay_traits;

//! global symbols for arrangement computation

Arr_on_surface_2 arr_instance;

Points_3 arr_points; //! set of spatial arrangement event points and
XArcs_3 arr_xarcs;   //! x-monotone arcs

Base_surfaces base_surfaces; //! set of base surfaces
Surface_set surface_set[2];     //! set of surfaces intersecting the base surf

typedef std::map<std::size_t, int> Color_map;

bool XSurface_arrangements::base_surface_valid() const {
    return (base_index < base_surfaces.size());
}

// inserts to an empty arrangement
void arr_insert_empty(Arr_on_surface_2& arr, Geo_traits& geo, 
        const Surface_set& sset) {

    CGAL::set_pretty_mode(std::cout);
    XArcs_3 xcurves;
    Points_3 pts;

    CGAL::Timer mxm_timer, sweep_timer;
    mxm_timer.start();

    std::cout << "Doing make_x_monotone: ..\n\n";

    CGAL::make_x_monotone(sset.begin(), sset.end(),
        std::back_inserter(xcurves), std::back_inserter(pts), &geo);
    
    mxm_timer.stop();
    std::cout << "Time for make_x_monotone: " << mxm_timer.time() <<
        std::endl;    

    std::cout << "Doing sweep on surface: " << Geo_traits::instance().base()
        << "\n\n";
    
    arr.clear();
    arr = Arr_on_surface_2(&geo);
    
    sweep_timer.start();
    CGAL::insert_empty(arr, xcurves.begin(), xcurves.end(), pts.begin(),
            pts.end());
        
    sweep_timer.stop();
    std::cout << "Time for sweep on surface: " << sweep_timer.time() <<
        std::endl;    

    std::cout << "#V=" << arr.number_of_vertices() << ", ";
    std::cout << "#E=" << arr.number_of_edges() << ", ";
    std::cout << "#F=" << arr.number_of_faces() << std::endl;

}

void fill_containers(const Arr_on_surface_2& arr) {

    arr_points.clear();
    arr_xarcs.clear();
    arr_points.reserve(arr_instance.number_of_vertices());
    arr_xarcs.reserve(arr_instance.number_of_edges());

    Arr_on_surface_2::Vertex_const_iterator vit;
    for(vit = arr.vertices_begin(); vit != arr.vertices_end(); vit++) {
        Point_3 p = vit->point();
        arr_points.push_back(p);
    }

    Arr_on_surface_2::Edge_const_iterator eit;
    for(eit = arr.edges_begin(); eit != arr.edges_end(); eit++) 
        arr_xarcs.push_back(eit->curve()); 
} 

// picks up a color index based on id parameter
int pick_color(std::size_t id, Color_map& color_map) {

    int cindex = color_map.size();
    std::pair<Color_map::iterator, bool> ret =
        color_map.insert(std::make_pair(id, cindex));
    if(!ret.second) // already exists in the map
        cindex = ret.first->second;
    return cindex;
}

void XSurface_arrangements::compute_arr_on_surface() {

    if(base_index >= base_surfaces.size()) {
        std::cerr << "wrong base surface index\n";
        valid_flag = false;
        return;
    }

    Geo_traits geo(base_surfaces[base_index]);
    Geo_traits::set_instance(geo);
    arr_insert_empty(arr_instance, geo, surface_set[0]);

    fill_containers(arr_instance);
        
    xarcs_color_index.clear();
    xarcs_color_index.reserve(arr_xarcs.size());
    points_color_index.clear();
    
    Color_map cmap;
    for(XArcs_3::const_iterator xit = arr_xarcs.begin(); xit != 
            arr_xarcs.end(); xit++) 
        xarcs_color_index.push_back(pick_color(xit->curve().id(), cmap));
    
    valid_flag = true;
}

void XSurface_arrangements::compute_overlay_arr() {

    if(base_index >= base_surfaces.size()) {
        std::cerr << "wrong base surface index\n";
        valid_flag = false;
        return;
    }

    Geo_traits geo(base_surfaces[base_index]);
    Geo_traits::set_instance(geo);
    
    std::cout << "Computing 2nd arrangement..\n";
    Arr_on_surface_2 arr2;
    arr_insert_empty(arr2, geo, surface_set[1]);
    
    Arr_on_surface_ext_2 arr_overlay;
    Arr_overlay_traits ovl_traits;

    CGAL::Timer ovl_timer;
    std::cout << "Doing overlay arrangements..\n";
    ovl_timer.start();

    CGAL::overlay(arr_instance, arr2, arr_overlay, ovl_traits);
    
    ovl_timer.stop();
    std::cout << "Time for overlay: " << ovl_timer.time() << std::endl;   
    
    std::cout << "#V(O)=" << arr_overlay.number_of_vertices() << ", ";
    std::cout << "#E(O)=" << arr_overlay.number_of_edges() << ", ";
    std::cout << "#F(O)=" << arr_overlay.number_of_faces() << std::endl;
        
    //fill_containers(arr_overlay);
    // and assign colors
    // 0 - 1st arr, 1 - 2nd arr, 2 - both arrs
    unsigned color_map[] = {0, 3, 5};
    
    arr_points.clear();
    arr_points.reserve(arr_instance.number_of_vertices());

    points_approx.clear();
    points_approx.reserve(arr_instance.number_of_vertices());
    
    points_color_index.clear();
    points_color_index.reserve(arr_points.size());
    
    Arr_on_surface_ext_2::Vertex_const_iterator vit;
    for(vit = arr_overlay.vertices_begin(); 
        vit != arr_overlay.vertices_end(); vit++) {
        
        Point_3 p = vit->point();
        arr_points.push_back(p);
        unsigned idx = vit->data().first + (vit->data().second << 1) - 1;
        if(idx > 2) {
            std::cerr << "warning: incorrect overlay index: " << idx << "\n";
        }
        points_color_index.push_back(color_map[idx % 3]);
    }

    arr_xarcs.clear();
    arr_xarcs.reserve(arr_instance.number_of_edges());
    
    xarcs_color_index.clear();
    xarcs_color_index.reserve(arr_xarcs.size());
    
    Arr_on_surface_ext_2::Edge_const_iterator eit;
    for(eit = arr_overlay.edges_begin(); 
            eit != arr_overlay.edges_end(); eit++) {
        arr_xarcs.push_back(eit->curve());

        unsigned idx = eit->data().first + (eit->data().second << 1) - 1;
        if(idx > 2) {
            std::cerr << "warning: incorrect overlay index: " << idx << "\n";
        }
        xarcs_color_index.push_back(color_map[idx % 3]);
    }   
    valid_flag = true;
}

bool XSurface_arrangements::read_base_surfaces(const char *filename) {
    
    std::vector< Polynomial_3 > polys;
    base_surfaces.clear();
    valid_flag = false;

    typedef Base_surface_3::Matrix Matrix;
    typedef Base_surface_3::Vector Vector;
    
    std::ifstream ifstr(filename);
    if(!ifstr) {
        std::cerr << "Could not read cyclide bases file" << std::endl;
        return false;
    }
    std::string curr_tag;
    bool curr_parsing = false;
    Integer mu, base_radius_a, base_radius_b;

    Matrix basis(3);
    Vector center(3);
    
    while(ifstr >> curr_tag) {
        if (curr_tag == "BEGIN_CYCLIDE") {
            CGAL_assertion(!curr_parsing);
            base_radius_a = base_radius_b=2;
            mu = 1;
            
            // stupid way to define identity matrix
            const Base_surface_3::Coefficient help[] = {0, 0, 1, 0, 0};
            basis[0] = Vector(help+2, help+5);
            basis[1] = Vector(help+1, help+4);
            basis[2] = Vector(help, help+3);
            center = Vector(3);
            curr_parsing = true;
        }
        if(curr_tag == "END_CYCLIDE") {
            CGAL_assertion(curr_parsing);
            std::cout << "New cyclide with\n a = " << base_radius_a
              << "\n b = " << base_radius_b << "\n mu= " << mu
              << "\nbase plane =\n" << basis[0][0] << " " 
              << basis[0][1] << " " << basis[0][2] << "\n"
              << basis[1][0] << " " << basis[1][1] << " " 
              << basis[1][2] << "\n" << basis[2][0] << " " 
              << basis[2][1] << " " << basis[2][2] << "\n" << "center = " 
              << center[0] << " " << center[1] << " " << center[2]
                 << std::endl;
            base_surfaces.push_back(Base_surface_3(base_radius_a,
                 base_radius_b, mu, basis, center));
            curr_parsing = false;
        }
        if (curr_tag == "mu") {
            CGAL_assertion(curr_parsing);
            ifstr >> mu;
        }
        if (curr_tag == "a") {
            CGAL_assertion(curr_parsing);
            ifstr >> base_radius_a;
        }
        if (curr_tag == "b") {
            CGAL_assertion(curr_parsing);
            ifstr >> base_radius_b;
        }
        if (curr_tag == "base_plane") {
            CGAL_assertion(curr_parsing);
            for(int i = 0;i < 9; i++) {
                ifstr >> basis[i/3][i%3];
            }
        }
        if (curr_tag == "center") {
            CGAL_assertion(curr_parsing);
            for(int i = 0; i < 3; i++) {
                ifstr >> center[i];
            }
        }
    }

    std::cout << base_surfaces.size() << " base surfaces found in file.\n";  
    base_index = 0;  
    return (base_surfaces.size() > 0);
}

bool XSurface_arrangements::read_surface_set(const char *filename, 
        unsigned which, bool clear_flag) {

    CGAL_precondition(which < 2);
    if(clear_flag) {
        valid_flag = false;
        surface_set[which].clear();
    }
    
    typedef Arithmetic_kernel::Rational Rational;

    std::vector< Polynomial_3 > polys;
    if(!CGAL::read_file< Arithmetic_kernel >(filename,
              std::back_inserter(polys)) || polys.size() == 0) {
    
        std::cerr << "Trying other format" << std::endl;
       
        std::ifstream file(filename);
        if(!file) {
            std::cerr << "Could not read file" << std::endl;
            return false;
        }
        
        typedef CGAL::Polynomial_type_generator< Rational, 3 >::Type
          Poly_rat_3;

        CGAL::Polynomial_parser_d< 
              Poly_rat_3,
              CGAL::Mixed_rational_parser_policy< Poly_rat_3 > > 
          parser;
        
        while(!file.eof()) {    
          std::string str;
          file >> str;
          
          Poly_rat_3 fr;
          
          if(!parser(str, fr)) {
            std::cerr << "Could not parse .." << std::endl;
            return false;        
          }
          
          std::cout << fr << "\n\n";    
          
          Polynomial_3 f;
          typedef CGAL::Fraction_traits< Poly_rat_3 > FTraits;
          FTraits::Denominator_type det(1);
          FTraits::Decompose decompose;
          decompose(fr, f, det);
          
          polys.push_back(f);
        }
    }

//     if(polys.size() == 0) {
//         std::cerr << "Trying another format" << std::endl;
//         if(1/*!QdX::read_file< Arithmetic_kernel >(filename,
//              std::back_inserter(polys))*/) {
//             
//             std::cerr << "Could not read file" << std::endl;
//             return false;
//         }
//     }
    
    for(unsigned i = 0; i < polys.size(); i++) {
        if(polys[i] == base_surfaces[base_index].f()) {
            std::cerr << " duplicate surface found, skipping..\n";
            continue;
        }    
        surface_set[which].push_back(SURFACE_CONSTRUCTION(polys[i]));
    }
    std::cout << surface_set[which].size() << " total surfaces.\n";
    return true;
}


