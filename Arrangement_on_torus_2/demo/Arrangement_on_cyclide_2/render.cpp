// ============================================================================
//
// Copyright (c) 2001-2007 Max-Planck-Institut Saarbruecken (Germany).
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
// File          : demos/xsurface/render.C
// QdX_release   : $Name:  $
// Revision      : $Revision$
// Revision_date : $Date$
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

//!@file render.C
//! collects everything which relates to approximating spatial points and arcs

#define NDEBUG 1
#define CGAL_CKVA_COMPILE_RENDERER 
#define CGAL_CKVA_RENDER_WITH_REFINEMENT

#ifndef CGAL_AK_ENABLE_DEPRECATED_INTERFACE
#define CGAL_AK_ENABLE_DEPRECATED_INTERFACE 1
#endif

#include "include/arrangements.h"
#include "include/includes_ckva.h"

#include <boost/format.hpp>

#include <CGAL/Timer.h>

extern Points_3 arr_points; 
extern XArcs_3 arr_xarcs;  

extern Base_surfaces base_surfaces;
extern Surface_set surface_set;    

Point_3::Poly_double_2 cv[4];

static void compute_equation() {

    typedef Point_3::Poly_double_2 Poly_2;
    //CGAL::Polynomial_traits_d< Poly_2 >::Derivative der;
    
    const Poly_2* params = Point_3::param_surface();

//     Poly_2 pwsqr = pw * pw, pwq = pwsqr * pwsqr, pw6 = pwq * pwsqr,
//         dws = der(pw, 0), dwt = der(pw, 1);
//     
//     Poly_2 pxds = der(px, 0) * pw - dws * px,
//            pyds = der(py, 0) * pw - dws * py,
//            pzds = der(pz, 0) * pw - dws * pz;
//     
//     Poly_2 pxdt = der(px, 1) * pw - dwt * px,
//            pydt = der(py, 1) * pw - dwt * py,
//            pzdt = der(pz, 1) * pw - dwt * pz;
// 
//     Poly_2 dwss = der(dws, 0), dwst = der(dws, 1), dwtt = der(dwt, 1);
// 
//     Poly_2 pxdss = der(pxds, 0) * dws - dwss * pxds,
//            pydss = der(pyds, 0) * dws - dwss * pyds,
//            pzdss = der(pzds, 0) * dws - dwss * pzds;
//     
//     Poly_2 pxdst = der(pxds, 1) * dws - dwst * pxds,
//            pydst = der(pyds, 1) * dws - dwst * pyds, 
//            pzdst = der(pzds, 1) * dws - dwst * pzds;
// 
//     Poly_2 pxdtt = der(pxdt, 1) * dwt - dwtt * pxdt,
//            pydtt = der(pydt, 1) * dwt - dwtt * pydt, 
//            pzdtt = der(pzdt, 1) * dwt - dwtt * pzdt;

//     cvx = (pxds + pydt)*pwsqr - pxdst;
//     cvy = (pyds + pzdt)*pwsqr - pydst;
//     cvz = (pzds + pxdt)*pwsqr - pzdst;
//     cvw = pwq;

//     cvx = pxds + pzdt + px*pw;// * pzdt - pzds * pydt,
//     cvy = pyds + pxdt+py*pw;// * pxdt - pxds * pzdt,
//     cvz = pzds + pydt+pz*pw;// * pydt - pyds * pxdt;
//     cvw = pwsqr;

    // curvature: dvt^2 * dvss - 2 * dvs * dvt * dvst + dvs^2 * dvtt
    //Poly_2::NT _2(2.0, 0);

    //cvx = pxdt*pxdt*pxdss - _2*pxds*pxdt*pxdst + pxds*pxds*pxdtt,
    //cvy = pydt*pydt*pydss - _2*pyds*pydt*pydst + pyds*pyds*pydtt,
    //cvz = pzdt*pzdt*pzdss - _2*pzds*pzdt*pzdst + pzds*pzds*pzdtt; 

    //cvw = pw6;*/
    cv[0] = params[0], cv[1] = params[1], cv[2] = params[2], cv[3] = params[3];
}

/*!\brief
 * given a set point coordinates in parameter space (3rd coordinate ignored),
 * this function computes coordinates of points lying on the surface
 *
 * \c in and \c verts can point to the same container (in-place
 *  modification)
 */
void XSurface_arrangements::compute_spatial_coords(
    Point_vec_3f& in, Point_vec_3f& verts, Point_vec_3f *normals) const {
    
    Point_3::setup_parameterization(base_surfaces[base_index]);
    
    compute_equation();
    typedef Point_3::Poly_double_2 Poly_2;    
    const Poly_2 *params = cv;
    Poly_2 pnx, pny, pnz, pwq;    

    if(normals != 0) { // compute normals
        typedef CGAL::Polynomial_traits_d< Point_3::Poly_double_2 > P_traits;
        P_traits::Differentiate der;
    
        // assume pv(s, t) = (px(s, t); py(s, t); pz(s, t)) / pw(s, t)
        // compute dpv / ds and dpv / dt (s - innermost, t - outermost)
        Poly_2 pwsqr = params[3] * params[3], dws = der(params[3], 0), 
            dwt = der(params[3], 1);
        pwq = pwsqr * pwsqr;
    
        Poly_2 pxds = der(params[0], 0) * params[3] - dws * params[0],
               pyds = der(params[1], 0) * params[3] - dws * params[1],
               pzds = der(params[2], 0) * params[3] - dws * params[2];
    
        Poly_2 pxdt = der(params[0], 1) * params[3] - dwt * params[0],
               pydt = der(params[1], 1) * params[3] - dwt * params[1],
               pzdt = der(params[2], 1) * params[3] - dwt * params[2];
        // now compute cross product
        // pn = (dpv / ds) x (dpv / dt)
        pnx = pyds * pzdt - pzds * pydt,
        pny = pzds * pxdt - pxds * pzdt,
        pnz = pxds * pydt - pyds * pxdt;

        normals->resize(in.size());
    }

    CGAL::Polynomial_traits_d< Poly_2 >::Substitute subst;
    double coords[2];

    typedef Point_3f::value_type NT;
    verts.resize(in.size());
    for(unsigned i = 0; i < in.size(); i++) {
        
        coords[0] = in[i][0], coords[1] = in[i][1];
        
        double invw = 1.0 / subst(params[3], coords, coords+2);
        verts[i][0] = subst(params[0], coords, coords+2) * invw;
        verts[i][1] = subst(params[1], coords, coords+2) * invw;
        verts[i][2] = subst(params[2], coords, coords+2) * invw;

        if(normals != 0) {
            double xx = subst(pnx, coords, coords+2),
                yy = subst(pny, coords, coords+2),
                zz = subst(pnz, coords, coords+2);
            double invdet = 1.0/std::sqrt(xx*xx+yy*yy+zz*zz);

            (*normals)[i][0] = xx * invdet;
            (*normals)[i][1] = yy * invdet;
            (*normals)[i][2] = zz * invdet;
        }
    }
}

//! computes parameterization of tube / outer circles
//! \c which = 0: outer (x-coordinate taken), \c which = 1: tube (y-coordinate)
void XSurface_arrangements::compute_aux_coords(
    Point_vec_3f& in, Point_vec_3f& verts, int which) const {

    Point_3::setup_parameterization(base_surfaces[base_index]);

    const Point_3::Poly_double_1 *params = 
        (which == 0 ? Point_3::param_outer_circle() : 
            Point_3::param_tube_circle());

    typedef Point_3f::value_type NT;
    verts.resize(in.size());
    for(unsigned i = 0; i < in.size(); i++) {

        NT c = (which == 0 ? in[i][0] : in[i][1]);
        double invw = 1.0 / params[3].evaluate(c);
        verts[i][0] = params[0].evaluate(c) * invw,
        verts[i][1] = params[1].evaluate(c) * invw, 
        verts[i][2] = params[2].evaluate(c) * invw;
    }
}

void XSurface_arrangements::get_pole_coords(Point_3f& res) const {

    Point_3::setup_parameterization(base_surfaces[base_index]);
    const Point_3::Point_double_4& pole = Point_3::param_pole();
    
    double inv = 1.0 / pole[3];
    res[0] = pole[0] * inv, res[1] = pole[1] * inv, res[2] = pole[2] * inv;
}

double XSurface_arrangements::get_base_surf_extent() const {
    if(base_index < base_surfaces.size())
        return CGAL::to_double(base_surfaces[base_index].radius_a());
    return 0;
}

void XSurface_arrangements::get_base_surf_origin(Point_3d& origin) const {
    if(base_index < base_surfaces.size()) {
        // stupid way to store point coordinates in a vector...
        Base_surface_3::Vector vvv = base_surfaces[base_index].torus_center();
        origin[0] = CGAL::to_double(vvv[0]);
        origin[1] = CGAL::to_double(vvv[1]);
        origin[2] = CGAL::to_double(vvv[2]);
    }
}

typedef CGAL::Curve_renderer_facade< Geo_traits > CR_facade;

typedef std::pair< double, double > Coord_2;
typedef std::vector< Coord_2 > Coord_vec_2;

int n_pts_rendered = 0, n_pts_dropped = 0;

    /*!\brief
     * computes arc's approximation using preset window and resolution
     *
     * @note: call Dupin_cyclide_point_2::setup_renderer() and
     * setup_parameterization() before computing approximation
     */
void compute_approximation(Point_vec_3d& vec,
        const std::list< Coord_vec_2 >& points, double threshold)  {
    
    const Point_3::Poly_double_2* params = Point_3::param_surface();
    double invw0, cur_x = -1e100;
    
        std::list<Coord_vec_2>::const_iterator lit = points.begin();  

        int ii = 0;
        Point_3d *prev = 0, cur;
        bool last_dropped = false;
        double coords[2];
        CGAL::Polynomial_traits_d< Point_3::Poly_double_2 >::Substitute subst;

        while(lit != points.end()) {

            const Coord_vec_2& tmp = *lit++;
            for(Coord_vec_2::const_iterator cit = tmp.begin(); 
                cit != tmp.end(); cit++, ii++) {
                
                coords[0] = cit->first, coords[1] = cit->second;
                if(coords[0] < cur_x) {
                    // ignore discrepancies if drawing windows do
                    continue;   // not match
                }
                n_pts_rendered++;
                cur_x = coords[0];
                
                invw0 = 1.0 / subst(params[3], coords, coords+2);
                cur[0] = subst(params[0], coords, coords+2) * invw0,
                cur[1] = subst(params[1], coords, coords+2) * invw0,
                cur[2] = subst(params[2], coords, coords+2) * invw0;
                
                if(prev != 0) {
                    double dx = cur[0] - (*prev)[0], dy = cur[1] - (*prev)[1], 
                           dz = cur[2] - (*prev)[2], dist = dx*dx+dy*dy+dz*dz;
                    if(dist < threshold) {
                        last_dropped = true;
                        n_pts_dropped++;
                        continue;
                    }
                }

                //std::cerr << "[" << coords[0] << "; " << coords[1] << "]\n";
                last_dropped = false;
                vec.push_back(cur);
                prev = &vec.back();
            }
        }
        if(prev != 0 && last_dropped) {
            
            if(vec.size() < 2u) // if less than 2 points -> add the last one
                vec.push_back(cur);
            else // otherwise replace previous point by the last
                vec.back() = cur;
        }
        //std::cerr << "points dropped: " << n_drops << "\n";
}
 
/*!\brief
 * renders spatial points and arcs in the rectangular box \c bbox with
 * resolution \c res_w by \c res_h
 */
void XSurface_arrangements::render(const CGAL::Bbox_2& bbox,
        int res_w, int res_h) {

    if(base_index >= base_surfaces.size())
        return;

    // collects ids of already rasterized points
    std::map< std::size_t, unsigned > pt_id_map; // pt.id -> approx index
    xarcs_approx.clear();
    points_approx.clear();

    Point_3::setup_parameterization(base_surfaces[base_index]);
    double threshold = get_base_surf_extent() / 400.0;

    double extx = 10, exty = 10; // that is 3w x 3h
    double w = bbox.xmax() - bbox.xmin(), h = bbox.ymax() - bbox.ymin(),
           aw = (extx-1)*w/2, bh = (exty-1)*h/2;
    
    const int n_windows = 5;
    CGAL::Bbox_2 bboxes[5] = { bbox,
        CGAL::Bbox_2(bbox.xmin() - aw, bbox.ymin() - bh,
            bbox.xmin() + w, bbox.ymin()),
        CGAL::Bbox_2(bbox.xmin() + w, bbox.ymin() - bh,
            bbox.xmin() + w + aw, bbox.ymin() + h),
        CGAL::Bbox_2(bbox.xmin(), bbox.ymin() + h,
            bbox.xmin() + w + aw, bbox.ymin() + h + bh),
        CGAL::Bbox_2(bbox.xmin() - aw, bbox.ymin(),
            bbox.xmin(), bbox.ymin() + h + bh)
    };
    
    xarcs_approx.reserve(arr_xarcs.size());
    points_approx.reserve(arr_points.size());
    
    typedef std::list< Coord_vec_2 > Coord_vec_2_list; 
    std::vector< Coord_vec_2_list > lumps(arr_xarcs.size());

    CGAL::set_pretty_mode(std::cerr);
    std::cerr << "rendering.. density threshold: " << threshold << "\n";
    std::cerr << "narcs: " << arr_xarcs.size() << "\n";
    std::cerr << "npts : " << arr_points.size() << "\n";

    CGAL::Timer tm;

    tm.start();

    n_pts_rendered = 0, n_pts_dropped = 0;
    double res_factor = 10; // resolution factor
    for(int j = 0; j < n_windows; j++) {

    unsigned rw = res_w, rh = res_h;
    if(j > 0) {
        rw = (unsigned)((bboxes[j].xmax() - bboxes[j].xmin()) * rw /
                (w * res_factor));
        rh = (unsigned)((bboxes[j].ymax() - bboxes[j].ymin()) * rh / 
                (h * res_factor));
    }

    std::cerr << "\nrendering with box: " << bboxes[j] << "; resolution: " <<
        rw << " x " << rh << std::endl;
    Point_3::setup_renderer(bboxes[j], rw, rh);

    if(j == 0) { // render points with main resolution
        std::cerr << "processing event points..\n";

        //int idx;
        for(Points_3::const_iterator pcit = arr_points.begin(); 
                pcit != arr_points.end(); pcit++) {

            Point_3d pt;
            if(pcit->compute_approximation(pt.begin()) != pt.begin()) {
                pt_id_map.insert(std::make_pair(pcit->id(), 
                    points_approx.size()));
                points_approx.push_back(pt);
            }
        }
    }

    typeof(&compute_approximation)_(compute_approximation);

    std::set< int > iddds;

    int i = 0;
    for(XArcs_3::const_iterator xit = arr_xarcs.begin(); 
            xit != arr_xarcs.end(); xit++, i++) {

//         if(i != 17 ) continue;
        
        boost::optional< Coord_2 > ept1, ept2;
        Coord_vec_2_list points;
        
        if(j == n_windows - 1)
            CR_facade::instance().draw(*xit, points, &ept1, &ept2);
        else
            CR_facade::instance().draw(*xit, points);

        Coord_vec_2_list& lump = lumps[i];
        if(!points.empty()) {
            // merge lump & points sorted by x-coordinate
            Coord_vec_2_list::iterator cvl1 = lump.begin(), 
                cvl2 = points.begin();    
            while(cvl1 != lump.end() && cvl2 != points.end()) {
                if(cvl2->front().first < cvl1->front().first) {
                    lump.insert(cvl1, *cvl2); // insert before cvl1
                    cvl2++;
                } else
                    cvl1++;
            }
            std::copy(cvl2, points.end(), std::back_inserter(lump));
        }

        // boost::str( boost::format("[V%1%:%2%]") %  % (p) )
// std::cerr << boost::format("size of lump: %1% = %2%\n") % i % lump.size();
        
        if(j == n_windows - 1) { // last iteration: handle endpoints, colors
                                 // and draw the arc
            Point_vec_3d vec;

            Point_3 xitend = xit->curve_end(CGAL::ARR_MIN_END);
            if(xitend.location() != CGAL::ARR_INTERIOR) {
                
                Point_3d pt;
                typeof(pt_id_map.begin()) lookfor;
                if((lookfor = pt_id_map.find(xitend.id()))!= pt_id_map.end()) 
                    pt = points_approx[lookfor->second];
                else
                    xitend.compute_approximation(pt.begin());
                vec.push_back(pt);

            } else if(ept1) {
                if(lump.empty()) 
                    lump.push_back(Coord_vec_2(1));
                lump.front().front() = *ept1;
            }

            xitend = xit->curve_end(CGAL::ARR_MAX_END);
            if(xitend.location() == CGAL::ARR_INTERIOR && ept2) {
                if(lump.empty()) 
                    lump.push_back(Coord_vec_2());
                lump.back().push_back(*ept2); 
            }
            
            if(iddds.insert(xit->id()).second) {
                //CGAL::set_ascii_mode(std::cerr);
                 std::cerr << *xit << "\n\n";
            }
            if(!lump.empty()) {
                _(vec, lump, threshold);
            } else {
//                 std::cerr << "arc = " << i << "; empty! " <<
//                     *xit << "\n";
            }

            //std::cerr << "arc: " << *xit << "\n\n";

            if(xitend.location() != CGAL::ARR_INTERIOR) {
                Point_3d pt;
                typeof(pt_id_map.begin()) lookfor;
                
                if((lookfor = pt_id_map.find(xitend.id())) != pt_id_map.end())
                    pt = points_approx[lookfor->second];
                else
                    xitend.compute_approximation(pt.begin());
                vec.push_back(pt);
            }
            xarcs_approx.push_back(vec);
        }
    }
    }

    tm.stop();
    
    double efc = 100.0*(1.0 - (double)n_pts_dropped / n_pts_rendered);

    std::cerr << "\ntotal points rendered: " << n_pts_rendered << 
        "; points dropped: " << n_pts_dropped << "; EFFICIENCY: " <<
         efc << " %\n" << tm.time() << "\n";;
}



