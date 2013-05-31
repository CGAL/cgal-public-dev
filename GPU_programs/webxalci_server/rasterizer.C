// ============================================================================
//
// Copyright (c) 2001-2006 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of EXACUS (http://www.mpi-inf.mpg.de/projects/EXACUS/);
// you may redistribute it under the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with EXACUS.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// ----------------------------------------------------------------------------
//
// Library       : AlciX
// File          : demos/webxalci/rasterizer.C
// SoX_release   : $Name:  $
// Revision      : $Revision: 1.14 $
// Revision_date : $Date: 2007/02/06 12:03:04 $
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

#define NDEBUG 1

#include <pthread.h>
#include <qapplication.h>
#include <qtimer.h>
#include <qpixmap.h>
#include <qpainter.h>
#include <qimage.h>

#include <arpa/inet.h>
#include <netdb.h>

#include <iostream>
#include <vector>

// #define CGAL_NO_LEDA

#include "include/IPC.h"
#include "include/CGAL_includes.h"
#include "include/rasterizer.h"

#include <CGAL/Curved_kernel_via_analysis_2/Curve_renderer_facade.h>

//#include <SoX/GAPS/Gfx_wrapper.h>
//#include <SoX/GAPS/Subdivision_1.h>

#include <CGAL/Arr_walk_along_line_point_location.h>
//#include <CGAL/Arr_naive_point_location.h>

typedef CGAL::Curve_renderer_interface<CKvA_2, CGAL::Interval_nt<true> >
         Gfx_wrapper_inst;

struct Coord_2 {
    Coord_2() { }
    Coord_2(int x_, int y_) : x(x_), y(y_) { }
    
    int x, y;
};

// typedef std::pair< int, int > Coord_2;
typedef std::vector< Coord_2 > Coord_vec_2;
// arc's end-point coordinates
typedef std::pair<Coord_2, Coord_2> Point_pair_2;
// a list of vectors of point coordinates
typedef std::list<Coord_vec_2> List_vector_2;
    
const static QColor Highlight_color = QColor(40, 40, 128);

const static QColor ccv(QRgb(0x123123));

static QColor rasterizer_colors[] = {
    QColor(QRgb(0xcc0000)),
    QColor(QRgb(0x0099ff)),
    QColor(QRgb(0xff0099)),
    QColor(QRgb(0xff3300)),
    QColor(QRgb(0xb3b300)),
    QColor(QRgb(0x00cc00)),
    QColor(QRgb(0x009999)),
    QColor(QRgb(0x0000ff)),
    QColor(QRgb(0x9900cc))
};
static int n_rast_colors = sizeof(rasterizer_colors) / sizeof(QColor);

static QColor aux_colors[] = {
    QColor(QRgb(0x6495ed)),
    QColor(QRgb(0xcd5c5c)),
    QColor(QRgb(0x00bfff)),
    QColor(QRgb(0xfa8072)),
    QColor(QRgb(0x7b68ee)),
    QColor(QRgb(0xff1493)),
    QColor(QRgb(0x483d8b)),
    QColor(QRgb(0xff69b4)),
    QColor(QRgb(0x1e90ff)),
    QColor(QRgb(0xffa500)),
    QColor(QRgb(0x87cefa)),
    QColor(QRgb(0x00fa9a)),
    QColor(QRgb(0xadd8e6)),
    QColor(QRgb(0x2e8b57))
};
static int n_aux_colors = sizeof(aux_colors) / sizeof(QColor);

// QColor(189, 0, 0),
//     QColor(138, 43, 226),
//     QColor(238, 18, 137),
//     QColor(95, 158, 160),
//     QColor(205, 55, 0),
//     QColor(0, 206, 209),
//     QColor(255, 127, 0),
//     QColor(85, 107, 47),
//     QColor(139, 0, 139),
//     QColor(0, 189, 0),
//     QColor(205, 149, 12),
//     QColor(0, 0, 189),
//     QColor(238, 99, 99)

// representation of infinity symbols
const static char *minus_inf = "-oo";
const static char *plus_inf = "%2Boo";

extern pthread_mutex_t painter_mtx;
extern pthread_mutex_t arcs_cache_mtx;

void plot_arc(const Arc_2& arc, Gfx_wrapper_inst& renderer,
        Rasterize_mode mode, QColor color, void *dst);
        
void plot_point(const Point_2& pt, Gfx_wrapper_inst& renderer,
        Rasterize_mode mode, QColor color, void *dst);

void plot_arc_impl(const List_vector_2& points, const Point_pair_2& end_points,
        Rasterize_mode mode, QColor color, QPixmap *plot);

void print_arc_impl(const List_vector_2& points, const Point_pair_2& end_points,
        Rasterize_mode mode, std::string *print_out);

void plot_point_impl(const Coord_2& coord, Rasterize_mode mode, QColor color,
    QPixmap *plot);

void print_point_impl(const Coord_2& coord, Rasterize_mode mode,
    std::string *print_out);

int Rasterizer::width = 0;
int Rasterizer::height = 0;

static Kernel_2 *kernel_2;

Rasterizer::Rasterizer() {

    // some horribly terrified sh*t...
    // somehow affects to ACK_2 caching
    kernel_2 = &CKvA_2::instance().kernel();
}

//! calls curve analysis routine
void Rasterizer::analyse_curve(Analysis_entry& analysis, 
    const Poly_int_vector& poly_vec) {

    CGAL_precondition("switch off preconditions!!!");

    //Kernel_2::Construct_curve_2 cc_2;
    typedef std::vector< Curve_analysis_2 > Curve_vector;
    Curve_vector curves(poly_vec.size());

//     Kernel_2 *kernel_2 = &CKvA_2::instance().kernel();

    int i = 0;
    typedef typeof(kernel_2->curve_cache_2()) Curve_cache_2;
    const Curve_cache_2& ccache = kernel_2->curve_cache_2();

    for(Poly_int_vector::const_iterator pit = poly_vec.begin();
        pit != poly_vec.end(); pit++, i++) {

        Curve_analysis_2 ccc;

        Curve_cache_2::Canonicalizer canonicalize;
        Poly_int2 key = canonicalize(*pit);
        pthread_mutex_lock(&arcs_cache_mtx);
        Curve_cache_2::Find_result p = ccache.find(key);
        pthread_mutex_unlock(&arcs_cache_mtx);
        if(!p.second) {
            Poly_int2 ff = key;
Lbegin:            
            try {
                Curve_cache_2::Creator create(kernel_2);
                ccc = create(ff);
                // this is to force lazy stuff get processed..
                volatile unsigned x = ccc.number_of_status_lines_with_event();

                pthread_mutex_lock(&arcs_cache_mtx);
                ccache.insert(Curve_cache_2::Data_type(ff, ccc));
                pthread_mutex_unlock(&arcs_cache_mtx);
                curves[i] = ccc;
                continue; 
            } catch(...) {
                Kernel_2::Decompose_2 d2 = kernel_2->decompose_2_object();
                ff = d2(ff);
                std::cerr << "restarting.. squarefree part: " << ff << "\n";
            }
            goto Lbegin;

        } else 
            curves[i] = (p.first)->second;
    }

    CGAL::insert(analysis.arr, curves.begin(), curves.end(),
            boost::false_type());
    analysis.arcs.clear();
    analysis.arcs.reserve(analysis.arr.number_of_isolated_vertices() +
            analysis.arr.number_of_edges());
    analysis.points.clear();
    analysis.points.reserve(analysis.arr.number_of_vertices());

    for(CGAL_Arrangement_2::Vertex_const_iterator vit =
        analysis.arr.vertices_begin(); vit != analysis.arr.vertices_end();
            vit++) {
        
        if(vit->is_isolated()) // collect isolated points only
            analysis.iso_points.push_back(vit->point());
        else
            analysis.points.push_back(vit->point());
    }

    for(CGAL_Arrangement_2::Edge_const_iterator eit =
        analysis.arr.edges_begin(); eit != analysis.arr.edges_end(); eit++) {
        
        analysis.arcs.push_back(eit->curve());
    }

    std::stringstream fbuf;
    typedef CGAL_Arrangement_2::Ccb_halfedge_const_circulator
        Ccb_circulator;

    CGAL_Arrangement_2::Face_const_iterator fit;
    for(fit = analysis.arr.faces_begin(), i = 0;
        fit != analysis.arr.faces_end(); fit++, i++) {
        if(i != 0)
            fbuf << '\n';
        if(fit->is_unbounded()) {
            int ccb;
            ccb = (fit->number_of_outer_ccbs() > 0 ? CGAL::circulator_size(
                Ccb_circulator(*fit->outer_ccbs_begin())) : 0);
            fbuf << "unbounded," << ccb;
        } else {
            fbuf << "bounded," << CGAL::circulator_size(fit->outer_ccb());
        }
        if(i >= MAX_FACES_PRINT_OUT)
            break;
    }
    analysis.faces_print = fbuf.str();
}

bool Rasterizer::locate_point(const Analysis_entry& analysis,
    const CGAL::Bbox_2& box, const Point_coords& coords,
        SHM_Point_query_reply *reply) {

    //std::cerr << "starting point location query: "
      //  << coords.x << " and " << coords.y << "\n";
    CGAL::Arr_walk_along_line_point_location<CGAL_Arrangement_2>
        point_query(analysis.arr);

    //CGAL::Arr_naive_point_location<CGAL_Arrangement_2>
      //   point_query(analysis.arr);

    if(coords.x < 0 || coords.x >= Rasterizer::width ||
            coords.y < 0 || coords.y >= Rasterizer::height)
        return false;

    double xf = box.xmin() + double(coords.x)*(box.xmax()-box.xmin())/
                    double(Rasterizer::width),
           yf = box.ymin() + double(Rasterizer::height-coords.y)*
                (box.ymax() - box.ymin())/double(Rasterizer::height);

    /*std::cerr << "box: [" << box.x_min << "; "<< box.x_max << "]x[" <<
        box.y_min << "; " << box.y_max << std::endl;
    std::cerr << Rasterizer::width << " and " << Rasterizer::height << "\n";
    std::cerr << "xf = " << xf << "; yf = " << yf << "\n";
    */            
    typedef CGAL::Fraction_traits<Rational> FTraits;
    Rational rat_x(xf), rat_y(yf);

    Integer num_y, den_y;
    FTraits::Decompose()(rat_y, num_y, den_y);

    // construct horiz line y*den_y - num_y = 0
    Poly_int2 poly(Poly_int1(-num_y), Poly_int1(den_y));

    Curve_analysis_2 curve(kernel_2, poly);
    Point_2 pt(X_coordinate_1(rat_x), curve, 0);
    CGAL::Object obj = point_query.locate(pt);

    CGAL_Arrangement_2::Vertex_const_handle    v;
    CGAL_Arrangement_2::Halfedge_const_handle  e;
    CGAL_Arrangement_2::Face_const_handle      f;

    reply->index = -1;
    reply->type = -1;
    if(CGAL::assign (f, obj)) {
        // q is located inside a face:
        if(f->is_unbounded()) {
            //std::cerr << "Inside the unbounded face.";
        } else {
            //std::cerr << "Inside a bounded face.";
        }
        reply->type = QUERY_ON_FACE;

        int i = 0;
        for(CGAL_Arrangement_2::Face_const_iterator fit =
            analysis.arr.faces_begin(); fit != analysis.arr.faces_end();
            fit++, i++) {
            if(fit == f) {
                //std::cerr << "got it: " << i << "\n";
                reply->index = i;
                break;
            }
        }
          
    } else if (CGAL::assign (e, obj)) {
        // q is located on an edge:
        //std::cerr << "On an edge.";
        reply->type = QUERY_ON_EDGE;
        
    } else if (CGAL::assign (v, obj)) {
        // q is located on a vertex:
        if (v->is_isolated()) {
            //std::cerr << "On an isolated vertex.";
        } else {
            //std::cerr << "On a vertex.";
        }
        reply->type = QUERY_ON_POINT;
        
    } else {
        //std::cerr << "Invalid object.\n";
        return false;
    }

    return true;
}

//! renders a set of curve segments to the \c plot;
//! \c box defines drawing window, \c indices specifies a set of arcs to be
//! rasterized
void Rasterizer::plot_arcs(const Analysis_entry& analysis, uint *indices,
    uint n_indices, const CGAL::Bbox_2& box, Rasterize_mode mode,
    void *dst)
{
    if(n_indices == 0) // draw only axes
        return;

    Gfx_wrapper_inst renderer;

    typedef std::map<std::size_t, int> Color_map;
    Color_map color_map;
    QColor color;

    if(mode != PRINT_COORDS)
        renderer.setup(box, width, height);
    else
        renderer.setup(box, width/2, height/2);
    
    // first come isolated vertices and then normal edges
    uint n_isolated = analysis.iso_points.size(),
         n_arcs = analysis.arcs.size(), ssize = n_isolated + n_arcs, n, i,
         idx, cindex;

    n = (indices[0] != -1u ? n_indices : ssize);
    color = Highlight_color;
    
    // n is either #of indices or total number of edges+isolated
    for(i = 0; i < n; i++) {
        idx = (indices[0] == -1u ? i : indices[i]);
        if(idx >= ssize) 
            continue;
        
        int id;
        if(idx < n_isolated) 
            id = analysis.iso_points[idx].curve().id();
        else
            id = analysis.arcs[idx-n_isolated].curve().id();

        if(mode == DRAW_DEFAULT || mode == DRAW_FLAT) {
            if(!analysis.one_curve) {
                // pick up the first available color index
                cindex = color_map.size();
                std::pair<Color_map::iterator, bool> ret =
                    color_map.insert(std::make_pair(id, cindex));
                if(!ret.second) // already exists in the map
                    cindex = ret.first->second;
            } else
                cindex = (mode == DRAW_DEFAULT ? idx : id);
            color = rasterizer_colors[cindex % n_rast_colors];
        }

        if(idx < n_isolated) 
            plot_point(analysis.iso_points[idx], renderer, mode, color, dst);
        else 
            plot_arc(analysis.arcs[idx-n_isolated], renderer, mode, color,
                 dst);
    }
}

void Rasterizer::plot_points(const Analysis_entry& analysis, uint *indices,
    uint n_indices, const CGAL::Bbox_2& box, Rasterize_mode mode,
    QPixmap *plot) {
    
    if(n_indices == 0) // draw only axes
        return;

    Gfx_wrapper_inst renderer;
    renderer.setup(box, width, height);
    
    // first come isolated vertices and then normal edges
    uint n_verts = analysis.points.size(), ssize = n_verts, n, i, idx;
    n = (indices[0] != -1u ? n_indices : ssize);

    for(i = 0; i < n; i++) {
        idx = (indices[0] == -1u ? i : indices[i]);
        if(idx >= ssize) 
            continue;
        plot_point(analysis.points[idx], renderer, mode, Qt::black, plot);
    }
}

void Rasterizer::plot_faces(const Analysis_entry& analysis, uint *indices,
    uint n_indices, const CGAL::Bbox_2& box, Rasterize_mode mode,
        QPixmap *plot) {

    if(n_indices == 0) // draw only axes
        return;

    Gfx_wrapper_inst renderer;
    renderer.setup(box, width, height);
    QColor color = Highlight_color;
    
    uint n_faces = analysis.arr.number_of_faces(), ssize = n_faces, n, i,
          idx;
    n = (indices[0] != -1u ? n_indices : ssize);

    CGAL_Arrangement_2::Inner_ccb_const_iterator icit, icend;
    CGAL_Arrangement_2::Outer_ccb_const_iterator ocit, ocend;
    
    typedef CGAL_Arrangement_2::Ccb_halfedge_const_circulator
        Ccb_circulator;
    Ccb_circulator circ, curr;
    CGAL_Arrangement_2::Face_const_iterator fit, fend =
        analysis.arr.faces_end();
        
    for(fit = analysis.arr.faces_begin(), i = 0; i < n; i++) {

        idx = (indices[0] == -1u ? i : indices[i]);
        if(idx >= ssize)
            continue;
        while(i < idx && fit != fend)
            i++, fit++;
        if(fit == fend)
            break;
        
//         if(fit->is_unbounded()) {
            /*std::cerr << "unbounded face, inner: " <<
                fit->number_of_inner_ccbs() << "\n";
            std::cerr << "outer: " << fit->number_of_outer_ccbs() << "\n";*/
            for(icit = fit->inner_ccbs_begin(),
                icend = fit->inner_ccbs_end();
                    icit != icend; icit++) {
                circ = Ccb_circulator(*icit);
                curr = circ;
                do {
                    if(!curr->is_fictitious())
                     plot_arc(curr->curve(), renderer, mode, color, plot);
                //std::cout << "arc: " << curr->curve().support().id() << "\n";
                } while(++curr != circ);
            }
//             continue;
//         }

        for(ocit = fit->outer_ccbs_begin(), ocend = fit->outer_ccbs_end();
                    ocit != ocend; ocit++) {
            circ = Ccb_circulator(*ocit);
            curr = circ;
            do {
                if(!curr->is_fictitious())
                  plot_arc(curr->curve(), renderer, mode, color, plot);
                //std::cout << "arc: " << curr->curve().support().id() << "\n";
            } while(++curr != circ);
        }
    }
}

void plot_arc(const Arc_2& arc, Gfx_wrapper_inst& renderer,
    Rasterize_mode mode, QColor color, void *dst) {

    List_vector_2 points;
    
    boost::optional < Coord_2 > p1, p2;
    renderer.draw(arc, points, &p1, &p2);
   
    Point_pair_2 end_points(*p1, *p2);

    if(points.empty())  // to ensure that colors are not alternated
        return;

    if(mode == PRINT_COORDS)
        print_arc_impl(points, end_points, mode, (std::string *)dst);
    else
        plot_arc_impl(points, end_points, mode, color, (QPixmap *)dst);
}

void plot_point(const Point_2& pt, Gfx_wrapper_inst& renderer,
     Rasterize_mode mode, QColor color, void *dst) {

    Coord_2 coord;
    if(!renderer.draw(pt, coord))
        return;
    
    if(mode == PRINT_COORDS)
        print_point_impl(coord, mode, (std::string *)dst);
    else
        plot_point_impl(coord, mode, color, (QPixmap *)dst);
}

void plot_arc_impl(const List_vector_2& points, const Point_pair_2& end_points,
    Rasterize_mode mode, QColor color, QPixmap *plot) {

    QPainter pnt;
    pthread_mutex_lock(&painter_mtx);
    pnt.begin(plot);

    int pen_w = (mode == HIGHLIGHT_ARCS || mode == HIGHLIGHT_FACES ? 4 : 2);
    int h = Rasterizer::height;
    pnt.setPen(QPen(color, pen_w, Qt::SolidLine, Qt::RoundCap,
        Qt::MiterJoin));
                
    List_vector_2::const_iterator lit = points.begin();
    while(lit != points.end()) {
        const Coord_vec_2& pvec = *lit++;
        Coord_vec_2::const_iterator vit = pvec.begin();

        pnt.moveTo(vit->x, h - vit->y);  
        if(pvec.size() == 2) {
            vit++;
            pnt.lineTo(vit->x, h - vit->y);
        } else {
            while(vit != pvec.end()) {
                pnt.lineTo(vit->x, h - vit->y);
                vit++;
            }
        }
    }
            
    // do not draw end-points in one-color mode
    if(mode == DRAW_DEFAULT) {
        QBrush b1(Qt::black);
        pnt.setBrush(b1);
        pnt.setPen(QPen(Qt::NoPen));
        pnt.drawEllipse(end_points.first.x - 3,
             h - end_points.first.y - 3, 6, 6);
        pnt.drawEllipse(end_points.second.x - 3,
             h - end_points.second.y - 3, 6, 6);
    }
    pnt.end();
    pthread_mutex_unlock(&painter_mtx);
}

// just prints out point coordinates
void print_arc_impl(const List_vector_2& points,
    const Point_pair_2& end_points, Rasterize_mode mode,
    std::string *print_out) {

    if(print_out->length() >= MAX_COORDS_PRINTOUT)
        return;

    std::stringstream sbuf;
    int h = Rasterizer::height/2;
    sbuf.put('{');
    List_vector_2::const_iterator lit = points.begin();
    while(lit != points.end()) {
        const Coord_vec_2& pvec = *lit;
        Coord_vec_2::const_iterator vit = pvec.begin();

        sbuf << (vit->x) << '|' << (h - vit->y);
        vit++;
        while(vit != pvec.end()) {
            sbuf << ',' << (vit->x) << '|' << (h - vit->y);
            vit++;
        }
        lit++;
    }
    sbuf.put('}');
    sbuf.put('\n');
    *print_out += sbuf.str();
}

void plot_point_impl(const Coord_2& coord, Rasterize_mode mode, QColor color,
        QPixmap *plot) {
    
    QPainter pnt;
    pthread_mutex_lock(&painter_mtx);
    pnt.begin(plot);

    uint sz = (mode == HIGHLIGHT_ARCS || mode == HIGHLIGHT_VERTS ? 5 : 3);
    int h = Rasterizer::height;
    QBrush b1(color);
    pnt.setBrush(b1);
    pnt.setPen(QPen(Qt::NoPen));
    pnt.drawEllipse(coord.x - sz, h - coord.y - sz, sz*2, sz*2);
    pnt.end();
    pthread_mutex_unlock(&painter_mtx);
}

void print_point_impl(const Coord_2& coord, Rasterize_mode mode,
    std::string *print_out) {

    if(print_out->length() >= MAX_COORDS_PRINTOUT)
        return;

    std::stringstream sbuf;
    int h = Rasterizer::height/2;
    sbuf << '{' << (coord.x) << '|' << (h - coord.y) << '}';
    *print_out += sbuf.str();
}

void Rasterizer::print_endpoint(const Arc_2& arc,
    CGAL::Arr_curve_end end, std::ostream& os) {

    CGAL::Arr_parameter_space loc = arc.location(end);

    if(loc == CGAL::ARR_INTERIOR) {
        std::pair< double, double > xy =
            arc.curve_end(end).xy().to_double();
        os << xy.first << ',' << xy.second;
        return;
    }

    if(loc == CGAL::ARR_LEFT_BOUNDARY)
        os << minus_inf;
    else if(loc == CGAL::ARR_RIGHT_BOUNDARY)
        os << plus_inf;
    else 
        os << CGAL::to_double(arc.curve_end_x(end)) << " (asympt)";

    os << ',';
    switch(loc) {
    case CGAL::ARR_BOTTOM_BOUNDARY:
        os << minus_inf;
        return;

    case CGAL::ARR_TOP_BOUNDARY:
        os << plus_inf;
        return;
        
    case CGAL::ARR_LEFT_BOUNDARY:
    case CGAL::ARR_RIGHT_BOUNDARY: {

        CGAL::Object obj =
            arc.curve().asymptotic_value_of_arc(loc, arc.arcno());
            
        CGAL::Arr_parameter_space loc2;
        if(CGAL::assign(loc2, obj)) 
            os << (loc2 == CGAL::ARR_BOTTOM_BOUNDARY ? minus_inf : plus_inf);
        else {
            X_coordinate_1 y;
            if(CGAL::assign(y, obj)) 
                os << CGAL::to_double(y) << " (asympt)";
            else {
                os << ' ';
                std::cerr << "Ill-typed object found..";
            }
        }
        return;
        }    
    default:
        os << arc.curve_end(end).arcno();
    }
}

//! prints out a set of arcs, vertices and faces
void Rasterizer::print_out(Analysis_entry& analysis) {
    
    typedef std::map<int, int> Index_map;
    Index_map index_map;
    Index_map::iterator imap;
    std::vector<std::pair<double, int> > approx;
    int i, len;
    double fx;
    int fy = 0;
    std::stringstream vbuf, abuf;
    Points_2::const_iterator vit;

    for(vit = analysis.points.begin(), i = 0; vit != analysis.points.end();
        vit++, i++) {

        //std::cerr << "printing out point: " << i << "\n";
        if(i != 0)
            vbuf << '\n';
        // need to save points in an array since they will be required later
//         fx = CGAL::to_double(vit->x());
//         event.refine_to(vit->arcno(), Rational(1,20000));
//         fy = vit->arcno();
        std::pair< double, double > xy = vit->xy().to_double();

        //approx.push_back(std::make_pair(fx, fy));
        //index_map.insert(std::make_pair(vit->id(), i));

        //vbuf << approx[i].first << ',' << approx[i].second;
        vbuf << xy.first << ',' << xy.second;
        if(i >= MAX_VERTICES_PRINT_OUT || vbuf.str().length() >=
                MAX_SEG_SIZE-200)
            break;
    }
    analysis.verts_print = vbuf.str();
    len = vbuf.str().length();

    // no isolated points
    for(vit = analysis.iso_points.begin(), i = 0;
        vit != analysis.iso_points.end(); vit++, i++) {

        if(i != 0)
            abuf << '\n';
        // need to save points in an array since they will be required later
        std::pair< double, double > xy = vit->xy().to_double();

        // treat this point as being isolated
        abuf << xy.first << ',' << xy.second;
        // << ',' << fx << ',' << fy << ','
        if(i >= MAX_ARCS_PRINT_OUT || len + abuf.str().length() >=
                (int)MAX_SEG_SIZE-200)
            break;
    }

    Arcs_2::const_iterator eit, eend = analysis.arcs.end();
    // put arcs to the same buffer as isolated points
    for(eit = analysis.arcs.begin(), i = 0; eit != eend; eit++, i++) {

        // either not the 1st arc or isolated != 0
        if(i != 0 || analysis.iso_points.size() != 0)
             abuf << '\n';

        print_endpoint(*eit, CGAL::ARR_MIN_END, abuf);
        abuf << ',';
        print_endpoint(*eit, CGAL::ARR_MAX_END, abuf);
        abuf << ',';
        if(!eit->is_vertical())
            abuf << eit->arcno();
        else
            abuf << "-1";
        // ensure that we do not run out of segment's size
        if(i >= MAX_ARCS_PRINT_OUT ||
                len + abuf.str().length() >= (int)MAX_SEG_SIZE-200)
            break;
    }
    analysis.arcs_print = abuf.str();
   /* std::cout << "verts: " << analysis.verts_print << "\n\n" <<
        analysis.arcs_print << "\n";*/
}

struct Approx_point {
    Approx_point() { }
    Approx_point(const std::pair< double, double >& p_,
            int flag_, int cindex_) : xy(p_), flag(flag_),
                cindex(cindex_) { }

    std::pair< double, double > xy;
    int flag; // 00 - both interior; 10 - x special, y interior;
              // 01 - x interior, y special; 11 - both special ??
    int cindex; // color index
};

inline void set_minmax(bool& set,
        const double& coord, double& cmin, double& cmax) {

    if(!set) {
        cmin = cmax = coord;
        set = true;
    } else if(cmin > coord)
        cmin = coord;
    else if(cmax < coord)
        cmax = coord;
}

void graph_get_endpoint(const Arc_2& arc, CGAL::Arr_curve_end end,
        Approx_point& p, int w, int h) {

    CGAL::Arr_parameter_space loc = arc.location(end);

    p.flag = 0; // neither are special
    switch(loc) {

    case CGAL::ARR_BOTTOM_BOUNDARY:
    case CGAL::ARR_TOP_BOUNDARY:
        p.xy.first = arc.curve_end_x(end).to_double();
        p.xy.second = (loc == CGAL::ARR_BOTTOM_BOUNDARY ? 0 : h);
        p.flag = 1; // y special
        break;
    case CGAL::ARR_LEFT_BOUNDARY:
    case CGAL::ARR_RIGHT_BOUNDARY: {
        p.xy.first = (loc == CGAL::ARR_LEFT_BOUNDARY ? 0 : w);
        p.flag = 2; // x special
        CGAL::Object obj =
            arc.curve().asymptotic_value_of_arc(loc, arc.arcno());
            
        CGAL::Arr_parameter_space loc2;
        if(CGAL::assign(loc2, obj)) {
            p.xy.second = (loc2 == CGAL::ARR_BOTTOM_BOUNDARY ? 0 : h);
            p.flag = 3; // both special
        } else {
            Kernel_2::Coordinate_1 y;
            if(CGAL::assign(y, obj)) {
                p.xy.second = CGAL::to_double(y);
            } else
                CGAL_error_msg("Ill-typed object returned..\n");
        }
        break;
        }
//     case CGAL::ARR_INTERIOR:
    default:
        p.xy = arc.curve_end(end).xy().to_double();
        break;
    }
}

void Rasterizer::draw_topology_graph(const Analysis_entry& analysis,
        const CGAL::Bbox_2& box, void *dst) {

    int w = Rasterizer::width, h = Rasterizer::height;
    uint n_isolated = analysis.iso_points.size(),
         n_arcs = analysis.arcs.size();
    
    double xmin, ymin, xmax, ymax;
    bool x_set = false, y_set = false;;

    Curve_analysis_2::Status_line_1 cv_line;
    typedef Curve_analysis_2::Coordinate_2 Coordinate_2;

    // every 3 consecutive points correspond to one arc
    std::vector< Approx_point > arcs_xy;
    // coordinates of isolated points
    std::vector< Approx_point > pts_xy;

    AK_1::Algebraic_real_traits::Lower_bound lbound_x;
    AK_1::Algebraic_real_traits::Upper_bound ubound_x;

    int i, j;

    typedef std::map<std::size_t, int> Color_map;
    Color_map color_map;

    for(i = 0; i < n_isolated; i++) {

        const Point_2& pt = analysis.iso_points[i];
        Approx_point p;
        int id = pt.curve().id();
        p.cindex = color_map.size();
        std::pair<Color_map::iterator, bool> ret =
            color_map.insert(std::make_pair(id, p.cindex));
        if(!ret.second) // already exists in the map
            p.cindex = ret.first->second;

        p.xy = pt.xy().to_double();
        set_minmax(x_set, p.xy.first, xmin, xmax);
        set_minmax(y_set, p.xy.second, ymin, ymax);
        pts_xy.push_back(p);
    }

    for(i = 0; i < n_arcs; i++) {

        Approx_point p;
        const Arc_2& arc = analysis.arcs[i];
        int id = arc.curve().id();
        p.cindex = color_map.size();
        std::pair<Color_map::iterator, bool> ret =
            color_map.insert(std::make_pair(id, p.cindex));
        if(!ret.second) // already exists in the map
            p.cindex = ret.first->second;

        // handle arc's end-points
        CGAL::Arr_curve_end end = CGAL::ARR_MIN_END;
        for(j = 0; j < 3; j++) {

            if(j % 2 == 0) {
                graph_get_endpoint(arc, end, p, w, h);
            // j == 1: handle intermediate point
            } else if(arc.is_vertical()) {
                    // need special processing here !!!
// distinguish cases: whether one or both end-points are finite or not..
                p.xy.second = h * 0.5; // leave the first coordinate unchanged
                p.flag = 0x1; // set y-coordinate special
            } else {

                // this branch only for middle points
                CGAL::Arr_parameter_space loc =
                        arc.location(CGAL::ARR_MAX_END);
                bool max_end_x = (loc == CGAL::ARR_INTERIOR ||
                        loc == CGAL::ARR_BOTTOM_BOUNDARY ||
                        loc == CGAL::ARR_TOP_BOUNDARY);

                if((p.flag & 0x2) == 0) { // indicates that point is finite
                    const X_coordinate_1& c1 =
                        arc.curve_end_x(CGAL::ARR_MIN_END);

                    if(/*!analysis.one_curve &&*/ max_end_x) {
                        const X_coordinate_1&
                          c2 = arc.curve_end_x(CGAL::ARR_MAX_END);
                        Rational bnd = Kernel_2::Bound_between_1()(c1, c2);
                        cv_line = arc.curve().status_line_at_exact_x(bnd);
                    } else {
                        // when one_curve || !max_end_x
                        cv_line = arc.curve().status_line_at_exact_x(
                                ubound_x(c1) + Rational(1));
                    }
                // otherwise try the right boundary
                } else if(max_end_x) {
                    const X_coordinate_1&
                          c2 = arc.curve_end_x(CGAL::ARR_MAX_END);
                    cv_line = arc.curve().status_line_at_exact_x(
                        lbound_x(c2) - Rational(1));
                } else 
                    cv_line = arc.curve().status_line_of_interval(0);
                    
                p.xy = cv_line.xy_coordinate_2(arc.arcno()).to_double();
                p.flag = 0;
            }
            if((p.flag & 0x2) == 0)
                set_minmax(x_set, p.xy.first, xmin, xmax);
            if((p.flag & 0x1) == 0)
                set_minmax(y_set, p.xy.second, ymin, ymax);
            arcs_xy.push_back(p);
            end = CGAL::ARR_MAX_END;
        }
    }

// TODO: what if some boundaries are undefined ??
    printf("[%.4f; %.4f] -[%.4f; %.4f]\n", xmin, ymin, xmax, ymax);
    printf("x_set: %d; y_set: %d\n", x_set, y_set);

    if(!x_set) {
        printf("FATAL: x-range is not set! bailing out..\n");
        return;
    }
    if(!y_set) { // when a single vertical line, y-range is not set
        ymin = ymax = 0.0; 
    }

    double diffx = xmax - xmin, diffy = ymax - ymin;
    if(diffx <= 1e-10) {
        if(diffy <= 1e-10) {
            diffy = 2.0;
            ymax = ymin + 1.0;
            ymin = ymin - 1.0;
        }
        xmax = xmin + diffy * 0.5;
        xmin = xmin - diffy * 0.5;
    } else if(diffy <= 1e-10) {
        ymax = ymin + diffx * 0.5;
        ymin = ymin - diffx * 0.5;
    }
    printf(" [%.4f; %.4f] -[%.4f; %.4f]\n", xmin, ymin, xmax, ymax);

    double enlarge_f = 0.05; // enlarge window by 5 %
    double xinc = enlarge_f*(xmax - xmin), yinc = enlarge_f*(ymax - ymin);
    xmin -= xinc, xmax += xinc;
    ymin -= yinc, ymax += yinc;

    double x_min = box.xmin(), x_max = box.xmax(),
            y_min = box.ymin(), y_max = box.ymax();
        
//     double inv_w = w / (xmax - xmin), inv_h = h / (ymax - ymin);
    double inv_w = w / (x_max - x_min),
            inv_h = h / (y_max - y_min);

    QPainter pnt;
    QPixmap *plot = (QPixmap *)dst;
    pthread_mutex_lock(&painter_mtx);
    pnt.begin(plot);

//     pnt.setPen(QPen(QColor(163, 209, 255), 3, Qt::SolidLine, Qt::RoundCap,
//             Qt::MiterJoin));

    for(i = 0; i < arcs_xy.size(); i++) {
        int x, y, flag = arcs_xy[i].flag;        

         const std::pair< double, double >& xy =
                arcs_xy[i].xy;

        if(flag & 0x2) // x is special
            x = (int)xy.first;
        else
            x = (int)((xy.first - x_min) * inv_w);

        if(flag & 0x1) // y is special
            y = (int)xy.second;
        else
            y = (int)((xy.second - y_min) * inv_h);

        pnt.setPen(QPen(aux_colors[arcs_xy[i].cindex % n_aux_colors], 3,
            Qt::SolidLine, Qt::RoundCap, Qt::MiterJoin));

        if(i % 3 == 0) // 3 points for each arc
            pnt.moveTo(x, h - y);
        else
            pnt.lineTo(x, h - y);
    }

    pnt.setPen(QPen(QColor(QRgb(0x446644)), 2, Qt::SolidLine, Qt::RoundCap,
            Qt::MiterJoin));
//     pnt.setBrush(QBrush(QColor(226, 221, 165)));

    int rad = 5;
    for(i = 0; i < arcs_xy.size(); i++) {

        if((i-1) % 3 == 0)  // skip middle points
            continue;

        int x, y, flag = arcs_xy[i].flag;

        if(flag & 0x3) // skip special points
            continue;

         const std::pair< double, double >& xy =
                arcs_xy[i].xy;

        x = (int)((xy.first - x_min) * inv_w);
        y = (int)((xy.second - y_min) * inv_h);

        pnt.setBrush(QBrush(aux_colors[arcs_xy[i].cindex % n_aux_colors]));
        pnt.drawEllipse(x - rad, h - y - rad, rad*2, rad*2);
    }

    for(i = 0; i < pts_xy.size(); i++) {
        int x, y;
        x = (int)((pts_xy[i].xy.first - x_min) * inv_w);
        y = (int)((pts_xy[i].xy.second - y_min) * inv_h);

        pnt.setBrush(QBrush(aux_colors[pts_xy[i].cindex % n_aux_colors]));
        pnt.drawEllipse(x - rad, h - y - rad, rad*2, rad*2);
    }

    pnt.end();
    pthread_mutex_unlock(&painter_mtx);
}

//! plots a curve equation using 1D subdivision method
void Rasterizer::plot_subdivision(const Poly_int_vector& poly_vec,
    const CGAL::Bbox_2& box, QPixmap *plot)
{
    /*Subdivision_1<NiX::Interval, AC2> subdiv_renderer;
    QPainter pnt;
    
    subdiv_renderer.setup(box.x_min, box.y_min, box.x_max, box.y_max, width,
        height);
    //////// currently only the 1st poly is rasterized //////////////////
    subdiv_renderer.set_polynomial(poly_vec[0]);
                
    typedef std::pair<int, int> Int_pair;
    std::list<Int_pair> points;
           
    subdiv_renderer.draw(std::back_inserter(points));
    if(points.empty()) 
         return;
    
    pthread_mutex_lock(&painter_mtx);
    pnt.begin(plot);
    pnt.setPen(QPen(Qt::black ,1));
        
    std::list<Int_pair>::iterator it = points.begin();
    while(it != points.end()) {
        pnt.drawPoint(it->first, it->second);
        it++;
    }
    pnt.end();
    pthread_mutex_unlock(&painter_mtx);*/
}

//! sets up rasterizer parameters
void Rasterizer::setup(int w, int h)
{
    // setup CORE precision settings
/*    CORE::CORE_init(2);
    CORE::setDefaultPrecision(70, CORE::extLong::getPosInfty());*/
    width = w;
    height = h;
}

//! draws coordinate axis onto the bitmap \c plot
void Rasterizer::draw_axis(const CGAL::Bbox_2& box, QPixmap *plot)
{
    return;
    /*QPainter pnt;
    
    pthread_mutex_lock(&painter_mtx);
    pnt.begin(plot);
    pnt.setPen(QPen(QColor(20,20,200),1,Qt::DashDotDotLine));
        
    double x_scale = (double)width/(box.x_max-box.x_min),
           y_scale = (double)height/(box.y_max-box.y_min);
            
    pnt.moveTo(0,(int)(box.y_max*y_scale));           
    pnt.lineTo(width,(int)(box.y_max*y_scale));
    pnt.moveTo((int)(-box.x_min*x_scale),0);           
    pnt.lineTo((int)(-box.x_min*x_scale),height);
    pnt.end();
    pthread_mutex_unlock(&painter_mtx);*/
}


//      QRgb rgb = colors[index%n_colors].rgb();
//      while(++it != points.end())
//      {
//          int x = (*it).x, y = height - (*it).y;
//          if(x<0||x>=width||y<0||y>=height)
//              continue;
//          uint *ptr = (uint *)plot->scanLine(y)+x;
//          *ptr = rgb;
//      }
