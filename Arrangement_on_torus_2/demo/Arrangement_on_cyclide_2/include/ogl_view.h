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
// Library       : QdX
// File          : demos/xsurface/include/ogl_view.h
// QdX_release   : $Name:  $
// Revision      : $Revision$
// Revision_date : $Date$
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

#ifndef OGL_VIEW_H 
#define OGL_VIEW_H

#define GL_GLEXT_PROTOTYPES
#include <qgl.h>
#include <qframe.h>
#include <qpoint.h>

#include "math_routines.h"

typedef boost::array< unsigned char, 4 > Color_ub_4;
typedef boost::array< unsigned int, 3 > Tri_index;

struct VBO_context { // describes a context for rendering from a vertex buffer

    VBO_context() :
        verts_idx(-1u), normals_idx(-1u), colors_idx(-1u),
            /*tex_idx(-1u),*/ index_idx(-1u), index_size(0), min_vidx(0),
            max_vidx(0) {
    }

    VBO_context(unsigned _vidx, unsigned _nidx, unsigned _cidx, 
        /*unsigned _tidx,*/ unsigned _iidx, unsigned _idx_sz, 
        unsigned _min_vidx, unsigned _max_vidx) :
            verts_idx(_vidx), normals_idx(_nidx), colors_idx(_cidx),
            /*tex_idx(_tidx),*/ index_idx(_iidx), index_size(_idx_sz),
            min_vidx(_min_vidx), max_vidx(_max_vidx) {
    } 

    unsigned verts_idx;     // buffer index storing vertices coordinates
    unsigned normals_idx;   // normal coordinates
    unsigned colors_idx;    // colors 
    //unsigned tex_idx;       // texture coordinates
    unsigned index_idx;     // indexed primitives
    unsigned index_size;    // size of indexed array

    // lower and upper limits of vertex indices in index_idx as a hint for
    // glDrawRangeElements
    unsigned min_vidx, max_vidx;
};

class XSurface_view : public QGLWidget
{
    Q_OBJECT

public:

    XSurface_view(QWidget*, const char*);

    virtual ~XSurface_view();

    virtual void reset_view();

    void take_screenshot(const char *fname);
   
    bool surf_mesh_computed;  //! base surface mesh computed
    bool tube_mesh_computed;  //! intersection curves & points mesh computed
    bool aux_mesh_computed;   //! cut circles and pole mesh computed

    bool lighting;
    bool transparency;
    bool bumpy_surf;
    bool solid_fill;    //! wireframe / solid objects rendering
    bool black_bg_color;
    
    bool show_grid;     //! whether to display grid
    bool draw_base_surf; //! whether to draw base surface
    bool draw_curves;
    bool draw_points;
    bool flat_color;    //! draw intersection curves in a flat color
   
    bool draw_outer_circle;
    bool draw_tube_circle;
    bool draw_pole;

protected:
    
    void initializeGL();
    void paintGL();
    void resizeGL(int, int);
    void setup_view(); //! sets up viewport & projection
       
    QImage readFramebuf(int w, int h, bool withAlpha = false);
    void render_scene(); 
    void render_curves();       //! render spatial curves 
    void render_base_surf();    //! render base surface 
    void render_aux();          //! render auxiliary cut/outer cirles & pole

    GLhandleARB load_shader(const char* vs_filename,
        const char* fs_filename);

    void draw_spline_curve();

    void check_extensions();

    void depth_peeled_render(); // depth-peeled transparency
    void draw_grid();
    void draw_axis();

    void mouseMoveEvent(QMouseEvent* me);
    void mousePressEvent(QMouseEvent* me); 
    void mouseReleaseEvent(QMouseEvent* me);
    void wheelEvent(QWheelEvent * we);

    void world_coord(int x, int y, Point_2d& res);

    void construct_plane(const Point_2d& left_top, 
        const Point_2d& extent, unsigned nx, unsigned ny,
            Point_vec_3f& verts, std::vector< Tri_index > *indices,
            std::vector< Color_ub_4 > *colors);

    void draw_tube(const Point_vec_3d& v, double rad, 
        unsigned n_slices, unsigned& base_vidx, unsigned& base_tri_idx,
         Point_vec_3f& verts, Point_vec_3f& normals,  
             std::vector< unsigned >& indices);

    void bind_VBO(const VBO_context& ctx);

    void load_VBO(const VBO_context& ctx, Point_vec_3f* verts, 
        Point_vec_3f *normals, std::vector< Color_ub_4 > *colors, 
        unsigned *indices);

    bool is_extension_supported(const char *);
    
    //!\name data
    //!@{ 
    unsigned int width, height;
    bool render_to_framebuf;

    Quaternion q_rotate, q_start; //! current rotation state
    Vector_3d transform;  //! xy-translation and uniform zooming (z component)
    //! translating in z direction (depending on clip planes)
    double z_translate; 
    Point_3d base_surf_origin;
    double base_surf_scale;
    
    QPoint start_pt, last_pt;
    ButtonState which_button;

    unsigned framebuf[1], renderbuf[1];
    unsigned n_depth_layers; //! # of passes for depth-peeling
    unsigned *depth_tex, *image_tex; //! handles of depth & image textures
    
    unsigned n_VBOs, *VBOs; //! vertex buffer handles
    unsigned list_base_idx; //! starting index of a range of display lists

    VBO_context mesh_ctx, tube_ctx, aux_ctx;
    unsigned tube_vidx_start, tube_iidx_start;

    std::vector< unsigned > tube_idx_counts, tube_idx_shifts;

    //! shader hanldes
    GLhandleARB lighting_shdr, antialias_shdr, depth_peel_shdr;
    //!@}
};

#endif // OGL_VIEW_H
