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
// File          : demos/xtriangulate/include/ogl_view.h
// QdX_release   : $Name:  $
// Revision      : $Revision: 1.3 $
// Revision_date : $Date: 2009-07-24 13:21:30 $
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

#ifndef OGL_VIEW_H 
#define OGL_VIEW_H

#include <stddef.h>

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

class OGL_view : public QGLWidget
{
    Q_OBJECT

public:
    OGL_view(QWidget*, const char*);

    virtual ~OGL_view();

    void reset_view();
    void write_STL(std::ostream& os);
   
    bool surf_mesh_computed;  //! base surface mesh computed

    bool lighting;
    bool transparency;
    bool bumpy_surf;
    bool solid_fill;    //! wireframe / solid objects rendering
    bool black_bg_color;
    bool inverse_normals;
        
    bool show_grid;     //! whether to display grid
    bool show_triangles;
    bool show_normals;
    
protected:
    
    void initializeGL();
    void paintGL();
    void resizeGL(int, int);
    void setup_view(); //! sets up viewport & projection
       
    void depth_peeled_render(); // depth-peeled transparency
    void render_scene(); 
    void render_surf();    //! render base surface 

    GLhandleARB load_shader(const char* vs_filename,
        const char* fs_filename);

    void check_extensions();

    void draw_grid();
    void draw_axis();

    void mouseMoveEvent(QMouseEvent* me);
    void mousePressEvent(QMouseEvent* me); 
    void mouseReleaseEvent(QMouseEvent* me);
    void wheelEvent(QWheelEvent * we);

    void world_coord(int x, int y, Point_2d& res);

    void bind_VBO(const VBO_context& ctx);

    void load_VBO(const VBO_context& ctx, const Point_vec_3f* verts, 
        const Point_vec_3f *normals, const std::vector< Color_ub_4 > *colors, 
        const unsigned *indices);

    bool is_extension_supported(const char *);
    
    //!\name data
    //!@{ 
    unsigned int width, height;

    Quaternion q_rotate, q_start; //! current rotation state
    Vector_3d transform;  //! xy-translation and uniform zooming (z component)
    //! translating in z direction (depending on clip planes)
    double z_translate; 
    Point_3d surf_origin;
    double surf_scale;
    
    QPoint start_pt, last_pt;
    ButtonState which_button;

    unsigned framebuf[1], renderbuf[1];
    unsigned n_depth_layers; //! # of passes for depth-peeling
    unsigned *depth_tex, *image_tex; //! handles of depth & image textures
    
    unsigned n_VBOs, *VBOs; //! vertex buffer handles
    unsigned list_base_idx; //! starting index of a range of display lists

    VBO_context mesh_ctx, sil_triangles_ctx, norms_ctx;

    //! shader hanldes
    GLhandleARB lighting_shdr, antialias_shdr, depth_peel_shdr;
    //!@}
};

#endif // OGL_VIEW_H
