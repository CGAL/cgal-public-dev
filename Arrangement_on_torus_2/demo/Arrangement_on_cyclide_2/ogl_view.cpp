//============================================================================
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
// File          : demos/xsurface/ogl_view.C
// QdX_release   : $Name:  $
// Revision      : $Revision$
// Revision_date : $Date$
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

#include "include/ogl_view.h"
#include "include/arrangements.h"

#include "include/mainwnd.h"

//#ifdef CGAL_USE_QT

#include <GL/glut.h>

#include <qapplication.h>
#include <qcursor.h>
#include <qimage.h>

//#include <CGAL/algorithm.h>
#include <boost/format.hpp>

#define MESH_VERTS   0
#define MESH_NORMALS 1
#define MESH_COLORS  2
#define MESH_INDICES 3

#define TUBE_VERTS   4
#define TUBE_NORMALS 5
#define TUBE_COLORS  6
#define TUBE_INDICES 7

static const QColor rasterize_colors[] = {
    QColor(139, 0, 0),
    QColor(138, 43, 226),
    QColor(95, 158, 160),
    QColor(0, 0, 139),
    QColor(205, 149, 12),
    QColor(0, 100, 0),
    QColor(139, 0, 139),
    QColor(85, 107, 47),
    QColor(255, 127, 0),
    QColor(0, 206, 209),
    QColor(238, 18, 137),
    QColor(238, 99, 99), 
    QColor(205, 55, 0)
};

static const int n_rast_colors = 13;
extern XSurface_arrangements *engine; 

XSurface_view::XSurface_view(QWidget* parent, const char* name) :
            QGLWidget(parent, name) {

    surf_mesh_computed = false;
    tube_mesh_computed = false;
    aux_mesh_computed = false;

    lighting = true;
    transparency = false;
    bumpy_surf = false;
    solid_fill = true;
    black_bg_color = true;

    show_grid = false;
    draw_base_surf = true;    
    draw_curves = true;
    draw_points = true;   
    flat_color = false;

    draw_outer_circle = false;
    draw_tube_circle = false;
    draw_pole = false;

    render_to_framebuf = false; // rendering to framebuffer (for screenshot)

    base_surf_origin.assign(0);
    base_surf_scale = 1.0;    

    n_depth_layers = 5;
    depth_tex = new unsigned[n_depth_layers * 2];
    image_tex = depth_tex + n_depth_layers;

    n_VBOs = 12;
    VBOs = new unsigned[n_VBOs]; 

    z_translate = 0;
    width = parent->width();
    height = parent->height();
    reset_view();
}

unsigned texID[3]; QImage tex;

void XSurface_view::initializeGL() {

    char *vendor = (char*)glGetString(GL_VENDOR);
    char *version = (char*)glGetString(GL_VERSION);
    char *renderer = (char*)glGetString(GL_RENDERER);

    std::cout << "Vendor: " << vendor << "\nVersion: " << version <<
        "\nRenderer: " << renderer << std::endl;

    int nidxs, nverts; 
    glGetIntegerv(GL_MAX_ELEMENTS_INDICES, &nidxs);
    glGetIntegerv(GL_MAX_ELEMENTS_VERTICES, &nverts);
    std::cout << "Max # of indexes per VBO: " << nidxs << 
            "; max vertices: " << nverts << "\n";
        
    qglClearColor(black);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_RESCALE_NORMAL);
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnable(GL_VERTEX_PROGRAM_TWO_SIDE);

    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    glFrontFace(GL_CCW);

    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1); 
    glEnable(GL_COLOR_MATERIAL);

    float ambient[] = {0.1f, 0.1f, 0.1f, 0.0f};
    glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
   
    GLfloat pos[4] = {0.7f, 0, -1.0f, 1.0f};
    glLightfv(GL_LIGHT0, GL_POSITION, pos);
    //glEnable(GL_LIGHT0);

    glShadeModel(GL_SMOOTH);

    check_extensions();
    
    glGenBuffers(n_VBOs, VBOs);
    list_base_idx = glGenLists(4);

    glGenFramebuffersEXT(1, framebuf);
    glGenTextures(n_depth_layers, depth_tex);
    glGenTextures(n_depth_layers, image_tex);

    lighting_shdr = load_shader("sample.vert", "sample.frag");
    antialias_shdr = load_shader("antialias.vert", "antialias.frag");

    //glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
    //glEnable(GL_MULTISAMPLE);

    //glEnable(GL_POLYGON_OFFSET_FILL);
    //glPolygonOffset(1.0f, 1.0f);
}

void XSurface_view::reset_view() {
   
    transform = Vector_3d(0.0, 0.0, 1.0);
    quaternion_from_axis(q_rotate, Vector_3d(1, 0, 0), 45.0);
}

void XSurface_view::setup_view() {

    glViewport(0, 0, width, height);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    
    double ratio = (double)width / height;
    //glOrtho(-ratio, ratio, -1, 1, 0.01, 10);
    gluPerspective(60, ratio, 0.01, 40.0);
    z_translate = -4.0;    // translate according to near/far clip planes

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
   
    for(unsigned i = 0; i < n_depth_layers; i++) {

        glBindTexture(GL_TEXTURE_RECTANGLE_ARB, depth_tex[i]);
        glTexImage2D(GL_TEXTURE_RECTANGLE_ARB, 0, GL_DEPTH_COMPONENT, 
            width, height, 0, GL_DEPTH_COMPONENT, GL_UNSIGNED_INT, 0L);

        glTexParameterf(GL_TEXTURE_RECTANGLE_ARB,
                        GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameterf(GL_TEXTURE_RECTANGLE_ARB,
                        GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameterf(GL_TEXTURE_RECTANGLE_ARB,
                        GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameterf(GL_TEXTURE_RECTANGLE_ARB,
                        GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameterf(GL_TEXTURE_RECTANGLE_ARB,
                        GL_TEXTURE_COMPARE_MODE_ARB,
                        GL_COMPARE_R_TO_TEXTURE_ARB); // Crucial!
        glTexParameterf(GL_TEXTURE_RECTANGLE_ARB,
                        GL_TEXTURE_COMPARE_FUNC_ARB, GL_LEQUAL);
    
        glBindTexture(GL_TEXTURE_RECTANGLE_ARB, image_tex[i]);
        glTexImage2D(GL_TEXTURE_RECTANGLE_ARB, 0, GL_RGBA, width, height, 0,
                 GL_RGBA, GL_UNSIGNED_BYTE, 0L);
        glTexParameteri(GL_TEXTURE_RECTANGLE_ARB,
                        GL_TEXTURE_MIN_FILTER,  GL_NEAREST);
        glTexParameteri(GL_TEXTURE_RECTANGLE_ARB,
                        GL_TEXTURE_MAG_FILTER,  GL_NEAREST);
        glTexParameteri(GL_TEXTURE_RECTANGLE_ARB,
                        GL_TEXTURE_WRAP_S,  GL_CLAMP);
        glTexParameteri(GL_TEXTURE_RECTANGLE_ARB,
                        GL_TEXTURE_WRAP_T,  GL_CLAMP);
    }
}

void XSurface_view::paintGL() {

    if(transparency) {
        depth_peeled_render();
    } else {
        glEnable(GL_DEPTH_TEST);
        render_scene();
    }
    draw_axis();
}

void XSurface_view::take_screenshot(const char *fname) {

    unsigned oldw = width, oldh = height;
    width = 2500, height = 2500;
    
    setup_view();

    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, framebuf[0]);
    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT,
        GL_TEXTURE_RECTANGLE_ARB, image_tex[0], 0);
    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT,
        GL_TEXTURE_RECTANGLE_ARB, depth_tex[0], 0);

    GLenum status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
    if(status != GL_FRAMEBUFFER_COMPLETE_EXT) {
        std::cerr << "incomplete framebuf ---------\n";
    }

    render_to_framebuf = true;
    paintGL();
    render_to_framebuf = false;

    QImage img = readFramebuf(width, height, true);
    img.save(QString("images/") + fname, "PNG");

    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
    
    width = oldw, height = oldh;
    setup_view();
}

void XSurface_view::depth_peeled_render() {

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);

    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, framebuf[0]);
    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT,
        GL_TEXTURE_RECTANGLE_ARB, image_tex[0], 0);
    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT,
        GL_TEXTURE_RECTANGLE_ARB, depth_tex[0], 0);

    GLenum status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
    if(status != GL_FRAMEBUFFER_COMPLETE_EXT) {
        std::cerr << "incomplete framebuf\n";
    }
    
    transparency = false; // draw the first layer (everything)
    render_scene();
    glFlush();

    transparency = true;
    
    glUseProgramObjectARB(lighting_shdr);
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, framebuf[0]);
    glEnable(GL_TEXTURE_RECTANGLE_ARB);

    // Now do additional passes for each other layer with the pixel shader
    // attached and the previous layer's depth map bound as input through
    // "uniform sampler2DRectShadow DT" in the shader.
    for(unsigned i = 1; i < n_depth_layers; i++) {

        glBindTexture(GL_TEXTURE_RECTANGLE_ARB, depth_tex[i-1]);
        
        glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT,
              GL_TEXTURE_RECTANGLE_ARB, image_tex[i], 0);
        glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT,
              GL_TEXTURE_RECTANGLE_ARB, depth_tex[i], 0);

        render_scene();
        glFlush();
    }

    if(!render_to_framebuf)
        glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
    else
        glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, framebuf[0]);
    glUseProgramObjectARB(0);

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(-width/2, width/2, height/2, -height/2, -10, 200);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    glActiveTexture(GL_TEXTURE0);
    glEnable(GL_TEXTURE_RECTANGLE_ARB);
    glColor4f(1.0, 1.0, 1.0, 1.0);
    
    glEnable(GL_BLEND);
    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
    
    glDisable(GL_DEPTH_TEST);
    //glUseProgramObjectARB(antialias_shdr);
//     glCullFace(GL_BACK);
//     glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    // Draw each layer over each other in reverse order
    for(unsigned i = n_depth_layers - 1; (int)i >= 0; i--) {
        
        glBindTexture(GL_TEXTURE_RECTANGLE_ARB, image_tex[i]);
        glBegin(GL_QUADS);
        glTexCoord2f(0.0, 0.0); glVertex3f(-width/2.0, height/2.0, 0.0f);
        glTexCoord2f(width, 0.0); glVertex3f(width/2.0, height/2.0, 0.0f);
        glTexCoord2f(width, height); glVertex3f(width/2.0, -height/2.0, 0.0f);
        glTexCoord2f(0.0, height);  glVertex3f(-width/2.0, -height/2.0, 0.0f);
        glEnd();
    }
    glDisable(GL_TEXTURE_RECTANGLE_ARB);
    glDisable(GL_BLEND);

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
}

void XSurface_view::render_scene() {

    if(black_bg_color)
        glClearColor(0.0, 0.0, 0.0, 1.0);
    else
        glClearColor(1.0, 1.0, 1.0, 1.0);
    
    glClearDepth(1.0);
    
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    //std::cerr << "base_surf_origin = " << base_surf_origin << "\n";
    glTranslated(transform.x()*fabs(z_translate - base_surf_origin[0]), 
        transform.y()*fabs(z_translate) - base_surf_origin[1], 
        z_translate - base_surf_origin[2]);

    //std::cerr << "transform = " << transform << "\n";

    glScaled(transform.z(), transform.z(), transform.z());
               
    double mat[16];
    quaternion_2_matrix(q_rotate, mat);
    glMultMatrixd(mat);
    draw_grid();
    
    glPolygonMode(GL_FRONT_AND_BACK, (solid_fill ? GL_FILL : GL_LINE));
    
    if(lighting && !transparency)
        glUseProgramObjectARB(lighting_shdr);
    
    GLint alpha_loc = glGetUniformLocation(lighting_shdr, "enable_alpha"),
        bumpy_loc = glGetUniformLocation(lighting_shdr, "bumpy_surf"),
        depth_peel_loc = glGetUniformLocation(lighting_shdr,
             "do_depth_peeling");
        
    glUniform1i(depth_peel_loc, transparency);

    glUniform1i(bumpy_loc, false);
    glUniform1i(alpha_loc, false);
    render_aux();
    render_curves();
    glUniform1i(bumpy_loc, bumpy_surf);
    glUniform1i(alpha_loc, true);

//     if(transparency) {
//         glBlendEquation(GL_FUNC_ADD);
//         glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
//       //glBlendFunc(GL_ONE_MINUS_DST_ALPHA, GL_DST_ALPHA);
//         
//         glEnable(GL_BLEND);
//         glCullFace(GL_FRONT);
//         glUniform1f(rescale_loc, -1);       
//         //glColor4ub(250, 40, 250, 150);
//         
//         glColor4ub(200, 100, 22, 150);
//         //glColor4ub(115, 232, 255, 150);
//         render_base_surf();
//     }
//   
//     glCullFace(GL_BACK);
    //if(transparency);
        glDisable(GL_CULL_FACE);  // otherwise the inner surface is not drawn
    
#define GOLD 255, 215, 0
#define EMERALD 0x50, 0xc8, 0x78
#define NICE_BLUE 40, 80, 200
#define VIOLET  200, 10, 80
#define ORANGE 253, 192, 24      
#define GREEN 50, 200, 100
#define BLUE 40, 100, 230

    // NOTE: alpha component is set through the shader
    glColor3ub(VIOLET);      // front-face color
  
    glSecondaryColor3ub(EMERALD);  
    render_base_surf();

    //if(transparency)
      //  glEnable(GL_CULL_FACE);  
}

#if 1
static void draw_int(const Point_3d& pos, int n) {

    static char buf[64];
    int len = sprintf(buf, "%d", n);
    glColor3ub(255, 255, 255);
    glRasterPos3d(pos[0], pos[1], pos[2]);
    for(int i = 0; i < len; i++) {
        glutBitmapCharacter(GLUT_BITMAP_9_BY_15, buf[i]);
    }
    
}
#endif

void XSurface_view::render_curves() {
       
    if(!draw_curves && !draw_points)
        return;

    if(engine->get_arcs_approx().size() > 0 && !tube_mesh_computed) {

        Point_vec_3f verts, normals; 
        std::vector< unsigned > idxs;
        std::vector< Color_ub_4 > colors;

        int i;
        const std::vector<Point_vec_3d>& xarcs = engine->get_arcs_approx();
        const std::vector< unsigned >& cindex = engine->get_color_index(),
            cpindex = engine->get_points_color_index();
        
        std::vector<Point_vec_3d>::const_iterator vit;
        std::vector< unsigned >::const_iterator cit;

        tube_idx_counts.resize(xarcs.size());
        tube_idx_shifts.resize(xarcs.size());
        
        glNewList(list_base_idx+3, GL_COMPILE);

        double rad = engine->get_base_surf_extent()*3 / 200;
        unsigned base_vidx = 0, base_iidx = 0;
        for(vit = xarcs.begin(), cit = cindex.begin(), i = 0; 
            vit != xarcs.end(); vit++, cit++, i++) {

           //if(i != 119 && i != 48) continue; // 53

            const QColor& cc = rasterize_colors[*cit % n_rast_colors];

#if 0            
            glColor3ub(cc.red(), cc.green(), cc.blue()); 
            const Point_vec_3d& tmp = *vit;
            Point_vec_3d::const_iterator pit;
            int r = 10, g = 0, b = 0;
            glPointSize(4);
            glLineWidth(3);

            unsigned sz = tmp.size();
            Point_3d where;
            glBegin(GL_LINE_STRIP);  // GL_LINE_STRIP GL_POINTS
            for(pit = tmp.begin(); pit != tmp.end(); pit++, sz--) {
                glColor3ub(255, g, b); g+=15, b+=15;                
                glVertex3d((*pit)[0], (*pit)[1], (*pit)[2]);
                
                if(sz == tmp.size()/2) 
                    where = *pit;
            }
            glEnd();
            draw_int(where, i);
#endif             
            tube_idx_shifts[i] = base_iidx * sizeof(unsigned);
            unsigned vstart = base_vidx, istart = base_iidx;
            draw_tube(*vit, rad, 10, base_vidx, base_iidx, verts, normals,
                 idxs);

            tube_idx_counts[i] = base_iidx - istart;
            colors.resize(verts.size());
            for(; vstart < base_vidx; vstart++) {
                colors[vstart][0] = cc.red(),
                colors[vstart][1] = cc.green(), 
                colors[vstart][2] = cc.blue(),
                colors[vstart][3] = 1;
            }
        }
        glEndList();

          // _vidx,  _nidx,  _cidx, _iidx,  _idx_sz)
        tube_ctx = VBO_context(4, 5, 6, 7, idxs.size(), 0, verts.size()-1);
        load_VBO(tube_ctx, &verts, &normals, &colors, idxs.data());
        
        std::cerr << "\ntube mesh computed\n# of vertices = " << verts.size()
             << "; # of normals = " << normals.size() << " # of indices = " <<
                 idxs.size() << "\n";

//         glColor4ub(100, 100, 100, 40);
//         glBegin(GL_QUADS);
//         std::vector< CGAL::Bbox_2 >::const_iterator bbix = boxxes.begin();
//         while(bbix != boxxes.end()) {
//             glVertex3d(bbix->xmin(), bbix->ymin(),0);
//             glVertex3d(bbix->xmax(), bbix->ymin(),0);
//             glVertex3d(bbix->xmax(), bbix->ymax(),0);
//             glVertex3d(bbix->xmin(), bbix->ymax(),0);
//             bbix++;
//         }
//         glEnd();    
       
        rad = engine->get_base_surf_extent()*5 / 200;
        glNewList(list_base_idx, GL_COMPILE);
        glutSolidSphere(rad, 15, 15);
        glEndList();

        glNewList(list_base_idx+1, GL_COMPILE);
        glColor3ub(250, 40, 40);
        //glBegin(GL_POINTS);
               
        const Point_vec_3d& pts = engine->get_points_approx(); 
        Point_vec_3d::const_iterator pit;
        for(pit = pts.begin(), cit = cpindex.begin(); pit != pts.end();
                 pit++) {
            
            if(!cpindex.empty()) {
                const QColor& cc = rasterize_colors[(*cit++) % n_rast_colors];
                glColor3ub(cc.red()+30, cc.green()+30, cc.blue()+30); 
            }
            //glPushMatrix(); 
            glTranslatef((float)(*pit)[0], (float)(*pit)[1], (float)(*pit)[2]);
            glCallList(list_base_idx);
            glTranslatef((float)-(*pit)[0], (float)-(*pit)[1],
                 (float)-(*pit)[2]);
            //glutSolidSphere(rad, 15, 15);
            //glPopMatrix();
        }
        glEndList();
        tube_mesh_computed = true;
    }
        
    if(tube_mesh_computed) {
        if(draw_curves) {
            unsigned save = tube_ctx.colors_idx;
            if(flat_color) {
                glColor3ub(120, 250, 40);
                tube_ctx.colors_idx = -1;
            }

#if 1
            bind_VBO(tube_ctx);
            glMultiDrawElements(GL_TRIANGLE_STRIP, 
                (int *) tube_idx_counts.data(), GL_UNSIGNED_INT, 
                    (const void **)tube_idx_shifts.data(),
                        tube_idx_counts.size());
            glPopClientAttrib();
#else            
            glCallList(list_base_idx+3);
#endif
            tube_ctx.colors_idx = save;
        }

        if(draw_points) 
            glCallList(list_base_idx+1);
    }
}

void XSurface_view::render_base_surf() {

    if(engine->base_surface_valid() && !surf_mesh_computed) {
        
        Point_2d lt, extent;
        lt[0] = -35, lt[1] = -35;
        extent[0] = 70, extent[1] = 70;

        Point_vec_3f verts, normals; 
        std::vector< Tri_index > indices;

        construct_plane(lt, extent, 100, 100, verts, &indices, 0);
        //CGAL::output_range(std::cerr, verts.begin(), verts.end());
        engine->compute_spatial_coords(verts, verts, &normals);
        
         // _vidx,  _nidx,  _cidx,  _tidx, _iidx,  _idx_sz
        mesh_ctx = VBO_context(0, 1, -1, 3, indices.size()*3,
            0, verts.size()-1);
        load_VBO(mesh_ctx, &verts, &normals, 0, (unsigned *)indices.data());

        engine->get_base_surf_origin(base_surf_origin);
        base_surf_scale = 2.0 / engine->get_base_surf_extent();
        transform[2] = base_surf_scale;

        surf_mesh_computed = true;
    }
    
    if(surf_mesh_computed && draw_base_surf) {
        bind_VBO(mesh_ctx);
        glDrawRangeElements(GL_TRIANGLES, mesh_ctx.min_vidx, mesh_ctx.max_vidx,
            mesh_ctx.index_size, GL_UNSIGNED_INT, 0);
        
        glPopClientAttrib();
    }
}

void XSurface_view::render_aux() {
    
    if(!draw_tube_circle && !draw_outer_circle && !draw_pole)
        return;
   
    if(engine->base_surface_valid() && !aux_mesh_computed) {
        
        Point_vec_3f vtmp;
        Point_2d lt, extent;
        lt[0] = -35, lt[1] = 0;
        extent[0] = 70, extent[1] = 0;

        construct_plane(lt, extent, 100, 0, vtmp, 0, 0);
        engine->compute_aux_coords(vtmp, vtmp, 0); // outer circle
        
        Point_vec_3d v3d;
        v3d.resize(vtmp.size());
        for(unsigned i = 0; i < vtmp.size(); i++) {
            v3d[i][0] = vtmp[i][0], v3d[i][1] = vtmp[i][1], 
            v3d[i][2] = vtmp[i][2];
        }            

        Point_vec_3f verts, normals;
        std::vector< unsigned > idxs; 
        unsigned base_vidx = 0, base_iidx = 0;  
        double ref = engine->get_base_surf_extent() / 35;

        draw_tube(v3d, ref, 20, base_vidx, base_iidx, verts, normals, idxs);

        vtmp.clear(), v3d.clear();
        lt[0] = 0, lt[1] = -35;
        extent[0] = 0, extent[1] = 70;
        construct_plane(lt, extent, 0, 100, vtmp, 0, 0);
        engine->compute_aux_coords(vtmp, vtmp, 1); // tube
        
        v3d.resize(vtmp.size());
        for(unsigned i = 0; i < vtmp.size(); i++) {
            v3d[i][0] = vtmp[i][0], v3d[i][1] = vtmp[i][1], 
            v3d[i][2] = vtmp[i][2];
        }
        
        tube_vidx_start = base_vidx;
        tube_iidx_start = base_iidx;
        //std::cerr << "base_vidx = " << base_vidx << "\n";
        draw_tube(v3d, ref, 20, base_vidx, base_iidx, verts, normals, idxs);

         // _vidx,  _nidx,  _cidx, _iidx,  _idx_sz)
        aux_ctx = VBO_context(8, 9, -1, 11, idxs.size(), 0, verts.size()-1);
        load_VBO(aux_ctx, &verts, &normals, 0, idxs.data());

        double rad = engine->get_base_surf_extent() / 15;

        Point_3f pole;
        engine->get_pole_coords(pole);

        glNewList(list_base_idx+2, GL_COMPILE);
        glTranslatef(pole[0], pole[1], pole[2]);
        glutSolidSphere(rad, 15, 15);
        glTranslatef(-pole[0], -pole[1], -pole[2]);
        glEndList();

        aux_mesh_computed = true;
    }
    
    if(aux_mesh_computed) {
        if(draw_tube_circle || draw_outer_circle)
            bind_VBO(aux_ctx);

        if(draw_outer_circle) {
            glColor3ub(120, 120, 120);
            glDrawRangeElements(GL_TRIANGLE_STRIP, 0, tube_vidx_start-1,
                tube_iidx_start, GL_UNSIGNED_INT, 0);
        }
        if(draw_tube_circle) {
            glColor3ub(150, 80, 80);
            glDrawRangeElements(GL_TRIANGLE_STRIP, tube_vidx_start,
                 aux_ctx.max_vidx, aux_ctx.index_size - tube_iidx_start,
                 GL_UNSIGNED_INT, (void *)(tube_iidx_start *
                        sizeof(unsigned)));            
        }

        if(draw_tube_circle || draw_outer_circle)
            glPopClientAttrib();

        if(draw_pole) {
            glColor3ub(120, 50, 120);
            glCallList(list_base_idx+2);
        }
    }
}

void XSurface_view::draw_tube(const Point_vec_3d& v, double rad, 
        unsigned n_slices, unsigned& base_vidx, unsigned& base_iidx,
        Point_vec_3f& verts, Point_vec_3f& normals,  
             std::vector< unsigned >& indices) {

    if(v.size() < 2u)    
        return;

    verts.resize(base_vidx + v.size()*n_slices);
    normals.resize(base_vidx + v.size()*n_slices);
    //indices.resize(base_tri_idx + (v.size()-1)*n_slices*2);
    indices.resize(base_iidx + (v.size() - 1)*(n_slices+1)*2);
    
    unsigned vidx = base_vidx, iidx = base_iidx;

    Point_vec_3d::const_iterator vit = v.begin();
    Vector_3d cur((*vit)[0], (*vit)[1], (*vit)[2]),
        next((*(vit+1))[0], (*(vit+1))[1], (*(vit+1))[2]);

    Vector_3d v0, v1 = next - cur, up(0, -v1.z(), v1.y());
    Vector_3d dir = v1 * 1;
    if(fabs(v1.y()) <= 1e-12 && fabs(v1.z()) <= 1e-12) {
        up[0] = -v1.y();
        up[1] = v1.x();
        up[2] = 0;
    }
    up.normalize();
    Vector_3d n = dir.cross_product(up);
    n.normalize();

    while(1) {

        double angle = 0, step = 2*M_PI/n_slices;
        for(unsigned i = 0; i < n_slices; i++, angle += step, vidx++) {
            
            Vector_3d t = n * cos(angle) + up * sin(angle);
            t.normalize();
            normals[vidx][0] = t[0], normals[vidx][1] = t[1],
            normals[vidx][2] = t[2];
            t = cur + t * rad;
            verts[vidx][0] = t[0], verts[vidx][1] = t[1],
            verts[vidx][2] = t[2];
            if(vit + 1 == v.end())
                continue;

            indices[iidx] = vidx;
            indices[iidx+1] = vidx+n_slices;
            iidx += 2;            
        }
        if(vit + 1 != v.end()) {
            indices[iidx] = vidx - n_slices;
            indices[iidx+1] = vidx;
            iidx += 2;
        }

        if(++vit == v.end())
            break;

        //std::cerr << "up = " << up << "; pos = " << *vit << "\n";
        v0 = next - cur;
        cur = next * 1.0; // cur.copy_on_write
       
        Vector_3d side, n2;
        bool collinear = true;
        if(vit + 1 != v.end()) {

            //next.copy_on_write();
            next = Vector_3d((*(vit+1))[0], (*(vit+1))[1], (*(vit+1))[2]);
            v1 = next - cur;
            Vector_3d diff = v1 - v0;
            //std::cerr << "diff = " << diff << "\n";

            if(fabs(diff.x()) <= 1e-12 && fabs(diff.y()) <= 1e-12 && 
                    fabs(diff.z()) <= 1e-12) {
                side = v1*1;
            } else {
                v1.normalize();
                v0.normalize();
//                 std::cerr << "prev = " << *(vit-1) << "; cur = " <<
//                     cur << "; next = " << next << "\n";
                Vector_3d mid = (v1 - v0);
                n2 = v0.cross_product(v1);
                side = mid.cross_product(n2);
                collinear = false;
            }
        } else 
            side = v0;
                
        up = n.cross_product(side);
        up.normalize();
        //std::cerr << "side = " << side << "; n = " << n << "; up = " << up << "\n\n";

        if(!collinear) {
            n = v1.cross_product(up);
            n.normalize();  
        }
    }   
    base_vidx = vidx;
    base_iidx = iidx;
}

void XSurface_view::construct_plane(const Point_2d& left_top, 
        const Point_2d& extent, unsigned nx, unsigned ny,
        Point_vec_3f& verts, std::vector< Tri_index > *indices,
            std::vector< Color_ub_4 > *colors) {
    
    double step_x = extent[0] / nx, step_y = extent[1] / ny, x, y;
    verts.resize((nx + 1)*(ny + 1));

    std::vector< Tri_index >::iterator iit;
    if(indices != 0) {
        indices->resize(nx*ny*2); // triangle indices   
        iit = indices->begin();
    }
    if(colors != 0)
        colors->resize(verts.size());
    //tex_coords.resize(verts.size());

    double ww = 1+extent[0]/2, hh = 1+extent[1]/2;
    unsigned i, j, vidx;
    for(j = 0, vidx = 0, y = left_top[1]; j <= ny; j++, y += step_y) {
        
        //unsigned tri_idx = j * nx * 2;
        for(i = 0, x = left_top[0]; i <= nx; i++, vidx++, x += step_x) { 
            
            double invx, invy;
            invx = (x < 0 ? -ww/(-x+1)+1 : ww/(x+1)-1);
            invy = (y < 0 ? -hh/(-y+1)+1 : hh/(y+1)-1);
            // 0 is mapped to w/2 (but should be both to -w/2 & w/2)
            verts[vidx][0] = invx, verts[vidx][1] = invy, 
            verts[vidx][2] = 0.0;
                   
//             double r = log(fabs(x)/extent[0]);
//             double g = log(fabs(y)/extent[1]);
//             double b = log(fabs(x+y)/extent[0]);
//             colors[vidx] = Color_ub_4((unsigned char)(r * 128.0),
//                (unsigned char)(g * 128.0), (unsigned char)(b * 128.0), 1u);

            if(indices != 0 && i < nx && j < ny) {
                unsigned t = vidx + nx + 1;
                (*iit)[0] = vidx, (*iit)[1] = vidx + 1, (*iit)[2] = t;
                iit++;
                (*iit)[0] = vidx+1, (*iit)[1] = t + 1, (*iit)[2] = t;
                iit++;
            }
        }
    }
}

void XSurface_view::draw_spline_curve() {

    /*const std::vector<Point_vec_4f>& xarcs = engine->get_arcs_approx();
    
    Point_4f pts[] = {
        Point_4f(0.0, 0.0, 0, 1),
        Point_4f(0.05, 0.1, 0.1, 1),
        Point_4f(0.2, 0.2, 0.3, 1),
        Point_4f(0.4, 0.0, 0.2, 1), 
        Point_4f(0.5, 0.1, 0.1, 1)};

    float knots[] = { 0, 0, 0.4, 0.4, 0.6, 0.6, 1, 1, 1};

    const int npts = 4;

    glPointSize(4);
    glColor3ub(0, 255, 0);
    glBegin(GL_POINTS);
    for(int i = 0; i < npts; i++) {
        glVertex4fv((float *)&pts[i]);
    }
    glEnd();

    GLUnurbsObj *rend = gluNewNurbsRenderer();

    glDisable(GL_LIGHTING);
    glColor3ub(255, 255, 255);
    gluBeginCurve(rend);
    gluNurbsCurve(rend, npts + 4, knots, 4, (float *)pts, 4, GL_MAP1_VERTEX_4);
    gluEndCurve(rend);

    gluDeleteNurbsRenderer(rend);*/
}

//! renders \c count primitives from \c ctx 
void XSurface_view::bind_VBO(const VBO_context& ctx) {
    //  GLenum primitive_type, unsigned count, unsigned shift) {
   
    glPushClientAttrib(GL_CLIENT_VERTEX_ARRAY_BIT);

    glBindBuffer(GL_ARRAY_BUFFER, VBOs[ctx.verts_idx]);
    glVertexPointer(3, GL_FLOAT, 0, 0);

    if(ctx.normals_idx != -1u) {
        glEnableClientState(GL_NORMAL_ARRAY);
        glBindBuffer(GL_ARRAY_BUFFER, VBOs[ctx.normals_idx]);
        glNormalPointer(GL_FLOAT, 0, 0);
    }

    if(ctx.colors_idx != -1u) {
        glEnableClientState(GL_COLOR_ARRAY);
        glBindBuffer(GL_ARRAY_BUFFER, VBOs[ctx.colors_idx]);
        glColorPointer(4, GL_UNSIGNED_BYTE, 0, 0);
    }
        
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, VBOs[ctx.index_idx]);
}

void XSurface_view::load_VBO(const VBO_context& ctx,
        Point_vec_3f* verts, Point_vec_3f *normals, 
        std::vector< Color_ub_4 > *colors, unsigned *indices) {
    
    glBindBuffer(GL_ARRAY_BUFFER, VBOs[ctx.verts_idx]);
    glBufferData(GL_ARRAY_BUFFER, verts->size()*3*sizeof(float),
            verts->data(), GL_STATIC_DRAW);

    if(ctx.normals_idx != -1u && normals != 0) {
        glBindBuffer(GL_ARRAY_BUFFER, VBOs[ctx.normals_idx]);
        glBufferData(GL_ARRAY_BUFFER, normals->size()*3*sizeof(float),
            normals->data(), GL_STATIC_DRAW);
    }

    if(ctx.colors_idx != -1u && colors != 0) {
        glBindBuffer(GL_ARRAY_BUFFER, VBOs[ctx.colors_idx]);
        glBufferData(GL_ARRAY_BUFFER, colors->size()*4*sizeof(unsigned char),
            colors->data(), GL_STATIC_DRAW);
    }

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, VBOs[ctx.index_idx]);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER,
        ctx.index_size*sizeof(unsigned), indices, GL_STATIC_DRAW);
}

// ripped from OGL::grabFrameBuffer()
QImage XSurface_view::readFramebuf(int w, int h, bool withAlpha) {
    
    makeCurrent();
    QImage res;
    if ( format().rgba() ) {
    res = QImage( w, h, 32 );
    glReadPixels( 0, 0, w, h, GL_RGBA, GL_UNSIGNED_BYTE, res.bits() );
    if ( QImage::systemByteOrder() == QImage::BigEndian ) {
        // OpenGL gives RGBA; Qt wants ARGB
        uint *p = (uint*)res.bits();
        uint *end = p + w*h;
        if ( withAlpha && format().alpha() ) {
        while ( p < end ) {
            uint a = *p << 24;
            *p = (*p >> 8) | a;
            p++;
        }
        }
        else {
        while ( p < end )
            *p++ >>= 8;
        }
    }
    else {
        // OpenGL gives ABGR (i.e. RGBA backwards); Qt wants ARGB
        res = res.swapRGB();
    }
    res.setAlphaBuffer( withAlpha && format().alpha() );
    }
    else {

    }

    return res.mirror();
}

void XSurface_view::world_coord(int x, int y, Point_2d& res) {

    res[0] = (2.0*x - width) / width, 
    res[1] = (height - 2.0*y) / height;
}

void XSurface_view::mousePressEvent(QMouseEvent* me) {
    start_pt = last_pt = me->pos();
    which_button = me->button();
    q_start = q_rotate;
}

void XSurface_view::mouseMoveEvent(QMouseEvent* me) {
    
    const QPoint& cur_pt = me->pos(), first = 
        (which_button == LeftButton ? start_pt : last_pt); 
    Point_2d pt1, pt2;    

    world_coord(first.x(), first.y(), pt1);
    world_coord(cur_pt.x(), cur_pt.y(), pt2);

    switch(which_button) {
    
    case MidButton: {
        double dir = (pt2[1] - pt1[1] > 0 ? 0.97 : 1.03);
        transform[2] *= dir;
        break;
    }
    case RightButton:
        setCursor(QCursor(Qt::SizeAllCursor));
        transform[0] += (pt2[0] - pt1[0]);
        transform[1] += (pt2[1] - pt1[1]);
        break; 

    case LeftButton: {
        setCursor(QCursor(Qt::SizeAllCursor));
        Quaternion tmp;
        quaternion_trackball(tmp, pt1, pt2);
        quaternion_multiply(q_rotate, q_start, tmp);
        quaternion_normalize(q_rotate);        
        break;
    }
    default:;
    }
    last_pt = cur_pt;
    updateGL();
}

void XSurface_view::mouseReleaseEvent(QMouseEvent* me) {
    setCursor(QCursor(Qt::ArrowCursor));
}

void XSurface_view::wheelEvent(QWheelEvent *we) {
    
    double dir = (we->delta() > 0 ? 0.95 : 1.05);
    transform[2] *= dir;
    //std::cerr << transform[2] << "\n";
    updateGL();
}

void XSurface_view::draw_grid() {

    if(!show_grid)
        return;    

    int n_cells = 10, i;
    float ext = engine->get_base_surf_extent(),
        grid_size = (ext == 0 ? 1.5f : ext*0.75f), 
        step = grid_size / n_cells, c;

    glLineWidth(1);
    glBegin(GL_LINES);
    // z-axis light blue
    glColor3f(0.3f, 0.3f, 0.9f);
    glVertex3f(0.0, 0.0, -grid_size);
    glVertex3f(0.0, 0.0, grid_size);

    for(c = -grid_size, i = 0; i <= 2*n_cells; c += step, i++) {
        if(i == n_cells) // // x-axis light red
            glColor3f(0.9f, 0.3f, 0.3f);
        else
            glColor3f(0.5f, 0.5f, 0.5f);
        glVertex2f(-grid_size, c);
        glVertex2f(grid_size, c);

        if(i == n_cells) // y-axis light green
            glColor3f(0.3f, 0.9f, 0.3f);
        glVertex2f(c, -grid_size);
        glVertex2f(c, grid_size);
    }
    glEnd();
}

void XSurface_view::draw_axis() {

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();

    glLoadIdentity();
    glOrtho(-1, 1, -1, 1, -1, 1);
    
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glTranslatef(-0.75f, -0.75f, 0);
    glScalef(0.2f, 0.2f * (float)width/height, 0.2f);
    
    double mat[16];
    quaternion_2_matrix(q_rotate, mat);
    glMultMatrixd(mat);
    glScalef(1, -1, -1); 

    glClear(GL_DEPTH_BUFFER_BIT);
 
    glLineWidth(1);
    if(black_bg_color)
        glColor3f(1.0f, 1.0f, 1.0f);
    else
        glColor3f(0.0f, 0.0f, 0.0f);
   
    // x-axis
    glPushMatrix();
    glBegin(GL_LINES);
    glVertex3f(-0.2, 0.0, 0.0);
    glVertex3f(0.9, 0.0, 0.0);
    glEnd();
    glTranslatef(0.9,0.0,0.0);
    glRotatef(90.0,0.0,1.0,0.0);
    glutSolidCone(0.1,0.1,10,1);
    glPopMatrix();
    
    // y-axis
    int d=1;
    glPushMatrix();
    glBegin(GL_LINES);
    glVertex3f( 0.0, d*0.2, 0.0);
    glVertex3f( 0.0, -d*0.9, 0.0);
    glEnd();
    glTranslatef(0.0,-d*0.9,0.0);
    glRotatef(d*90.0,1.0,0.0,0.0);
    glutSolidCone(0.1,0.1,10,1);
    glPopMatrix();

    // z-axis
    glPushMatrix();
    glBegin(GL_LINES);
    glVertex3f( 0.0, 0.0, 0.2);
    glVertex3f( 0.0, 0.0, -0.9);
    glEnd();
    glTranslatef(0.0,0.0,-0.9);
    glRotatef(180,1.0,0.0,0.0);
    glutSolidCone(0.1,0.1,10,1);
    glPopMatrix();
    
    glLineWidth(1);
    glDisable(GL_LIGHTING);
    glRasterPos3f(1.0f, 0, 0);
    glutBitmapCharacter( GLUT_BITMAP_9_BY_15, 'x' );
    glRasterPos3f(0, -d*1.0f, 0);
    glutBitmapCharacter( GLUT_BITMAP_9_BY_15, 'y' );
    glRasterPos3f(0, 0, -1.0f);
    glutBitmapCharacter( GLUT_BITMAP_9_BY_15, 'z' );
  
    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
}

void XSurface_view::resizeGL(int w, int h) {
    width = w;
    height = h;
    setup_view();
    updateGL();
}

void XSurface_view::check_extensions() {
    
    if(!is_extension_supported("GL_ARB_vertex_buffer_object")) {
        std::cerr << "No VBO support found, bailing out..\n";
        exit(1);
    }
    
    if(!is_extension_supported("GL_ARB_vertex_shader")) {
        std::cerr << "Extension GL_ARB_vertex_shader is required!"
                  << std::endl;
        return;
    }
    if(!is_extension_supported("GL_ARB_fragment_shader")) {
        std::cout << "Extension GL_ARB_fragment_shader is required!"
                  << std::endl;
        return;
    }

    if(!is_extension_supported("GL_EXT_multi_draw_arrays")) {
        std::cerr << "No glMultiDrawElements support, bailing out..\n";   
    }
}

XSurface_view::~XSurface_view() {

    if(glIsList(list_base_idx))
        glDeleteLists(list_base_idx, 4);
    glDeleteBuffers(n_VBOs, VBOs);
    
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
    glDeleteFramebuffersEXT(1, framebuf);

    if(depth_tex != 0) {
        glDeleteTextures(n_depth_layers, depth_tex);
        glDeleteTextures(n_depth_layers, image_tex);
        delete []depth_tex;
    }
}

//#endif // CGAL_USE_QT
