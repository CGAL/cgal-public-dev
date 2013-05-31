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
// File          : demos/xtriangulate/ogl_view.C
// QdX_release   : $Name:  $
// Revision      : $Revision: 1.8 $
// Revision_date : $Date: 2009-07-26 17:14:50 $
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

#include "include/ogl_view.h"
#include "include/skeletonizer_interface.h"
#include <fstream>

#include <GL/glut.h>

#include <qapplication.h>
#include <qcursor.h>
#include <qimage.h>

#include <boost/format.hpp>

#define MESH_VERTS   0
#define MESH_NORMALS 1
#define MESH_COLORS  2
#define MESH_INDICES 3
#define MESH_SIL_INDICES 4

#define NORM_VERTS 5

// #define CGAL_DEBUG_GL_LISTS 1

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
extern XSkeletonizer *skeletonizer; 

OGL_view::OGL_view(QWidget* parent, const char* name) :
            QGLWidget(parent, name) {

    surf_mesh_computed = false;

    lighting = true;
    transparency = false;
    bumpy_surf = false;
    solid_fill = true;
    black_bg_color = true;
    inverse_normals = false;
    show_triangles = true;
    show_grid = false;
    show_normals = false;

    surf_origin.assign(0);
    surf_scale = 1.0;    

    n_VBOs = 12;
    VBOs = new unsigned[n_VBOs]; 

    n_depth_layers = 5;
    depth_tex = new unsigned[n_depth_layers * 2];
    image_tex = depth_tex + n_depth_layers;

    z_translate = 0;
    width = parent->width();
    height = parent->height();
    reset_view();
}

unsigned texID[3]; QImage tex;

void OGL_view::initializeGL() {

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

    //glEnable(GL_CULL_FACE);
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

    lighting_shdr = load_shader("sample.vert", "sample.frag");
    antialias_shdr = load_shader("antialias.vert", "antialias.frag");
    //depth_peel_shdr = load_shader(0, "depth_peel.frag");
}

void OGL_view::reset_view() {
   
    transform = Vector_3d(0.0, 0.0, 1.0);
    quaternion_from_axis(q_rotate, Vector_3d(1, 0, 0), 45.0);
}

void OGL_view::setup_view() {

    glViewport(0, 0, width, height);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    
    double ratio = (double)width / height;
    //glOrtho(-ratio, ratio, -1, 1, 0.01, 10);
    gluPerspective(60, ratio, 0.01, 40.0);
    z_translate = -4.0;    // translate according to near/far clip planes

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
   
    glGenFramebuffersEXT(1, framebuf);
    glGenTextures(n_depth_layers, depth_tex);
    glGenTextures(n_depth_layers, image_tex);

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
    //glGenRenderbuffersEXT(1, &renderbuf);
    /*glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, renderbuf);
    glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_DEPTH_COMPONENT32,
        width, height);
    glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT,
        GL_RENDERBUFFER_EXT, renderbuf);*/
}

void OGL_view::paintGL() {

    if(transparency) {
        depth_peeled_render();
    } else {
        glEnable(GL_DEPTH_TEST);
        render_scene();
        
    }
    draw_axis();
}

void OGL_view::depth_peeled_render() {

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
    bool safe = show_triangles; // keep the old option
    show_triangles = false; // no triangles shown

    render_scene();
    glFlush();

    transparency = true, show_triangles = safe;
    
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

    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
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

#define GOLD 255, 215, 0
#define EMERALD 0x50, 0xc8, 0x78
#define NICE_BLUE 40, 80, 200
#define VIOLET  200, 10, 80
#define ORANGE 253, 192, 24
#define GREEN 50, 200, 100
#define BLUE 40, 100, 230

void OGL_view::render_scene() {

    if(black_bg_color)
        glClearColor(0.0, 0.0, 0.0, 1.0);
    else
        glClearColor(0.5, 0.5, 0.5, 1.0);
    
    glDepthMask(GL_TRUE);
    //glClearDepth(1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    glTranslated(transform.x()*fabs(z_translate - surf_origin[0]), 
        transform.y()*fabs(z_translate) - surf_origin[1], 
        z_translate - surf_origin[2]);

    glScaled(transform.z(), transform.z(), transform.z());
               
    double mat[16];
    quaternion_2_matrix(q_rotate, mat);
    glMultMatrixd(mat);
    draw_grid();
    
    glPolygonMode(GL_FRONT_AND_BACK, (solid_fill ? GL_FILL : GL_LINE));
    
    if(lighting && !transparency)
        glUseProgramObjectARB(lighting_shdr);
    
    GLint rescale_loc = glGetUniformLocation(lighting_shdr, "rescale_normals"),
        bumpy_loc = glGetUniformLocation(lighting_shdr, "bumpy_surf"),
        depth_peel_loc = glGetUniformLocation(lighting_shdr,
             "do_depth_peeling");
        
    glUniform1i(bumpy_loc, bumpy_surf);
    glUniform1i(depth_peel_loc, transparency);  

    //glDisable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    glUniform1f(rescale_loc, 1);

    glColor3ub(NICE_BLUE);
    glSecondaryColor3ub(ORANGE);

//     glColor3ub(VIOLET);
//     glSecondaryColor3ub(ORANGE);      // front-face color

//     glColor3ub(250, 120, 40);      // front-face color
//     glSecondaryColor3ub(200, 43, 100);  //
    
//     glSecondaryColor3ub(ORANGE);
    //glSecondaryColor3ub(40, 180, 220);  // ??
    render_surf();
     
//     if(transparency) {
//         glDisable(GL_BLEND);
//         glDisable(GL_CULL_FACE);
//     }
}

void OGL_view::render_surf() {

    if(skeletonizer->triangulation_valid() && 
            !surf_mesh_computed) {

        Point_vec_3d nd;
        Point_3d centroid;
        const Point_vec_3d& vd = skeletonizer->get_vertices();
        skeletonizer->compute_normals(vd, nd, centroid, inverse_normals);

        CGAL_assertion(vd.size() == nd.size());
        Point_vec_3f verts(vd.size()), norm_verts(vd.size()*2);
        Point_vec_3f normals(vd.size());

        float sc = 0.3f; // scale factor for normals
        // clamp data to single precision before loading to VBO
        for(unsigned i = 0; i < vd.size(); i++) {
            verts[i][0] = static_cast<float>(vd[i][0]);
            verts[i][1] = static_cast<float>(vd[i][1]);
            verts[i][2] = static_cast<float>(vd[i][2]);            

            normals[i][0] = static_cast<float>(nd[i][0]);
            normals[i][1] = static_cast<float>(nd[i][1]);
            normals[i][2] = static_cast<float>(nd[i][2]);

            norm_verts[2*i] = verts[i];
            norm_verts[2*i+1][0] = verts[i][0] + normals[i][0]*sc;
            norm_verts[2*i+1][1] = verts[i][1] + normals[i][1]*sc;
            norm_verts[2*i+1][2] = verts[i][2] + normals[i][2]*sc;
        }
         
        const Triangles_vector& indices = skeletonizer->get_triangles();
        const Triangles_vector& sil_indices =
                skeletonizer->get_sil_triangles();

        // _vidx,  _nidx,  _cidx,  _tidx, _iidx,  _idx_sz
        mesh_ctx = VBO_context(MESH_VERTS,  // vertex VBO index`
                               MESH_NORMALS,  // normals VBO index
                              -1,  // colors VBO index
                               MESH_INDICES,  // indices VBO index
                               indices.size() * 3,  // indices size
                               0, verts.size()-1);// lower/upper vertex indices

        sil_triangles_ctx = VBO_context(-1, // use MESH_VERTS VBO
                              -1,
                              -1,
                              MESH_SIL_INDICES,
                              sil_indices.size() * 3,
                              0, verts.size()-1);

        norms_ctx = VBO_context(NORM_VERTS,  // vertex VBO index`
                               -1,  // normals VBO index
                               -1,  // colors VBO index
                               -1,  // indices VBO index
                               0,  // indices size
                               0, norm_verts.size()-1);// lower/upper vertex indices


        load_VBO(mesh_ctx, &verts, &normals, 0, (unsigned *)indices.data());
        // indices already loaded with previous call
        load_VBO(sil_triangles_ctx, 0, 0, 0, (unsigned *)sil_indices.data());
        // drawing normals
        load_VBO(norms_ctx, &norm_verts, 0, 0, 0);

#if 0
        {
            typedef boost::array< unsigned short, 3 > Tri_short_3;
            std::vector< Tri_short_3 > idx_short(indices.size());
            
            for(unsigned i = 0; i < indices.size(); i++) {
                idx_short[i][0] = indices[i][0];
                idx_short[i][1] = indices[i][1];
                idx_short[i][2] = indices[i][2];
            }
         
            FILE *fp = fopen("output_model", "wb");

            printf("\n# verts: %d; # normals: %d; # tris: %d\n",
                verts.size(), normals.size(), indices.size());
            int nverts = verts.size(), ntris = indices.size();
//         fwrite(&verts.size(), sizeof(unsigned), 1, fp);
//         fwrite(&indices.size(), sizeof(unsigned), 1, fp);
//         fprintf(fp, "%d;%d", verts.size()*3, indices.size()*3);
            fwrite(&nverts, sizeof(int), 1, fp);
            fwrite(&ntris, sizeof(int), 1, fp);

            fwrite(verts.data(), verts.size()*3*sizeof(float), 1, fp);
            fwrite(tri_verts.data(), tri_verts.size()*3*sizeof(float), 1, fp);
            fwrite(normals.data(), normals.size()*3*sizeof(float), 1, fp);
            fwrite(idx_short.data(), idx_short.size()*3*sizeof(unsigned short), 1, fp);

            fclose(fp);
        }
#endif

        std::cout << "triangles loaded: " << indices.size() << "\n";
        surf_mesh_computed = true;
    }
    
    if(surf_mesh_computed) {

        bind_VBO(mesh_ctx); // attention: bind_VBO also pushes client attrib
                            // onto stack
        if(transparency || !show_triangles) {
            glDisable(GL_POLYGON_OFFSET_LINE);
            glDisable(GL_POLYGON_OFFSET_FILL);
        }

#if 0
    static unsigned char gl_colors[][3] =
        {{255, 215, 0},
         {0x50, 0xc8, 0x78},
        {40, 80, 200}, {200, 10, 80},
        {253, 192, 24}, {50, 200, 100},
        {40, 100, 230}};
    static unsigned n_gl_colors = sizeof(gl_colors) / 3;

    const std::vector< unsigned >& irange = skeletonizer->get_zlift_ranges();

    unsigned prev = 0;
    for(unsigned i = 0; i < irange.size(); i++) {

        const unsigned char *cl = gl_colors[i % n_gl_colors];
        glColor3ub(cl[0], cl[1], cl[2]);

        glDrawRangeElements(GL_TRIANGLES, mesh_ctx.min_vidx, mesh_ctx.max_vidx,
               (irange[i] - prev)*3, GL_UNSIGNED_INT,
                        (void *)(prev*3*4));
        prev = irange[i];
    }
#else        
        //glPolygonMode(GL_FRONT_AND_BACK, (solid_fill ? GL_FILL : GL_LINE));
        glDrawRangeElements(GL_TRIANGLES, mesh_ctx.min_vidx, mesh_ctx.max_vidx,
            mesh_ctx.index_size, GL_UNSIGNED_INT, 0);
#endif
                        
        if(transparency || !show_triangles)
            return;
            
        glUseProgramObjectARB(0);
        glPolygonMode(GL_FRONT_AND_BACK, (GL_LINE));
         
//         glEnable(GL_POLYGON_OFFSET_LINE);
//         glEnable(GL_POLYGON_OFFSET_FILL);
//         glPolygonOffset(10.0f, 1.0f);

        //glPopClientAttrib();
        glUseProgramObjectARB(0);
        
        glColor3ub(100, 100, 80);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glEnable(GL_LINE_SMOOTH);
        
        glDrawRangeElements(GL_TRIANGLES, mesh_ctx.min_vidx, mesh_ctx.max_vidx,
            mesh_ctx.index_size, GL_UNSIGNED_INT, 0);
        glDisable(GL_BLEND);
//         glDisable(GL_POLYGON_OFFSET_LINE);
//         glDisable(GL_POLYGON_OFFSET_FILL);

        glColor3f(0.6f, 0.8f, 0.6f);
        glLineWidth(0.7);

        bind_VBO(sil_triangles_ctx);
        glDrawRangeElements(GL_TRIANGLES, sil_triangles_ctx.min_vidx,
            sil_triangles_ctx.max_vidx,
            sil_triangles_ctx.index_size, GL_UNSIGNED_INT, 0);

        if(show_normals) {
            glColor3f(0.4f, 0.5f, 0.9f);
            bind_VBO(norms_ctx);
            glDrawArrays(GL_LINES, norms_ctx.min_vidx, norms_ctx.max_vidx);
        }
    }
}

//! saves solid model in STL format
void OGL_view::write_STL(std::ostream& os) {

    const Point_vec_3d& verts = skeletonizer->get_vertices();
    const Triangles_vector& indices = skeletonizer->get_triangles();

    os.precision(14);
    os.setf(std::ios::scientific, std::ios::floatfield);    
    os << "solid fluffy_bunny\n";
    for(unsigned i = 0; i < indices.size(); i++) {
        unsigned i1 = indices[i][0], i2 = indices[i][1], i3 = indices[i][2];

        Vector_3d v1(verts[i1]), v2(verts[i2]), v3(verts[i3]);
        Vector_3d v12 = v2 - v1, v13 = v3 - v1;       
        Vector_3d n = v12.cross_product(v13);
        n.normalize();

        os << "facet normal " << n.x() << " " << n.y() << " " << n.z() << 
            "\nouter loop\n";
        os << "vertex " << v1.x() << " " << v1.y() << " " << v1.z() << "\n";
        os << "vertex " << v2.x() << " " << v2.y() << " " << v2.z() << "\n";
        os << "vertex " << v3.x() << " " << v3.y() << " " << v3.z() << "\n";
        os << "endloop\nendfacet\n";
    }
    os << "endsolid fluffy_bunny";
}
//! renders \c count primitives from \c ctx 
void OGL_view::bind_VBO(const VBO_context& ctx) {
    
    if(ctx.verts_idx != -1u) {
        glPushClientAttrib(GL_CLIENT_VERTEX_ARRAY_BIT);
        glBindBuffer(GL_ARRAY_BUFFER, VBOs[ctx.verts_idx]);
        glVertexPointer(3, GL_FLOAT, 0, 0);
    }

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
      
    if(ctx.index_idx != -1u)   
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, VBOs[ctx.index_idx]);
}

void OGL_view::load_VBO(const VBO_context& ctx,
        const Point_vec_3f* verts, const Point_vec_3f *normals, 
        const std::vector< Color_ub_4 > *colors, const unsigned *indices) {

    if(ctx.verts_idx != -1u && verts != 0) {
        glBindBuffer(GL_ARRAY_BUFFER, VBOs[ctx.verts_idx]);
        glBufferData(GL_ARRAY_BUFFER, verts->size()*3*sizeof(float),
            verts->data(), GL_STATIC_DRAW);
    }

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

    if(ctx.index_idx != -1u && indices != 0) {
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, VBOs[ctx.index_idx]);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER,
            ctx.index_size*sizeof(unsigned), indices, GL_STATIC_DRAW);
    }
}

void OGL_view::world_coord(int x, int y, Point_2d& res) {

    res[0] = (2.0*x - width) / width, 
    res[1] = (height - 2.0*y) / height;
}

void OGL_view::mousePressEvent(QMouseEvent* me) {
    start_pt = last_pt = me->pos();
    which_button = me->button();
    q_start = q_rotate;
}

void OGL_view::mouseMoveEvent(QMouseEvent* me) {
    
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

void OGL_view::mouseReleaseEvent(QMouseEvent* me) {
    setCursor(QCursor(Qt::ArrowCursor));
}

void OGL_view::wheelEvent(QWheelEvent *we) {
    
    double dir = (we->delta() > 0 ? 0.95 : 1.05);
    transform[2] *= dir;
    //std::cerr << transform[2] << "\n";
    updateGL();
}

void OGL_view::draw_grid() {

    if(!show_grid)
        return;    

    int n_cells = 10, i;
    float ext = 0/*triangulator_engine->get_base_surf_extent()*/,
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

void OGL_view::draw_axis() {

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

void OGL_view::resizeGL(int w, int h) {
    width = w;
    height = h;
    setup_view();
    updateGL();
}

void OGL_view::check_extensions() {
    
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

OGL_view::~OGL_view() {

    if(glIsList(list_base_idx))
        glDeleteLists(list_base_idx, 4);
    glDeleteBuffers(n_VBOs, VBOs);

    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
    glDeleteFramebuffersEXT(1, framebuf);
        
    if(VBOs != 0)
        delete []VBOs;
    if(depth_tex != 0) {
        glDeleteTextures(n_depth_layers, depth_tex);
        glDeleteTextures(n_depth_layers, image_tex);
        delete []depth_tex;
    }
}

#include "ogl_view.moc"
