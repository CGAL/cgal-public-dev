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
// File          : demos/xsurface/ogl_exts.C
// QdX_release   : $Name:  $
// Revision      : $Revision$
// Revision_date : $Date$
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

/*!@file ogl_exts.C
 * extensions support for OpenGL
 */

#include "include/ogl_view.h"

// #ifdef CGAL_USE_QT

char *load_text(const char* filename) {

    std::ifstream is;
    is.open(filename, std::ios::binary);
    if(!is) 
        return 0;

    is.seekg(0, std::ios::end);
    unsigned len = is.tellg();
    is.seekg(0, std::ios::beg);
    if(len == 0) 
        return 0; 

    char *buf = new char[len+1];
    if(buf == 0) 
        return 0;
     
    is.read(buf, len);
    buf[len] = 0;
    is.close();
    return buf;
}

void shader_compiler_log(GLhandleARB prog_obj) {
    if(prog_obj == 0)
        return;    

    int len, slen;
    glGetObjectParameterivARB(prog_obj, GL_OBJECT_INFO_LOG_LENGTH_ARB, &len);
    if(len > 1) {
        char* compiler_log = new char[len];
        if(compiler_log == NULL) {
            std::cout << "Shader compiler log: out of memory!" << std::endl;
            return;
        }

        glGetInfoLogARB(prog_obj, len, &slen, compiler_log);
        if(compiler_log != 0) {
            std::cout << "Shader compiler errlog: " << std::endl;
            std::cout << compiler_log << std::endl;
        }
        delete []compiler_log;
    }
}

GLhandleARB XSurface_view::load_shader(const char* vs_filename,
        const char* fs_filename) {

    char* vs_src = 0, *fs_src = 0;
    if(vs_filename != 0)
        if((vs_src = load_text(vs_filename)) == 0)
            std::cerr << "Unable to load vertex shader: " << vs_filename <<
                    std::endl;

    if(fs_filename != 0)
        if((fs_src = load_text(fs_filename)) == 0)
            std::cerr << "Unable to load fragment shader: " << fs_filename <<
                    std::endl;

    GLhandleARB vs_obj = 0, fs_obj = 0, prog_obj = 0;
    GLint vs_compiled = 0, fs_compiled = 0, linked = 0;

    if(vs_src != 0) {
        vs_obj = glCreateShaderObjectARB(GL_VERTEX_SHADER_ARB);
        int len = strlen(vs_src);
        glShaderSourceARB(vs_obj, 1, (const GLcharARB**)&vs_src, &len);
        glCompileShaderARB(vs_obj);
        glGetObjectParameterivARB(vs_obj, GL_OBJECT_COMPILE_STATUS_ARB,
                                  &vs_compiled);
        if(!vs_compiled) 
            shader_compiler_log(vs_obj);
        delete []vs_src;
    } else {
        vs_obj = 0;
        vs_compiled = 1;
    }

    if(fs_src != 0) {
        fs_obj = glCreateShaderObjectARB(GL_FRAGMENT_SHADER_ARB);
        int len = strlen(fs_src);
        glShaderSourceARB(fs_obj, 1, (const GLcharARB**)&fs_src, &len);
        glCompileShaderARB(fs_obj);
        glGetObjectParameterivARB(fs_obj, GL_OBJECT_COMPILE_STATUS_ARB,
                                  &fs_compiled);
        if(!fs_compiled) 
            shader_compiler_log(fs_obj);
        delete []fs_src;
    } else {
        fs_obj = 0;
        fs_compiled = 1;
    }

    if(!vs_compiled || !fs_compiled) {
        if(vs_obj) 
            glDeleteObjectARB(vs_obj);
        if(fs_obj) 
            glDeleteObjectARB(fs_obj);
        return 0;
    }

    prog_obj = glCreateProgramObjectARB();
    if(vs_obj) {
        glAttachObjectARB(prog_obj, vs_obj);
        glDeleteObjectARB(vs_obj);
    }

    if(fs_obj) {
        glAttachObjectARB(prog_obj, fs_obj);
        glDeleteObjectARB(fs_obj);
    }

    glLinkProgramARB(prog_obj);
    glGetObjectParameterivARB(prog_obj, GL_OBJECT_LINK_STATUS_ARB, &linked);

    if(!linked) {
        shader_compiler_log(prog_obj);
        glDeleteObjectARB(prog_obj);
        prog_obj = 0;
    }

    return prog_obj;
}

bool XSurface_view::is_extension_supported(const char *name) {
    
    if(!name || *name == '\0') 
        return false;
 
    char *where;   
    const char *start = (const char *)glGetString(GL_EXTENSIONS);
    
    for (;;) {
        where = strstr(start, name);
        if(!where)
            break;
        char *next = where + strlen(name);
        if(where == start || *(where - 1) == ' ')
            if(*next == ' ' || *next == '\0')
                return true;
        start = next;
    }
    return false;
}

// #endif // CGAL_USE_QT
