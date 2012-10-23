// ============================================================================
//
// Copyright (c) 2001-2012 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of EXACUS (http://www.mpi-inf.mpg.de/projects/EXACUS/);
// you may redistribute it under the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with EXACUS.f
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// ----------------------------------------------------------------------------
//
// Library       : 
// File          : rendering_gl.js
//
// Author(s)     : Pavel Emeliyanenko: http://www.mpi-inf.mpg.de/~emeliyan
//
// ============================================================================

var webXtri = { };
var gl = null;

var framebuf0, framebuf1, mesh0, mesh1;
var verts_buf, normals_buf, tri_indices_buf, indices_buf, sil_indices_buf;
// var framebuf;
// var n_depth_layers = 5;
// var depth_tex, image_tex; // arrays for depth peeling

// projection, model-view and normal matrices
var proj_matrix, mview_matrix, normal_matrix;
// reference to shader uniform variable
var mview_matrix_loc, normal_matrix_loc;  
var verts_pos_loc, verts_color_loc,
  verts_normal_loc; // references to shader attribute variables

var color_loc, back_color_loc, lighting_loc, normals_loc;
var front_color, back_color, clear_color;

var mesh_loaded = 0, n_indices = 0, n_sil_indices = 0;
var show_triangles = false, inverse_normals = false,
    show_silhouette = false;

const MATH_ROUND = 10000.0;
webXtri.left = -2.0, webXtri.right = 2.0;
webXtri.btm = -2.0, webXtri.top = 2.0;
webXtri.centroid = [0, 0, 0];

webXtri.def_shader = null;
webXtri.box_shader = null;

var surfaceID = null;
var btri_disabled = true;
var n_z_patches = 0;
var z_patches = new Array;

function start() {

    webXtri.canvas = document.getElementById("glcanvas");
    gl = WebGLUtils.setupWebGL(webXtri.canvas);
    if(!gl) {
        return;
    }

    surfaceID = new Array(-1, -1, -1, -1);
    front_color = new Array(40.0/255.0, 80.0/255.0, 200.0/255.0);
    back_color = new Array(1.0, 0.3, 0.1);
    clear_color = new Array(0.0, 0.0, 0.0);

    webXtri.setup_user_interface();

    gl.viewport(0, 0, webXtri.canvas.width, webXtri.canvas.height);
    gl.clearColor(clear_color[0], clear_color[1], clear_color[2], 1.0);
    gl.clearDepth(1.0);
    gl.enable(gl.DEPTH_TEST);
    gl.depthFunc(gl.LEQUAL);

//     gl.enable(gl.CULL_FACE);
    gl.cullFace(gl.BACK);
    gl.frontFace(gl.CCW);
    
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

/*
RKgl.beginTFBO = function(tfbo) {
    gl.bindFramebuffer(gl.FRAMEBUFFER, tfbo.fbo);
}

RKgl.endTFBO = function(tfbo) {
   gl.bindTexture(gl.TEXTURE_2D, tfbo.texture);
   gl.generateMipmap(gl.TEXTURE_2D);
   gl.bindTexture(gl.TEXTURE_2D, null);
   gl.bindFramebuffer(gl.FRAMEBUFFER, null);
}*/
    var w = webXtri.canvas.width, h = webXtri.canvas.height;
    framebuf0 = new webXtri.TextureFrameBuffer(w, h);
    framebuf1 = new webXtri.TextureFrameBuffer(w, h);

    mesh0 = new webXtri.RectangleMesh
//             (0,0,799,599);
            (-1, -1, 2, 2);

    verts_buf = gl.createBuffer();
    normals_buf = gl.createBuffer();
    indices_buf = gl.createBuffer();
    sil_indices_buf = gl.createBuffer();
    tri_indices_buf = gl.createBuffer();

//     var canvas = document.getElementById("glcanvas");
//     var w = canvas.width, h = canvas.height;
    
    webXtri.add_external_shader("def_vert", "shader_def_vert.c?xx1");
    webXtri.add_external_shader("def_frag", "shader_def_frag.c?xx1");

    webXtri.add_external_shader("box_vert", "shader_box_vert.c");
    webXtri.add_external_shader("box_frag", "shader_box_frag.c");

    webXtri.add_external_shader("blur_vert", "shader_blur_vert.c");
    webXtri.add_external_shader("blur_frag", "shader_blur_frag.c");

    webXtri.load_external_shaders(function(){ init_shaders() });
}

function init_shaders() {

//    webXtri.def_shader = webXtri.embedded_shader_prog(gl,
//                               "shader-vs", "shader-fs");

    webXtri.def_shader = webXtri.external_shader_prog(gl,
                            "def_vert", "def_frag");

    webXtri.box_shader = webXtri.external_shader_prog(gl,
                            "box_vert", "box_frag");

    webXtri.blur_shader = webXtri.external_shader_prog(gl,
                            "blur_vert", "blur_frag");

    var prog = webXtri.def_shader;
    gl.useProgram(prog);

    proj_matrix = mat4.create();
    mat4.perspective(60, 1.333333, 0.1, 100, proj_matrix);
    
    gl.uniformMatrix4fv(gl.getUniformLocation(prog, "proj_matrix"),
       false, proj_matrix);

    mview_matrix = mat4.create(); 
    mview_matrix_loc = gl.getUniformLocation(prog, "mview_matrix");

    normal_matrix = mat4.create(); 
    normal_matrix_loc = gl.getUniformLocation(prog, "normal_matrix");

    q_rotate = quaternion_from_axis([1.0, 1.0, 1.0], 45);
    q_start = quat4.create(q_rotate);

    verts_pos_loc = gl.getAttribLocation(prog, "aPos");
    gl.enableVertexAttribArray(verts_pos_loc);
//     verts_color_loc = gl.getAttribLocation(prog, "aVertexColor");
//     gl.enableVertexAttribArray(verts_color_loc);
    verts_normal_loc = gl.getAttribLocation(prog, "aNormal");
    gl.enableVertexAttribArray(verts_normal_loc);

    color_loc = gl.getUniformLocation(prog, "uColor");
    back_color_loc = gl.getUniformLocation(prog, "uBackColor");
    lighting_loc = gl.getUniformLocation(prog, "lighting_mode");
    normals_loc = gl.getUniformLocation(prog, "inverse_normals");
}

webXtri.TextureFrameBuffer = function(w, h) {

    this.width = w;
    this.height = h;

    this.fbo = gl.createFramebuffer();
    gl.bindFramebuffer(gl.FRAMEBUFFER, this.fbo);

    this.fbo.width = this.width;
    this.fbo.height = this.height;
    this.texture = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, this.texture);

    //gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
    //gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER,
  //          gl.LINEAR_MIPMAP_NEAREST);
    //gl.generateMipmap(gl.TEXTURE_2D);

    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);

    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);

    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, this.fbo.width, this.fbo.height,
           0, gl.RGBA, gl.UNSIGNED_BYTE, null);

    var renderbuffer = gl.createRenderbuffer();
    gl.bindRenderbuffer(gl.RENDERBUFFER, renderbuffer);
    gl.renderbufferStorage(gl.RENDERBUFFER, gl.DEPTH_COMPONENT16,
           this.fbo.width, this.fbo.height);

    gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0,
                gl.TEXTURE_2D, this.texture, 0);
    gl.framebufferRenderbuffer(gl.FRAMEBUFFER, gl.DEPTH_ATTACHMENT,
                gl.RENDERBUFFER, renderbuffer);

    gl.bindTexture(gl.TEXTURE_2D, null);
    gl.bindRenderbuffer(gl.RENDERBUFFER, null);
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    return this;
}

webXtri.RectangleMesh = function(x0, y0, w, h, color, texture) {
    
//     var rectPrimitive = new RKPrimitives3.Rectangle(cx,cy,w,h,color);
    
    var verts = [x0, y0, 0,  x0+w, y0, 0,
                 x0+w, y0+h, 0,  x0, y0+h, 0];
    var idxs = [0, 1, 2,  0, 2, 3];
    var texcoords = [0.0, 0.0,  1.0, 0.0,  1.0, 1.0,  0.0, 1.0];

//     var colors = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
//                   1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
//     if (!color) {
//         color = {r:1.0,g:1.0,b:1.0,a:1.0}
//     }
// 
//     for(var i =0;i<4;i++) {
//        this.colors = this.colors.concat([color.r,color.g,color.b,color.a]);
// 
//     }

    this.vertex_buf = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, this.vertex_buf);
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(verts), gl.STATIC_DRAW);

    this.tex_buf = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, this.tex_buf);
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(texcoords),
                  gl.STATIC_DRAW);
    
    this.index_buf = gl.createBuffer();
    gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, this.index_buf);
    gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(idxs),
                  gl.STATIC_DRAW);
    this.n_indices = idxs.length;
}

/**


var draw = function() {

    gl.disable(gl.DEPTH_TEST);
   RKgl.beginTFBO(fb0);

   gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
   RKgl.clearColor(1.0, 1.0, 1.0, 1.0);
   RKgl.clear();
   RKgl.currentShaderProgram = boxProgSil.program;
   RKgl.drawMesh(modelMesh);
   RKgl.endTFBO(fb0);

   step=.005;
   RKgl.beginTFBO(fb1);
   RKgl.clearColor(0.0, 0.0, 0.0, 1.0);
   RKgl.clear();
   RKgl.currentShaderProgram = surfaceProg.program;
   RKgl.drawMesh(surface,function() {
        gl.uniform1f(RKgl.currentShaderProgram.uTime,time)
        gl.uniform1f(RKgl.currentShaderProgram.uStep,step)
   });
   RKgl.endTFBO(fb1);

   step=.01;
   RKgl.drawMesh(surface2,function() {
        gl.uniform1f(RKgl.currentShaderProgram.uTime,time)
        gl.uniform1f(RKgl.currentShaderProgram.uStep,step)
   });
   RKgl.currentShaderProgram = boxProg.program;
   //gl.enable(gl.DEPTH_TEST);
   gl.blendFunc(gl.ONE, gl.ONE);
   //gl.blendFunc(gl.CONSTANT_ALPHA, gl.ONE_MINUS_CONSTANT_COLOR);
   RKgl.drawMesh(modelMesh);

}

*/

webXtri.render_object = function(mesh, prog, texture) {

    gl.bindBuffer(gl.ARRAY_BUFFER, mesh.vertex_buf);

    var verts_pos = gl.getAttribLocation(prog, "aPos");
    gl.enableVertexAttribArray(verts_pos);
    gl.vertexAttribPointer(verts_pos, 3, gl.FLOAT, false, 0, 0);

    gl.bindBuffer(gl.ARRAY_BUFFER, mesh.tex_buf);
    var texs_pos = gl.getAttribLocation(prog, "aTexCoord");
    gl.enableVertexAttribArray(texs_pos);
    gl.vertexAttribPointer(texs_pos, 2, gl.FLOAT, false, 0, 0);

    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, texture);

    var sampler_pos = gl.getUniformLocation(prog, "uSampler");
    gl.uniform1i(sampler_pos, 0);

    gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, mesh.index_buf);

    var proj = mat4.create(); 
    mat4.identity(proj);
//     proj.ortho(0, 800, 0, 600, -1, 1);

    gl.uniformMatrix4fv(gl.getUniformLocation(prog, "proj_matrix"),
       false, proj);

    mat4.identity(mview_matrix);
    mat4.translate(mview_matrix, [0.0, 0.0, -0.5]);
    gl.uniformMatrix4fv(gl.getUniformLocation(prog, "mview_matrix"),
        false, mview_matrix);

    gl.drawElements(gl.TRIANGLES, mesh.n_indices, gl.UNSIGNED_SHORT, 0);

/*    document.getElementById("progress_msg").innerHTML = texs_pos + ' ' +
        verts_pos + ' ' + sampler_pos + ' ' +  texture + ' ' + prog;*/
//     RKgl.mvPopMatrix();
}
       
function draw_scene() {

    if(mesh_loaded == 0)
        return;

    gl.disable(gl.BLEND);
    if(is_any_z_patch_selected()) {

        gl.enable(gl.DEPTH_TEST);
//         gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
//     gl.enable(gl.CULL_FACE);

        gl.bindFramebuffer(gl.FRAMEBUFFER, framebuf0.fbo);
    
//     gl.clearColor(1.0, 1.0, 1.0, 1.0);
        gl.clearColor(0.0, 0.0, 0.0, 0.7);
    
        gl.useProgram(webXtri.def_shader);
        draw_surface_mesh(webXtri.def_shader, true);
    }
//     gl.bindTexture(gl.TEXTURE_2D, framebuf0.texture);
//     gl.generateMipmap(gl.TEXTURE_2D);

/*
    gl.disable(gl.DEPTH_TEST);

    // first blur pass
    gl.bindFramebuffer(gl.FRAMEBUFFER, framebuf1.fbo);
    gl.clearColor(1.0, 1.0, 1.0, 1.0);

    gl.useProgram(webXtri.blur_shader);
    gl.uniform1f(gl.getUniformLocation(webXtri.blur_shader, "uStep"), 0.005);
    webXtri.render_object(mesh0, webXtri.blur_shader, framebuf0.texture);*/
    
//     gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA)
    gl.bindTexture(gl.TEXTURE_2D, null);
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);

    gl.clearColor(clear_color[0], clear_color[1], clear_color[2], 1.0);
    gl.enable(gl.DEPTH_TEST);
    gl.useProgram(webXtri.def_shader);
    draw_surface_mesh(webXtri.def_shader, false);

    if(is_any_z_patch_selected()) {
        gl.disable(gl.DEPTH_TEST);
//     gl.clearColor(0.0, 0.0, 0.0, 1.0);
//     gl.clear(gl.COLOR_BUFFER_BIT);

        gl.enable(gl.BLEND);
        gl.blendFunc(gl.ONE, gl.ONE);
        gl.useProgram(webXtri.blur_shader);

        gl.uniform1f(gl.getUniformLocation(webXtri.blur_shader, "uStep"), 0.01);
        webXtri.render_object(mesh0, webXtri.blur_shader, framebuf0.texture);
    }
    
    gl.flush();
}

function draw_surface_mesh(prog, blur_image) {

    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

    gl.bindBuffer(gl.ARRAY_BUFFER, verts_buf);
    gl.vertexAttribPointer(verts_pos_loc, 3, gl.FLOAT, false, 0, 0);
  
//     gl.bindBuffer(gl.ARRAY_BUFFER, colors_buf);
//     gl.vertexAttribPointer(verts_color_loc, 3, gl.FLOAT, false, 0, 0);

    gl.bindBuffer(gl.ARRAY_BUFFER, normals_buf);
    gl.vertexAttribPointer(verts_normal_loc, 3, gl.FLOAT, false, 0, 0);

    gl.uniformMatrix4fv(gl.getUniformLocation(prog, "proj_matrix"),
       false, proj_matrix);

//     var mm = mat4.create();
//     mat4.identity(mm);
//     mat4.rotateX(mm, (xRot/5)/ 180 * Math.PI);
//     mat4.rotateY(mm, (yRot/5)/ 180 * Math.PI);
//     mat4.multiply(mm, rotate_mat, rotate_mat);

    var q = quaternion_trackball(vStart, vRot);
    quat4.multiply(q_start, q, q_rotate);
    quat4.normalize(q_rotate);

    var rotate_mat = mat4.create();
    quat4.toMat4(q_rotate, rotate_mat);
    
//     debug_msg("mv1: " + mat4.str(rotate_mat) + "<br />mv2: " +
//            (vStart[0] - vOffs[0]) + ' ' + (vStart[1] - vOffs[1]));

    mat4.identity(mview_matrix);
    mat4.translate(mview_matrix, translate);
    mat4.multiply(mview_matrix, rotate_mat);
    mat4.translate(mview_matrix, webXtri.centroid);

    gl.uniformMatrix4fv(mview_matrix_loc, false,
                       mview_matrix);

    mat4.set(rotate_mat, normal_matrix);
    mat4.inverse(normal_matrix);
    mat4.transpose(normal_matrix);
    gl.uniformMatrix4fv(normal_matrix_loc, false, //set true for transpose
                            normal_matrix);

    gl.uniform1i(normals_loc, inverse_normals);
    gl.uniform1i(lighting_loc, 0);
    gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, indices_buf);

    var obj = document.getElementById("z-patches");
    var sel = document.getElementById("show_selected").checked;
    
    if(!blur_image && !sel) {
        
        gl.uniform4f(color_loc, front_color[0], front_color[1],
                 front_color[2], 1.0);

        gl.uniform4f(back_color_loc, back_color[0], back_color[1],
                 back_color[2], 1.0);
        // here we give the actual array size (not the number of triangles!)
        gl.drawElements(gl.TRIANGLES, n_indices, gl.UNSIGNED_SHORT, 0);
        // offset: is a *byte* offset !!
    }

    if(blur_image || sel) {
    for(var i = 0; i <= n_z_patches; i++) {

        if(sel && !obj.options[i].selected) {
            continue;
        }
        
        var c_amp = 1.7; // color amplify factor
        if(!blur_image)
            c_amp = 1.0;
        
        if(!obj.options[i].selected) {
            continue;

        } else if(obj.options[i].selected) {
//             alert("selected!\n");
        }

        gl.uniform4f(color_loc, front_color[0]*c_amp, front_color[1]*c_amp,
                 front_color[2]*c_amp, 1.0);

        gl.uniform4f(back_color_loc, back_color[0]*c_amp, back_color[1]*c_amp,
                 back_color[2]*c_amp, 1.0);
        
        gl.drawElements(gl.TRIANGLES, (z_patches[i+1]-z_patches[i])*3,
                  gl.UNSIGNED_SHORT, z_patches[i]*3*2);
   }
       if(blur_image) 
            return;
   } // blur_image


    if(show_triangles) {
        gl.lineWidth(1.0);
        gl.uniform1i(lighting_loc, 1);
        gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, tri_indices_buf);
        gl.drawElements(gl.LINES, n_indices*2, gl.UNSIGNED_SHORT, 0);
    }

    if(show_silhouette) {
        gl.lineWidth(2.0);
        gl.uniform1i(lighting_loc, 2);
        gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, sil_indices_buf);
        gl.drawElements(gl.LINES, n_sil_indices*2, gl.UNSIGNED_SHORT, 0);
    }
}

// MD5 checksum: 4 words
// struct Triangulation_info {
//     uint n_verts, n_tris, n_sil_tris, n_z_patches;  // 4 words
//     double left, right, btm, top; // 8 words
//     double cx, cy, cz;            // 6 words
// };

function parse_model(stream) {
  
//     var tm_start = new Date().getTime();
    var len = stream.length;
    var header_size = 4 + 4 + 8 + 6; // header size in 32-bit words

    var buf_sz = len + 16;
    if(buf_sz < header_size * 4)
        buf_sz = header_size * 4;
        
    var buf = new ArrayBuffer(buf_sz);
    var bytes = new Uint8Array(buf, 0, buf_sz);
    var ints = new Int32Array(buf, 0, header_size);
    var dbls = new Float64Array(buf, 0, header_size/2);
    
    for(i = 0; i < len; i++) {
        bytes[i] = stream.charCodeAt(i) & 0xff;
    }
    surfaceID[0] = ints[0]; surfaceID[1] = ints[1];
    surfaceID[2] = ints[2]; surfaceID[3] = ints[3];

    if(surfaceID[0] == -1 && surfaceID[1] == -1 &&
        surfaceID[2] == -1 && surfaceID[3] == -1) {
        report_error(ints[4]);
        mesh_loaded = -1; // error model
        return;
    }
//     alert(dbls[4] + ' ' + dbls[5] + ' ' + dbls[6] + ' ' + dbls[7]);
    webXtri.left = dbls[4]; webXtri.right = dbls[5];
    webXtri.btm = dbls[6]; webXtri.top = dbls[7];
    webXtri.centroid[0] = -dbls[8];
    webXtri.centroid[1] = -dbls[9];
    webXtri.centroid[2] = -dbls[10];

/*    alert(webXtri.centroid[0] + ' ' + webXtri.centroid[1] + ' ' +
        webXtri.centroid[2]);*/
    
    var n_verts = ints[4], n_tris = ints[5], n_sil_tris = ints[6];
    n_z_patches = ints[7];

    if(isNaN(webXtri.left) || isNaN(webXtri.right) || isNaN(webXtri.btm) ||
        isNaN(webXtri.top) || webXtri.left >= webXtri.right ||
        webXtri.btm >= webXtri.top || isNaN(n_verts) ||
          isNaN(n_tris) || isNaN(n_sil_tris) || isNaN(n_z_patches)) {
        report_error(ERR_INVALID_MESH);
        mesh_loaded = -1; // error model
        return;
    }
        
    var bytes_skip = header_size * 4 + n_z_patches*4;
    var sz = (n_verts * 2)* 3 * 4 + (n_tris * 3 + n_sil_tris * 3) * 2;
    if(sz != len - bytes_skip) { // do not count header 
//         alert(sz + ' ' + len);
        report_error(ERR_INVALID_MESH);
        mesh_loaded = -1; // error model
        return;
    }
    report_error(ERR_OK);

    var ints2 = new Int32Array(buf, header_size * 4, n_z_patches);
    // pupolate z-patches select box
    set_z_patches(n_z_patches, ints2);
    
    // update absolute bounds respectively
    set_absolute_bounds();
    
    n_indices = n_tris*3;
    n_sil_indices = n_sil_tris*3;

    // ofs - byte offset in ArrayBuffer
    // length - # of elements (Floats!)
    var verts = new Float32Array(buf, bytes_skip, n_verts*3);
    var normals = new Float32Array(buf, bytes_skip + n_verts*12, n_verts*3);
    var indices = new Uint16Array(buf, bytes_skip + n_verts*24, n_indices);
    var sil_indices = new Uint16Array(buf, bytes_skip +
            n_verts*24 + n_indices*2, n_sil_indices);
    // n_indices * 2 * sizeof(short)
    var tri_buf = new ArrayBuffer(n_indices*2* 2);
    // use 3 lines for each triangle
    var tri_indices = new Uint16Array(tri_buf, 0, n_indices*2);

    for(i = 0; i < n_tris; i++) {
        var i1 = indices[3*i], i2 = indices[3*i+1], i3 = indices[3*i+2];
        var i6 = 6*i;
        tri_indices[i6] = i1; tri_indices[i6+1] = i2;
        tri_indices[i6+2] = i2; tri_indices[i6+3] = i3;
        tri_indices[i6+4] = i3; tri_indices[i6+5] = i1;
    }
    
    var sil_tri_buf = new ArrayBuffer(n_sil_indices*2* 2);
    // use 3 lines for each triangle
    var sil_tri_indices = new Uint16Array(sil_tri_buf, 0, n_sil_indices*2);
    for(i = 0; i < n_sil_tris; i++) {
        var i1 = sil_indices[3*i], i2 = sil_indices[3*i+1],
                i3 = sil_indices[3*i+2];
        var i6 = 6*i;
        sil_tri_indices[i6] = i1; sil_tri_indices[i6+1] = i2;
        sil_tri_indices[i6+2] = i2; sil_tri_indices[i6+3] = i3;
        sil_tri_indices[i6+4] = i3; sil_tri_indices[i6+5] = i1;
    }
    
    gl.bindBuffer(gl.ARRAY_BUFFER, verts_buf);
    gl.bufferData(gl.ARRAY_BUFFER, verts, gl.STATIC_DRAW);

    gl.bindBuffer(gl.ARRAY_BUFFER, normals_buf);
    gl.bufferData(gl.ARRAY_BUFFER, normals, gl.STATIC_DRAW);

    gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, indices_buf);
    gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, indices, gl.STATIC_DRAW);

    gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, sil_indices_buf);
    gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, sil_tri_indices, gl.STATIC_DRAW);

    gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, tri_indices_buf);
    gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, tri_indices, gl.STATIC_DRAW);

    mesh_loaded = 1;
    draw_scene();
}

