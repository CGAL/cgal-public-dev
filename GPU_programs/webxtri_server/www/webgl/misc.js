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
// File          : misc.js
//
// Author(s)     : Pavel Emeliyanenko: http://www.mpi-inf.mpg.de/~emeliyan
//
// ============================================================================

webXtri.currentShaderProgram = null;
 // collection of shaders to be loaded
webXtri.shader_sources = {};
 // callback for when all shaders are loaded
webXtri.onShadersLoaded = function() {};
webXtri.all_shaders_loaded = false;
webXtri.n_shaders_total = 0;
webXtri.n_shaders_loaded = 0;
webXtri.canvas = null;

var vOffs = [0,0], drag = 0;
var vStart = [0,0], vRot = [0,0];
var translate = [0,0,-8];
var server_id = -1;

var q_start, q_rotate;

const ERR_OK = 0,                   // no errors
      ERR_INVALID_POLYNOMIAL = -1,  // invalid polynomial format
      ERR_INVALID_DATA = -2,        // invalid client parameters
      ERR_TIMEOUT = -3,             // analyse/rasterize phase timeout
      ERR_INVALID_REFERRER = -6,    // an attempt to access the script from
                                   // outside the server
      ERR_SERVER_TIMEOUT = -7,      // no connection to the server
      ERR_REQUEST_PENDING = -8,     // a request from this IP is already being
                                   // processed by the server
      ERR_SERVER_OVERLOAD = -9,     // the server is overloaded (number of
                                   // requests processed >= MAX_CLIENTS
      ERR_TRIANGULATE_GENERIC = -10, // generic error during rasterize request
      ERR_INVALID_MESH = -11,        // unable to parse surface mesh
      ERR_GENERIC_ERROR = -12;    // unknown error


function set_surface_color(side, clr) {

    if(side == 0) {
        front_color[0] = clr.rgb[0];
        front_color[1] = clr.rgb[1];
        front_color[2] = clr.rgb[2];
    } else {
        back_color[0] = clr.rgb[0];
        back_color[1] = clr.rgb[1];
        back_color[2] = clr.rgb[2];
    }
    draw_scene();
}

function set_clear_color(clr) {
//     gl.clearColor(clr.rgb[0], clr.rgb[1], clr.rgb[2], 1.0);
//     gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
    clear_color = clr.rgb;
    draw_scene();
}

function toggle_show_triangles(checked) {
    show_triangles = checked;
    draw_scene();
}

function toggle_inverse_normals(checked) {
    inverse_normals = checked;
    draw_scene();
}

function toggle_show_silhouette(checked) {
    show_silhouette = checked;
    draw_scene();
}

function toggle_show_selected(checked) {
   draw_scene();
}

function toggle_use_abs_bounds(checked) {
    document.getElementById("left-id").disabled = checked;
    document.getElementById("right-id").disabled = checked;
    document.getElementById("btm-id").disabled = checked;
    document.getElementById("top-id").disabled = checked;

    document.getElementById("left-abs-id").disabled = !checked;
    document.getElementById("right-abs-id").disabled = !checked;
    document.getElementById("btm-abs-id").disabled = !checked;
    document.getElementById("top-abs-id").disabled = !checked;
}

function ajaxRequest() {
    // activeX versions to check for in IE
    var activexmodes=["Msxml2.XMLHTTP", "Microsoft.XMLHTTP"]
    // Test for support for ActiveXObject in IE first (as XMLHttpRequest in
    // IE7 is broken)
    if(window.ActiveXObject) {
        for(var i = 0; i < activexmodes.length; i++) {
        try {
            return new ActiveXObject(activexmodes[i]);
        }
        catch(e){
            //suppress error
        }
        }
    }
    else if (window.XMLHttpRequest) // if Mozilla, Safari etc
        return new XMLHttpRequest();
    else
        return false;
}

function check_prop(obj) {

    var s = "";
    for(var key in obj){
       s = s + key + ' = ' + obj[key] + '\n';
    }
    alert(s);
}

function get_window_size() {
    var winW = 630, winH = 460;
    if(document.body && document.body.offsetWidth) {
        winW = document.body.offsetWidth;
        winH = document.body.offsetHeight;
    }
    if(document.compatMode=='CSS1Compat' && document.documentElement &&
                document.documentElement.offsetWidth ) {
        winW = document.documentElement.offsetWidth;
        winH = document.documentElement.offsetHeight;
    }
    if(window.innerWidth && window.innerHeight) {
        winW = window.innerWidth;
        winH = window.innerHeight;
    }
    return winW;
}

webXtri.setup_user_interface = function() {

    document.getElementById('front_color').color.fromRGB(
        front_color[0], front_color[1], front_color[2]);

    document.getElementById('back_color').color.fromRGB(
        back_color[0], back_color[1], back_color[2]);

    document.getElementById('clear_color').color.fromRGB(
        clear_color[0], clear_color[1], clear_color[2]);

    document.getElementById("banalyse").disabled = false;
    document.getElementById("btriangulate").disabled = btri_disabled;

    show_triangles = document.getElementById("show_tris").checked;
    inverse_normals = document.getElementById("inv_normals").checked;
    show_silhouette = document.getElementById("show_sil").checked;
    document.getElementById("use_abs").checked = false;
    toggle_use_abs_bounds(false);

    document.getElementById("left-abs-id").value = webXtri.left;
    document.getElementById("right-abs-id").value = webXtri.right;
    document.getElementById("btm-abs-id").value = webXtri.btm;
    document.getElementById("top-abs-id").value = webXtri.top;

    webXtri.canvas.resize = function (){
        var size = Math.min(window.innerWidth, window.innerHeight) - 15;
        canvas.width =  size; canvas.height = size;
        gl.viewport(0, 0, size, size);
        draw_scene();
    }

    webXtri.canvas.onmousedown = function ( ev ){
        drag  = 1;

        if(ev.offsetX) {
            vOffs = [ev.offsetX, ev.offsetY];
        } else if(ev.layerX) {
            vOffs = [ev.layerX, ev.layerY];
        }
//         vOffs = [ev.clientX, ev.clientY];
        vStart = vOffs;
        vRot = vOffs;
        quat4.set(q_rotate, q_start);
        return false;
    }

    webXtri.canvas.onmouseup = function ( ev ){
        drag  = 0;
//         vOffs = [ev.clientX, ev.clientY];
        if(ev.offsetX) {
            vOffs = [ev.offsetX, ev.offsetY];
        } else if(ev.layerX) {
            vOffs = [ev.layerX, ev.layerY];
        }
    }

    document.onmouseup = function( ev ) {
        drag = 0; // stop dragging if out of the window
    }

    webXtri.canvas.onmousemove = function ( ev ){
        if(drag == 0)
            return;

        var xx = 0, yy = 0;
        if(ev.offsetX) {
            xx = ev.offsetX, yy = ev.offsetY;
        } else if(ev.layerX) {
            xx = ev.layerX, yy = ev.layerY;
        }
            
        if( ev.shiftKey) {
            translate[2] *= 1 + (yy - vOffs[1])/1000;
            yRot = - vOffs[0] + xx;
        } else if ( ev.ctrlKey ) {
            translate[0] += (xx - vOffs[0]) / 100;
            translate[1] += (vOffs[1] - yy) / 100;
        } else {
            vRot = [xx, yy];
/*            yRot = - vOffs[0] + xx;
            xRot = - vOffs[1] + yy;*/
        }
        vOffs = [xx, yy];
        draw_scene();
    }

    var wheelHandler = function(ev) {
        var del = 1.1;
        if (ev.shiftKey) del = 1.01;
        var ds = ((ev.detail || ev.wheelDelta) > 0) ? del : (1 / del);
        translate[2] *= ds;
        draw_scene();
        ev.preventDefault();
    };
    webXtri.canvas.addEventListener('DOMMouseScroll', wheelHandler, false);
    webXtri.canvas.addEventListener('mousewheel', wheelHandler, false);
}

webXtri.compile_shader = function(gl, s_type, source) {

    var shader;
    shader = gl.createShader(s_type);
    gl.shaderSource(shader, source);
    gl.compileShader(shader);
    if(gl.getShaderParameter(shader, gl.COMPILE_STATUS) == 0) {
        alert(gl.getShaderInfoLog(shader));
    }
    return shader;
}

webXtri.create_shader_program = function(gl, vert_src, frag_src) {

    var prog = gl.createProgram();
    gl.attachShader(prog, webXtri.compile_shader(gl, gl.VERTEX_SHADER,
                                  vert_src));
    gl.attachShader(prog, webXtri.compile_shader(gl, gl.FRAGMENT_SHADER,
                                  frag_src));
    gl.linkProgram(prog);
    return prog;
}

webXtri.external_shader_prog = function(gl, vert_name, frag_name) {

    return webXtri.create_shader_program(gl,
              webXtri.shader_sources[vert_name].str,
              webXtri.shader_sources[frag_name].str);
}

webXtri.embedded_shader_prog = function(gl, vert_id, frag_id) {
    
    var vert_script = document.getElementById(vert_id);
    var frag_script = document.getElementById(frag_id);
    var vert_src = "", frag_src = "";
    var k = vert_script.firstChild;
    while(k) {
        if(k.nodeType == 3)
            vert_src += k.textContent;
        k = k.nextSibling;
    }

    k = frag_script.firstChild;
    while(k) {
        if(k.nodeType == 3)
            frag_src += k.textContent;
        k = k.nextSibling;
    }

    return webXtri.create_shader_program(gl, vert_src, frag_src);
}

webXtri.add_external_shader = function(name, filename) {
    webXtri.shader_sources[name] = {
         str: null,
         loaded: false,
         filename: filename
    }
    webXtri.n_shaders_total +=1;
    return webXtri.shader_sources[name];
}

webXtri.load_shader = function(s) {
    var xhr = ajaxRequest();
    xhr.open("GET", webXtri.shader_sources[s].filename, true)
    xhr.onload = function () {
        webXtri.shader_sources[s].loaded = true;
        webXtri.shader_sources[s].str = xhr.responseText;
        webXtri.n_shaders_loaded += 1;
        if(webXtri.n_shaders_total == webXtri.n_shaders_loaded) {
            webXtri.onShadersLoaded();
        }
    }
    xhr.send();
}

webXtri.load_external_shaders = function(callback) {
    if(callback) {
        webXtri.onShadersLoaded = callback;
    }
    for(var s in webXtri.shader_sources) {
        webXtri.load_shader(s)
    }
}

function set_z_patches(n_z_patches, data) {

    var obj = document.getElementById("z-patches");
    obj.options.length = 0;
    z_patches.lenght = 0;
    z_patches[0] = 0; // put terminating symbol
    for(var i = 1; i <= n_z_patches; i++) {
        obj.options[obj.options.length] =
            new Option("patch " + i, i);
        z_patches[i] = data[i - 1];
    }
    obj.options[obj.options.length] =
            new Option("< unselect all >", n_z_patches + 1);
}

function select_z_patch(idx) {
    if(idx == "")
        return false;

//     if(idx == n_z_patches + 1)
//         alert("remove all");

    draw_scene();
    return true;
}

function is_any_z_patch_selected() {
    var obj = document.getElementById("z-patches");
    for(var i = 0; i < n_z_patches; i++) {
        if(obj.options[i].selected) {
            return true;
        }
    }
    return false;
}

function set_absolute_bounds() {

    if(document.getElementById("use_abs").checked)
        return;

    webXtri.left = Math.round(webXtri.left * MATH_ROUND) / MATH_ROUND;
    webXtri.right = Math.round(webXtri.right * MATH_ROUND) / MATH_ROUND;
    webXtri.btm = Math.round(webXtri.btm * MATH_ROUND) / MATH_ROUND;
    webXtri.top = Math.round(webXtri.top * MATH_ROUND) / MATH_ROUND;

    document.getElementById("left-abs-id").value = webXtri.left;
    document.getElementById("right-abs-id").value = webXtri.right;
    document.getElementById("btm-abs-id").value = webXtri.btm;
    document.getElementById("top-abs-id").value = webXtri.top;
}

function debug_msg(msg) {
    document.getElementById("progress_msg").innerHTML = msg;
}
   

function report_error(error_id) {

    btri_disabled = true;
    var txt;
    switch(error_id) {
    case ERR_OK:
        txt = "";
        btri_disabled = false; // enable recompute mesh button
        break;
    case ERR_INVALID_POLYNOMIAL:
        txt = "Incorrect format of the polynomial!";
        break;
    case ERR_INVALID_DATA:
        txt = "Invalid client parameters!";
        break;
    case ERR_TIMEOUT:
        txt = "Computation is timed out!";
        break;
    case ERR_INVALID_REFERRER:
        txt = "Invalid referrer!";
        break;
    case ERR_SERVER_TIMEOUT:
        txt = "No connection to the server!";
        break;
    case ERR_REQUEST_PENDING:
        txt = "You request is being processed!";
        break;
    case ERR_SERVER_OVERLOAD:
        txt = "Server is overloaded!";
        break;
    case ERR_TRIANGULATE_GENERIC:
        txt = "Triangulation algorithm error!";
        break;
    case ERR_INVALID_MESH:
        txt = "Unable to parse surface model!";
        break;
    default:
        txt = "Generic error!";
    }
    document.getElementById("banalyse").disabled = false;
    document.getElementById("btriangulate").disabled = btri_disabled;
    document.getElementById("progress_msg").innerHTML = "";
    if(error_id == ERR_OK)
        document.getElementById("error_msg").innerHTML = "";
    else
        document.getElementById("error_msg").innerHTML = "ERROR: " + txt;
}

var ping_mode;

function send_server_request(inURL, params) {

    document.getElementById("banalyse").disabled = true;
    document.getElementById("btriangulate").disabled = true;
    document.getElementById("error_msg").innerHTML = "";
    document.getElementById("progress_msg").innerHTML = "";

    http = ajaxRequest();
    var getURL = "../cgi-bin/xtri.cgi?id=4&sid=" + server_id;
    http.open("GET", getURL, true)
    http.send(null) //send GET request

    ping_mode = true;
    http.onreadystatechange = function() {

/*         document.getElementById("progress_msg").innerHTML = "response received.. " + http.readyState + ' ' + ping_mode;*/
        
    if(http.readyState == 4) {
 //           report_error(0); // unknown error occurred

// TODO: need to make sure that client further works with same server

    if(ping_mode) {
        var vars = http.responseText.slice(
             http.responseText.indexOf('?') + 1).split('&');

        var i, parsed = new Array;
        for(i = 0; i < vars.length; i++)  {
            var s = vars[i].split('=');
            parsed[s[0]] = s[1];
        }
        server_id = parseInt(parsed['sid']);
        if(isNaN(server_id))
            server_id = 0;
        params += "&sid=" + server_id

        var err_id = parseInt(parsed['error']);
        if(isNaN(err_id) || err_id != ERR_OK) {
            report_error(err_id);
            return false;
        }
        ping_mode = false;
        document.getElementById("progress_msg").innerHTML =
            "Computing, please wait..";

        http.open("POST", inURL, true);
        http.setRequestHeader("Content-type", "application/x-www-form-urlencoded");
        http.setRequestHeader("Content-length", params.length);
        http.setRequestHeader("Connection", "close");
        http.send(params);

    } else {

        document.getElementById("progress_msg").innerHTML = "Model received, parsing..";
        parse_model(http.responseText);

    }
    } // readyState == 4
    } // onreadystatechange
}

// type: 0 - integer; 1 - float
function check_val(id, minc, maxc, type) {
    var obj = document.getElementById(id), val;
    if(type == 0)
        val = parseInt(obj.value);
    else
        val = parseFloat(obj.value);

    if(isNaN(val) || val < minc || val > maxc) {
        obj.style.border = "2px solid #FF0000";
        return Number.NaN;
    } else {
        obj.style.border = "1px solid #7F9DB9";
        return val;
    }
}

function parse_parameters() {

    var params = "";
    var use_abs = document.getElementById("use_abs").checked;

    var left, right, btm, top;
    if(!use_abs) {
        left = check_val("left-id", 0, 200, 1) / 100.0;
        right = check_val("right-id", 0, 200, 1) / 100.0;
        btm = check_val("btm-id", 0, 200, 1) / 100.0;
        top = check_val("top-id", 0, 200, 1) / 100.0;
    } else {
        left = check_val("left-abs-id", -50, 50, 1);
        right = check_val("right-abs-id", -50, 50, 1);
        btm = check_val("btm-abs-id", -50, 50, 1);
        top = check_val("top-abs-id", -50, 50, 1);
    }

    var below = check_val("below-id", -50, 50, 1);
    var above = check_val("above-id", -50, 50, 1);
    var sx = check_val("sx-id", 1, 50, 0);
    var sy = check_val("sy-id", 1, 50, 0);

    if(isNaN(left) || isNaN(right) || isNaN(btm) || isNaN(top) ||
        isNaN(below) || isNaN(above) || isNaN(sx) || isNaN(sy) ||
            below >= above) {
        document.getElementById("error_msg").innerHTML=
                "Please correct the input parameters!";
        return params;
    }
    if(use_abs && (left >= right || btm >= top)) {
        document.getElementById("error_msg").innerHTML=
                "Please correct the input parameters!";
        return params;
    }

    // id = 1: surface=ascii_poly
    // id = 2: surfaceid=poly_id
    if(!use_abs)
        params = "r=0";
    else
        params = "r=1";

    params += "&sx=" + sx + "&sy=" + sy +
        "&left=" + left + "&right=" + right +
        "&btm=" + btm + "&top=" + top +
        "&above=" + above+ "&below=" + below;

    return params;
}

function click_analyse() {

    var poly = document.getElementById("surf").value;
    poly = encodeURIComponent(poly);
    if(poly == "")
        return false;

    var inURL = "../cgi-bin/xtri.cgi"
    var params = parse_parameters();
    if(params == "")
        return false;

    params += "&id=1&surface="+poly;
    send_server_request(inURL, params);
    return true;
}

function click_triangulate() {

    var inURL = "../cgi-bin/xtri.cgi"
    var params = parse_parameters();
    if(params == "")
        return false;

    params += "&id=2&surfaceid=" + surfaceID[0] + " " + surfaceID[1] + " " +
        surfaceID[2] + " " + surfaceID[3];

    send_server_request(inURL, params);
    return false;
}

function clear_text() {
    document.getElementById("surf").value = "";
}

function quaternion_from_axis(axis, angle) {

    var a = angle * Math.PI / 360.0;
    var cs = Math.sin(a);
    var v = axis;

    v[3] = Math.cos(a);
    var q = quat4.create(v);
    quat4.normalize(q);
    q[0] *= cs, q[1] *= cs, q[2] *= cs;
    return q;
//     alert(quat4.str(q));
//     q.w = cos(a);
//     q.v = axis;
//     q.copy_on_write();
//     q.v.normalize();
//     q.v *= sin(angle);
}

function quaternion_trackball(pp1, pp2) {

    var p1 = [pp1[0], pp1[1]],
        p2 = [pp2[0], pp2[1]];
    p1[0] = 2.0*p1[0] / webXtri.canvas.width - 1.0;
    p1[1] = 1.0 - 2.0*p1[1] / webXtri.canvas.height;

    p2[0] = 2.0*p2[0] / webXtri.canvas.width - 1.0;
    p2[1] = 1.0 - 2.0*p2[1] / webXtri.canvas.height;

/*    debug_msg("mv1: " + pp2[0] + ' ' + pp2[1] + "<br />" +
        p2[0] + ' ' + p2[1]);*/
    
    var v1 = vec3.create([p1[0], p1[1], 0]),
        v2 = vec3.create([p2[0], p2[1], 0]);

    var d = p1[0]*p1[0] + p1[1]*p1[1];
    if(d > 1.0) {
        d = Math.sqrt(d);
        vec3.scale(v1, 1.0/d);
    } else
        v1[2] = Math.sqrt(1.0 - d);

    d = p2[0]*p2[0] + p2[1]*p2[1];
    if(d > 1.0) {
        d = Math.sqrt(d);
        vec3.scale(v2, 1.0/d);
    } else
        v2[2] = Math.sqrt(1.0 - d);

    var dot = 1.0 + vec3.dot(v1, v2);
    vec3.cross(v1, v2);
    vec3.negate(v1);
    
    var q = quat4.create([v1[0],v1[1],v1[2],dot]);
    return q;
//     q.v = -v1.cross_product(v2);
//     q.w = 1.0 + v1.dot_product(v2);    
}

    /**
 c.addEventListener('webglcontextlost', handleContextLost, false);
 c.addEventListener('webglcontextrestored', handleContextRestored, false);

            framerate = new Framerate("framerate");
            var f = function() {
                drawPicture(ctx);
                requestId = window.requestAnimFrame(f, c);
            };
            f();

            function handleContextLost(e) {
                e.preventDefault();
                clearLoadingObjects();
                if (requestId !== undefined) {
                    window.cancelRequestAnimFrame(requestId);
                    requestId = undefined;
                }
            }

            function handleContextRestored() {
                loaded = false;
                init();
                f();
            }
        }
    **/
