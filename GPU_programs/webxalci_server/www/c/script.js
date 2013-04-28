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
// Author(s)     : Pavel Emeliyanenko (www.mpi-inf.mpg.de/~emeliyan)
//
// ============================================================================

var browser;

var check_group = new Array;
var n_check_btn = 0;

var TEXT_BTN_UP   = "#AABBDD";
var TEXT_BTN_OVER = "#AAAADD";
var TEXT_BTN_DOWN = "#7788CC";
var IMG_BTN_UP    = "#FFFFFF";
var IMG_BTN_DOWN  = "#CCCCCC";
var DISABLED_BGCOLOR = "#CCCCDD";

var TOOLTIP_DELAY = 500;
var TOOLTIP_BGCOLOR = "#E4E7FF";
var TOOLTIP_BORDER = "#002299";

function whatbrowser()
{
    if(document.layers)        
        browser="NN4";    
    if(document.all)        
        browser="ie";    
    if(!document.all && document.getElementById)        
        browser="NN6";
}

/*function xalci_DoFSCommand(cmd, args) {
 if(cmd == "messagebox") {
 //alert(args);
 //document.write("<a href="+args+">
 window.open(args, "Save plot");
 }
}*/

function init(){

    init_mini();
}

function init_mini()
{
    var pics = document.getElementById("picture_div");
    var obj = document.getElementById("main_content");
    if(obj != undefined)
        obj.style.left = (pics.offsetLeft+pics.offsetWidth+60)+'px';
    whatbrowser();  
}

// shows a menu with 'index' item selected and relative y-position: 'y_top'
function show_menu(index, y_top)
{
    var ids = ["hp-home", "hp-curve-gallery", "hp-surf-gallery", "hp-xalci", "hp-xtri",
        "hp-feedback", "hp-usage"];
    var paths = ["index.html", "gallery.html", "webgl/surface_gallery.html", "cgi-bin/xalci.cgi", 
        "webgl/", "feedback.html", "usage.html"];

    ids[index] = "hp-index";
    
    var rel = "../";
    if(index >= 2 && index <= 4) {
        for(var i = 0; i < paths.length; i++) {
            paths[i] = rel + paths[i];
        }
        
        if(index == 2)
           paths[2] = "#";
        else if(index == 3)
           paths[3] = "xalci.cgi";
        else 
           paths[4] = "#";
    }

    var text = '<div style="position: relative; top: '+y_top+'; width: 15em; z-index:4;"><div id="col1o2content"><ul class="linklist">' + 
'<li><a href="'+paths[0]+'" accesskey="h" title="XAlci web-demo home" id="'+ids[0]+'"><span class="ak">H</span>ome</a></li>' +
'<li><a href="'+paths[1]+'" accesskey="c" title="Gallery of implicit algebraic curves and their arrangements" id="'+ids[1]+'">Gallery of algebraic <span class="ak">c</span>urves</a></li>' +
'<li><a href="'+paths[2]+'" accesskey="s" title="Gallery of implicit algebraic surfaces" id="'+ids[2]+'">Gallery of algebraic <span class="ak">s</span>urfaces <span style="color: #FF2222;">(NEW!)</span></a></li>' +
'<li><a href="'+paths[3]+'" accesskey="a" title="XAlci web-demo" id="'+ids[3]+'">X<span class="ak">A</span>lci web-demo</a></li>' +
'<li><a href="'+paths[4]+'" accesskey="t" title="XTri web-demo" id="'+ids[4]+'">X<span class="ak">T</span>ri web-demo <span style="color: #FF2222;">(NEW!)</span> </a> </li>' +
'<li><a href="'+paths[5]+'" accesskey="f" title="Send us your anonynous comments and suggestions how to improve the program" id="'+ids[5]+'">Send <span class="ak">f</span>eedback</a></li>' +
'<li><a href="'+paths[6]+'" accesskey="u" title="Using the program" id="'+ids[6]+'"><span class="ak">U</span>sing the program</a></li>';
    
    /*if(index >= 3) // insert submenu
        text += '<ul style="list-style:none; margin-left:1.5em; padding-left:0;"><li><a href="'+paths[4]+'" accesskey="u" title="Using  interactive server" id="'+ids[4]+'"><span class="ak">U</span>sing interactive server</a></li>'+
'<li><a href="'+paths[5]+'" accesskey="o" title="Describes arrangement computation and visualization algorithms" id="'+ids[5]+'"><span class="ak">O</span>ur algorithms</a></li></ul>';*/
    
    text += '</ul></div></div>';
    document.write(text);
}

// show the authors layer with relative y-position 'y_top'
function show_authors(y_top)
{
//     document.write('<div style="position: relative; top: '+y_top+'; background-color: #E4E7FF; border: solid 1px #002299; padding: 5px;">'+
// '<span style="font-weight: bold;">Authors:</span><br /><br />' +
// '<a href="http://www.mpi-inf.mpg.de/~emeliyan">Pavel Emeliyanenko</a>:<br />'+
// 'Graphing implicit algebraic curves<br />and web-server design.<br /><br />' +
// '<a href="http://www.mpi-inf.mpg.de/~mkerber">Michael Kerber</a>: <br />' +
// 'Analysis of real algebraic curves <br />(<span style="font-variant: small-caps;">AlciX</span> library).<br /></div> ');

    document.write('<div style="position: relative; top: '+y_top+'; background-color: #E4E7FF; border: solid 1px #002299; padding: 5px;">'+
'<br />Web-application design<br /> and server software<br /> by ' +
'<a href="http://www.mpi-inf.mpg.de/~emeliyan">Pavel Emeliyanenko</a>:<br /><br />'+
'For relevant algorithms,<br /> see <a href="#">Web-demo home</a>.<br /> &nbsp;</div>');
}

/*
function start_drag(evt)
{
    if(is_focus_point||is_busy)
        return true;
    if(browser=="ie") {        
        event.cancelBubble = true;
        x = event.clientX + document.documentElement.scrollLeft; 
        y = event.clientY + document.documentElement.scrollTop;    
    } else if(browser=="NN6") {        
        evt.cancelBubble = true;
        x = evt.pageX;        
        y = evt.pageY;    
    }    
    x -= plot_x;
    y -= plot_y;
    if(x < 0 || x > plot_w || y < 0 || y > plot_h)
        return true;    
    
    is_drag = true;
    is_moved = false;
    
    var zoomer=document.getElementById("zoomer");    
    start_x = x + plot_x;
    start_y = y + plot_y;
    zoomer.style.left = (start_x - parent_x) + 'px';
    zoomer.style.top = (start_y - parent_y) + 'px';
    zoomer.style.width = '0px';
    zoomer.style.height = '0px';
    zoomer.style.visibility = "visible";
          
    if(browser=="ie") { 
        if(document.selection && document.selection.empty) 
            document.selection.empty();
        document.documentElement.onselectstart = function() {
            event.returnValue = false;
        };
    }
    return false;
}

function stop_drag(evt)
{
    if(is_focus_point||!is_drag||is_busy)
        return false;
    is_drag = false;
    if(!is_moved)
        return false;
    var zoomer=document.getElementById("zoomer");   
    zoomer.style.visibility = "hidden";
    if(browser=="ie") {        
        x = event.clientX + document.documentElement.scrollLeft; 
        y = event.clientY + document.documentElement.scrollTop;    
    } else if(browser=="NN6") {        
        x = evt.pageX;        
        y = evt.pageY;    
    }    
    x -= plot_x;
    y = plot_h - y + plot_y;
    start_x -= plot_x;
    start_y = plot_h - start_y + plot_y;
    if(start_x > x) { 
        tmp = start_x;
        start_x = x;
        x = tmp;
    }
    if(start_y > y) { 
        tmp = start_y;
        start_y = y;
        y = tmp;
    }
    if(x - start_x <= 3||y - start_y <= 3)
        return false;
    lx = (x_max - x_min)/plot_w;
    ly = (y_max - y_min)/plot_h;
    x_min_ = x_min + start_x*lx; 
    y_min_ = y_min + start_y*ly;
    x_max_ = x_min + x*lx; 
    y_max_ = y_min + y*ly;
    
    if(x_max_ - x_min_ < window_min||y_max_ - y_min_ < window_min) {    
        alert("Incorrect program parameters!");     
        return false;   
    }   
    xc = (x_min_ + x_max_) / 2;
    yc = (y_min_ + y_max_) / 2;
    lx = (x_max_ - x_min_) / 2;
    ly = (y_max_ - y_min_) / 2;
    x_min = xc - lx;    
    x_max = xc + lx;    
    y_min = yc - ly;    
    y_max = yc + ly;
    return check_plot();
}

function track_mouse(evt){    
    
    var x, y;    
    if(browser=="ie") {        
        x = event.clientX + document.documentElement.scrollLeft; 
        y = event.clientY + document.documentElement.scrollTop;    
    } else if(browser=="NN6") {        
        x = evt.pageX;        
        y = evt.pageY;    
    }    
    x -= plot_x;
    y -= plot_y;
    if(x < 0 || x > plot_w || y < 0 || y > plot_h)
        return false;    
   
    if(is_drag) {
        var zoomer=document.getElementById("zoomer");    
        w = (x + plot_x - start_x);
        h = (y + plot_y - start_y);
        if(w < 0) {
            w = -w;
            zoomer.style.left = (x + plot_x - parent_x) + 'px';
        }
        if(h < 0) {
            h = -h;
            zoomer.style.top = (y + plot_y - parent_y) + 'px';
        }
        zoomer.style.width = w + 'px';
        zoomer.style.height = h + 'px';
        is_moved = true;
    }
    y = plot_h - y;
    x = x_min + x*(x_max-x_min)/plot_w;    
    y = y_min + y*(y_max-y_min)/plot_h;    
    var coords=document.getElementById("coords");        
    coords.innerHTML="Coords: (" + x.toPrecision(7) + ", " + 
        y.toPrecision(7) + ")";   
    return false;
}

function focus_on_point(evt) {
    var x, y;
    if(is_busy||!is_focus_point)
        return false;
    if(browser=="ie") {
        x = event.offsetX;
        y = plot_h - event.offsetY;
    } else if(browser=="NN6") {
        x = evt.layerX;
        y = plot_h - evt.layerY;
    }
    var xc = x_min + x*(x_max-x_min)/plot_w; 
    var yc = y_min + y*(y_max-y_min)/plot_h; 
    var lx = (x_max-x_min)/2;
    var ly = (y_max-y_min)/2;
    x_min = xc - lx;
    x_max = xc + lx;
    y_min = yc - ly;
    y_max = yc + ly;
    return check_plot();
}

function check_plot() {    
    if(is_busy)
        return false;
    
    var src, rend_type = document.getElementById("renderer").selectedIndex;
    if(rend_type == 1) {
        if(!check_empty())
            return false;
        val = document.getElementById("cv").value;
        val = val.replace(/\s/gi, "");
        // encode '+' character otherwise it will be treated as space
        val = val.replace(/\+/gi, "\%2B"); 
        src = 'xalci.cgi?s=3&curve=' + val;
    } else {
        var idx = document.getElementById("idx");    
        src='xalci.cgi?s=2&curveid='+document.getElementById("curveid").value;
        if(document.getElementById("checkbox_all").checked)         
            src += '&all=1';     
        else if(document.getElementById("checkbox_flat").checked)
            src += '&all=2';
        else        
            for(i=0; i < idx.options.length; i++)           
                if(idx.options[i].selected)             
                    src += '&idx='+idx.options[i].value;     
    }
    src += '&x_min='+x_min+'&x_max='+x_max+'&y_min='+y_min+'&y_max='+y_max;  
    disable_controls(true);    
    document.getElementById("plot").src = src;    
    return false;
}

function zoom(ratio) {  
    var xc = (x_min+x_max)/2;   
    var yc = (y_min+y_max)/2;   
    var lx = (x_max-x_min)*ratio/2; 
    var ly = (y_max-y_min)*ratio/2; 
    x_min = xc - lx;    
    x_max = xc + lx;    
    y_min = yc - ly;    
    y_max = yc + ly;    
    return true;
}

function reset_view() { 
    if(is_busy)
        return false;
    x_min = -2.0;   
    x_max = 2.0;    
    y_min = -1.5;   
    y_max = 1.5;    
    return check_plot();
}

function zoomin() {
    
    if(is_busy)
        return false;
    if(x_max - x_min < window_min||y_max - y_min < window_min) {    
        alert("Incorrect program parameters!");     
        return false;   
    }   
    zoom(0.5);  
    return check_plot();
}

function zoomout() {
    if(is_busy)
        return false;
    if(x_max - x_min > window_max||y_max - y_min > window_max) {    
        alert("Incorrect program parameters!");     
        return false;   
    }   
    zoom(2);    
    return check_plot();
}

function toggle_all(obj)
{
    if(obj.checked)
        document.getElementById("checkbox_flat").checked = 0;
}

function toggle_flat(obj)
{
    if(obj.checked)
        document.getElementById("checkbox_all").checked = 0;
}

function renderer_changed(obj)
{
    var flag = false;
    if(obj.selectedIndex == 1) 
        flag = true;
    
    if(obj.selectedIndex == 1) {   
     // enable controls at the beginning when subdiv methods is  hosen
        is_busy = false;
        document.getElementById("plot").onload = 
            function() { disable_controls(false); }
        disable_controls(false);
        document.getElementById("checkbox_flat").disabled = true;
        document.getElementById("checkbox_all").disabled = true;        
        document.getElementById("idx").disabled = true;
    } else if(obj.selectedIndex == 0) {
      // in case Segment renderer is selected, enable rasterize button
      // if the number of segments is non-zero, i.e. curve was analysed
        if(document.getElementById("idx").options.length > 0) {
            document.getElementById("idx").disabled = false;
            document.getElementById("brasterize").disabled = false;
            document.getElementById("checkbox_flat").disabled = false;
            document.getElementById("checkbox_all").disabled = false;        
        } else { // nothing was analysed: disable controls
            disable_controls(true);
            is_busy = true;
        }
    }
    
    document.getElementById("banalyse").disabled = flag;
    document.getElementById("renderer").disabled = false;
}

function disable_controls(flag) 
{
    is_busy = flag;  
    document.getElementById("bzoomin").disabled = flag;  
    document.getElementById("bzoomout").disabled = flag;     
    document.getElementById("breset").disabled = flag;   
    document.getElementById("bpoint").disabled = flag;   
    document.getElementById("bregion").disabled = flag;  
    document.getElementById("brasterize").disabled = flag;   
    document.getElementById("banalyse").disabled = flag;
    document.getElementById("renderer").disabled = flag;
    
    //alert('dis controls');
}

function comment_change(obj)
{
    document.getElementById("n_symbols").innerHTML = 
        obj.value.length;
    return true;
}

function submit_comment(n_chars)
{
    var txt = document.getElementById("comment").value;
    if(txt.length == 0) {
        alert("Please enter a comment!");
        return false;
    }
    if(txt.length > n_chars) {
        alert("Your comment exceeds "+n_chars+" characters!");
        return false;
    }
    document.getElementById("id").value=3;
    document.forms[0].submit();    
    return true;
}
*/

// checked: 0 - normal button, 1 - unchecked, 2 - checked
function set_img_button(id, img_src, click, alt, checked, disabled)
{
    var btn_html = '<button id="'+id+'" class="img_btn" ';
    if(disabled)
       btn_html += 'disabled ';
//        'style="background-color: '+DISABLED_BGCOLOR+';" onclick="return false;" disabled>';
//     else
    btn_html += 'onmousedown="mouse_click(this,0);" onmouseup="mouse_click(this,1)" onmouseout="mouse_over_out(event,this,0)" onmouseover="mouse_over_out(event,this,1)" '+
       ' onclick="'+click+'">';
    btn_html += '<img src="'+img_src+'" class="inner_img" /></button>';
//     onmousemove="move_tip(event,this)"
//     btn_html += '<div id="'+id+'MyTiP" align="left" style="position: absolute; left: 0px; top: 0px; background-color: '+TOOLTIP_BGCOLOR+'; border: solid 1px '+TOOLTIP_BORDER+'; visibility: hidden; z-index: 1000; overflow: hidden;">'+alt+'</div>';
    document.write(btn_html);
//     var obj = document.getElementById(id);
//     obj.tip = document.getElementById(id+"MyTiP");
    if(checked == 0)
        return;
    check_group[n_check_btn++] = obj;
    obj.checked = checked - 1;
    if(obj.checked) {
        mouse_click(obj, 0);
        obj.style.borderStyle = "inset";
    }
}

function set_text_button(id, txt, click, disabled)
{
    var btn_html = '<button id="'+id+'" class="text_btn" ';
    if(disabled) {
//        btn_html += 'style="background-color: '+DISABLED_BGCOLOR+';" onclick="return false;" disabled><span style="color: #999999;">'+txt+'</span></button>';   
       btn_html += 'disabled ';   
    } 
    btn_html += 'onmousedown="mouse_click(this,2);" onmouseup="mouse_click(this,3)" onmouseout="mouse_over_out(event,this,2)" onmouseover="mouse_over_out(event,this,3)" onclick="'+click+'">'+txt+'</button>';
    document.write(btn_html);
}


function move_tip(evt, obj)
{
    var x, y;
    if(browser=="ie") {        
        x = event.clientX + document.documentElement.scrollLeft; 
        y = event.clientY + document.documentElement.scrollTop;    
    } else if(browser=="NN6") {        
        x = evt.pageX;        
        y = evt.pageY;    
    }    
    obj.tip.style.left = (x - obj.tip.offsetParent.offsetLeft + 15)+ 'px';
    obj.tip.style.top = (y - obj.tip.offsetParent.offsetTop + 8)+ 'px';
}

function mouse_over_out(evt, obj, mode)
{   
    var is_text = mode>>1;
    mode &= 1;
    if(obj.checked == undefined||obj.checked == 0) {
//         obj.style.borderStyle = (mode ? "inset" : "outset");
        if(mode == 0) {
            obj.style.backgroundColor = is_text ? TEXT_BTN_UP:
            IMG_BTN_UP;
            obj.style.borderStyle = "";
        }
        if(mode == 1 && is_text) {
            obj.style.borderStyle = "outset";
            obj.style.backgroundColor = TEXT_BTN_OVER;
        }
    }
   
//     if(obj.tip == undefined || obj.tip == "")
//         return false;
//     if(mode == 1) {
//         var x, y;
//         if(browser=="ie") {        
//             x = event.clientX + document.documentElement.scrollLeft; 
//             y = event.clientY + document.documentElement.scrollTop;    
//         } else if(browser=="NN6") {   
//             x = evt.pageX;        
//             y = evt.pageY;    
//         }    
//         obj.tip.style.left = (x - obj.tip.offsetParent.offsetLeft + 15)+ 'px';
//         obj.tip.style.top = (y - obj.tip.offsetParent.offsetTop + 8)+ 'px';
//         obj.tip.timer = window.setTimeout("show_id('"+obj.tip.id+"')", TOOLTIP_DELAY);
//     } else if(mode == 0) {
//         if(obj.tip.timer) {
//             window.clearTimeout(obj.tip.timer);
//             obj.tip.timer = 0;       
//          }
//          hide_id(obj.tip.id);
//     }
}

// mode: 0 - down, 1 - up
function mouse_click(obj, mode){    
    var is_text = mode>>1; 
    mode &= 1;
    if(is_text) {
       obj.style.backgroundColor = (mode ? TEXT_BTN_UP: TEXT_BTN_DOWN);
       if(mode == 0)  
           obj.style.borderStyle = "inset";
       else
           obj.style.borderStyle = "outset";
       return;
    }
    // image button
    if(mode == 0)
        obj.style.borderStyle = "inset";
    else
        obj.style.borderStyle = "outset";
    
    if(mode == 1 && obj.checked == 1)
        return; 
    obj.style.backgroundColor = (mode ? IMG_BTN_UP: IMG_BTN_DOWN);
    if(mode == 1 && obj.checked == 0) {
        obj.style.backgroundColor = IMG_BTN_DOWN;
        obj.checked = 1;
        for(var i = 0; i < n_check_btn; i++) {
            if(check_group[i] == obj)
                continue;
            check_group[i].checked = 0;
            check_group[i].style.backgroundColor = IMG_BTN_UP;
            check_group[i].style.borderStyle = "outset";
        }
    }
}

function show_id(id)
{   
    document.getElementById(id).style.visibility = "visible";
}

function hide_id(id)
{   
    document.getElementById(id).style.visibility = "hidden";
}
