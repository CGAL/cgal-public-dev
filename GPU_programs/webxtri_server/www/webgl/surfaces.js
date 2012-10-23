// this file contains a list of algebriac curve equations in the order they appear in gallery

function init_surfaces()
{
//     var pics = document.getElementById("picture_div");
//     document.getElementById("main_content").style.left = 
//         (pics.offsetLeft+pics.offsetWidth)+'px';
    //alert((pics.offsetLeft+pics.offsetWidth));
    for(var key in surfaces) {
        var obj = document.getElementById(key);
        if(obj != null)
            obj.value = surfaces[key];
    }
}

function submit_form(id)
{
    var val = document.getElementById(id).value;
    val = val.replace(/\s/gi, "");
    // encode '+' character otherwise it will be treated as space
    val = val.replace(/\+/gi, "\%2B");
    document.getElementById("surf").value = val;
    document.forms[0].submit(); 
}

//submit_form(\''+id+'\')
function set_surface(id, img, caption, degree)
{
    var str = '<tr><td align="center" valign="top">'+caption+
        '<br />degree:  <b>'+degree+'</b><br /><br /></td>'+
        '<td align="center" valign="top"><img class="image_btn"  onclick="submit_form(\''+id+'\')"  src="' +
            img+'" height="300" /></td><td align="center" valign="top"><textarea cols="45" style="border: 1px solid #7F9DB9;" rows="6" wrap="soft" readonly="readonly" id="'+
            id+'"></textarea></td></tr><tr><td colspan="3"><br /><hr /></td></tr>';

    document.write(str);
}

var surfaces = {
    "clebsch_diag_surf" : "x^3-3*x^2*y-3*x*y^2+y^3+1.5*x^2*z+1.5*y^2*z-1.9375*z^3-5.845671475544960865655131402583*z^2-4.5*z",
    "tangle_cube" : "z^4 + (-5)*z^2 + (y^4 + (-5)*y^2 + (x^4 + (-5)*x^2 + 10))",
    "KM43" : " z^3 + (((-3)*x + 3)*y^2 + (x^3 + 3*x^2 + (-4)))",
    "ffo" : "z + ((5*x)*y^4 + ((-10)*x^3)*y^2 + (x^5))",
    "vis-a-vis": "(-1)*z^4 + z^3 + (y^4 + y^2 + ((-1)*x^3 + x^2))"
 };

//  "apple3D" : "z^2*20+y^7 + (-3)*y^6 + (2*x^2 + (-1)*x + 2)*y^5 + (x^3 + (-6)*x^2 + x + 2)*y^4 + (x^4 + (-2)*x^3 + 2*x^2 + x + (-3))*y^3 + (2*x^5 + (-3)*x^4 + x^3 + 10*x^2 + (-1)*x + 1)*y^2 + ((-1)*x^5 + 3*x^4 + 4*x^3 + (-12)*x^2)*y + (x^7 + (-3)*x^5 + (-1)*x^4 + (-4)*x^3 + 4*x^2)",