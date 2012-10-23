
#ifdef GL_ES
precision highp float;
#endif

#define GOLD 255.0/255.0, 215.0/255.0, 0.0/255.0
#define EMERALD 0x50, 0xc8, 0x78
#define NICE_BLUE 40.0/255.0, 80.0/255.0, 200.0/255.0
#define VIOLET  200.0/255.0, 10.0/255.0, 80.0/255.0
#define ORANGE 253.0/255.0, 192.0/255.0, 24.0/255.0
#define GREEN 50.0/255.0, 200.0/255.0, 100.0/255.0
#define BLUE 40.0/255.0, 100.0/255.0, 230.0/255.0

// uniform float rescale_normals;
// uniform bool bumpy_surf;

attribute vec3 aPos;
attribute vec3 aNormal;

uniform mat4 normal_matrix;
uniform mat4 mview_matrix;
uniform mat4 proj_matrix;

//    varying mediump vec4 vColor;
varying vec3 vPos, vNormal;
//    varying vec3 vPos_original;

void main(void) {

    vec3 nn = aNormal;
    bool bumpy_surf = true;

/*    if(bumpy_surf) {
        nn += vec3(sin(aPos.x*40.0), sin(aPos.x*40.0), sin(aPos.x*40.0))*0.2;
        normalize(nn);
    }*/
    vec4 pos = mview_matrix * vec4(aPos, 1.0);
    vPos = vec3(pos); // pos.xyz
     // is this correct ? or should we take 3x3 matrix ?
    vNormal = (vec3(normal_matrix * vec4(nn, 0.0)));
    gl_Position = proj_matrix * pos;
}
   