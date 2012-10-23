
#ifdef GL_ES
precision highp float;
#endif

uniform int lighting_mode;
uniform int inverse_normals;
uniform vec4 uColor;
uniform vec4 uBackColor;
varying vec3 vPos, vNormal;

//     varying vec3 vPos_original;

const vec3 light1_pos = vec3(0.7, 0.0, -1.0);
const vec3 light1_color = vec3(0.3, 1.0, 0.4);
const vec3 light2_pos = vec3(-2.0, 0.0, 0.0);
const vec3 light2_color = vec3(1.0, 1.0, 1.0);
const vec3 ambient = vec3(0.1, 0.1, 0.1);

void main(void) {

    vec3 N = vNormal;

/*    N += vec3(sin(vPos_original.z*30.0), sin(vPos_original.z*30.0), cos(vPos_original.z*30.0))*0.1;*/

    vec3 color = uColor.rgb;
    if(!gl_FrontFacing) {
        N = -N;
        color = uBackColor.rgb;
    }
    if(inverse_normals == 1) {
        N = -N;
    }

    if(lighting_mode == 0) {
    // this is a position of light
    vec3 L = normalize(light1_pos - vPos);
    vec3 E = normalize(-vPos);
    vec3 R = normalize(-reflect(L, N));

    float ndotl = max(dot(N, L), 0.0);
    color = (color * ndotl);

    float ndoth = max(dot(R, E), 0.0);
    ndoth = pow(ndoth, 64.0);
    color += light1_color * ndoth;

    L = normalize(light2_pos - vPos);
    E = normalize(-vPos);
    R = normalize(-reflect(L, N));

    ndotl = max(dot(N, L), 0.0);
    color = (color * ndotl);

    ndoth = max(dot(R, E), 0.0);
    ndoth = pow(ndoth, 64.0);
    color += light2_color * ndoth;

    } else if(lighting_mode == 1) {
        color = vec3(100.0/255.0,100.0/255.0, 80.0/255.0);
    } else
        color = vec3(0.6, 0.8, 0.6);

    color += ambient;
    gl_FragColor = vec4(color, 1.0);
}
    