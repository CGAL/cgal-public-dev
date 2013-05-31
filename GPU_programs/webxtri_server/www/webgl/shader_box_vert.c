
attribute vec3 aPos;
// attribute vec4 aVertexColor;
attribute vec2 aTexCoord;

uniform mat4 normal_matrix;
uniform mat4 mview_matrix;
uniform mat4 proj_matrix;

varying vec2 vTexCoord;
// varying vec4 vColor;

void main(void) {
    gl_Position = proj_matrix * mview_matrix * vec4(aPos, 1.0);
//     vColor = aVertexColor;
    vTexCoord = aTexCoord;
}