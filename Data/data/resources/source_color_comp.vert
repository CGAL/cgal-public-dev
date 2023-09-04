varying highp vec4 vertex;
varying highp vec3 normal;
varying highp vec3 color;

uniform highp mat4 mvp_matrix;
uniform highp mat4 mv_matrix;
uniform highp float point_size;

varying highp vec4 fP;
varying highp vec3 fN;
varying highp vec4 fColor;

void main(void)
{
  fP = mv_matrix * vertex;
  highp mat3 mv_matrix_3;
  mv_matrix_3[0] = mv_matrix[0].xyz;
  mv_matrix_3[1] = mv_matrix[1].xyz;
  mv_matrix_3[2] = mv_matrix[2].xyz;
  fN = mv_matrix_3* normal;
  fColor = vec4(color, 1.0);
  gl_PointSize = point_size;

  gl_Position = mvp_matrix * vertex;
}