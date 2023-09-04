#version 150
in highp vec4 vertex;
in highp vec3 normal;
in highp vec3 color;

uniform highp mat4 mvp_matrix;
uniform highp mat4 mv_matrix;
uniform highp float point_size;

out highp vec4 fP;
out highp vec3 fN;
out highp vec4 fColor;
out highp vec4 m_vertex;

void main(void)
{
  fP = mv_matrix * vertex;
  fN = mat3(mv_matrix)* normal;
  fColor = vec4(color, 1.0);
  gl_PointSize = point_size;

  m_vertex = vertex;

  gl_Position = mvp_matrix * vertex;
}