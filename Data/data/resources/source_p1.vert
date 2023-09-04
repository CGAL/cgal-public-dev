#version 150
in highp vec4 vertex;
in highp vec3 color;

uniform highp mat4 mvp_matrix;
uniform highp float point_size;

out highp vec4 fColor;
out highp vec4 m_vertex;

void main(void)
{
  gl_PointSize = point_size;
  fColor = vec4(color, 1.0);
  m_vertex = vertex;
  gl_Position = mvp_matrix * vertex;
}