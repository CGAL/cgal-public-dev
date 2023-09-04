#version 150
in highp vec4 vertex;

uniform highp mat4 vp_matrix;
uniform highp mat4 m_matrix;

void main(void)
{
  gl_Position = vp_matrix * m_matrix * vertex;
}