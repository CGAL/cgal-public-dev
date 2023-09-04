#version 150
in highp vec4 fColor;
in highp vec4 m_vertex;

uniform highp vec4 clipPlane;
uniform highp vec4 pointPlane;
uniform highp float rendering_mode;

out highp vec4 out_color;

void main(void)
{
  // onPlane == 1: inside clipping plane, should be solid;
  // onPlane == -1: outside clipping plane, should be transparent;
  // onPlane == 0: on clipping plane, whatever;
  float onPlane = sign(dot((m_vertex.xyz-pointPlane.xyz), clipPlane.xyz));

  // rendering_mode == -1: draw both inside and outside;
  // rendering_mode == 0: draw inside only;
  // rendering_mode == 1: draw outside only;
  if (rendering_mode == (onPlane+1)/2) {
    // discard other than the corresponding half when rendering
    discard;
  }

  out_color = fColor;
}