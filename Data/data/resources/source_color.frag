#version 150
in highp vec4 fP;
in highp vec3 fN;
in highp vec4 fColor;
in highp vec4 m_vertex;

uniform highp vec4 light_pos;
uniform highp vec4 light_diff;
uniform highp vec4 light_spec;
uniform highp vec4 light_amb;
uniform highp float spec_power;

uniform highp vec4 clipPlane;
uniform highp vec4 pointPlane;
uniform highp float rendering_mode;
uniform highp float rendering_transparency;

out highp vec4 out_color;

void main(void)
{
  highp vec3 L = light_pos.xyz - fP.xyz;
  highp vec3 V = -fP.xyz;

  highp vec3 N = normalize(fN);
  L = normalize(L);
  V = normalize(V);

  highp vec3 R = reflect(-L, N);
  highp vec4 diffuse = vec4(max(dot(N,L), 0.0) * light_diff.rgb * fColor.rgb, 1.0);
  highp vec4 ambient = vec4(light_amb.rgb * fColor.rgb, 1.0);
  highp vec4 specular = pow(max(dot(R,V), 0.0), spec_power) * light_spec;

  // onPlane == 1: inside clipping plane, should be solid;
  // onPlane == -1: outside clipping plane, should be transparent;
  // onPlane == 0: on clipping plane, whatever;
  float onPlane = sign(dot((m_vertex.xyz-pointPlane.xyz), clipPlane.xyz));

  // rendering_mode == -1: draw all solid;
  // rendering_mode == 0: draw solid only;
  // rendering_mode == 1: draw transparent only;
  if (rendering_mode == (onPlane+1)/2) {
    // discard other than the corresponding half when rendering
    discard;
  }

  // draw corresponding part
  out_color = rendering_mode < 1 ? (diffuse + ambient) :
                      vec4(diffuse.rgb + ambient.rgb, rendering_transparency);
}