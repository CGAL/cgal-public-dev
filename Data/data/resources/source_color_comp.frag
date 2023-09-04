#version 440
layout(location = 0) in highp vec4 fP;
layout(location = 1) in highp vec3 fN;
layout(location = 2) in highp vec4 fColor;

layout(std140, binding = 0) uniform buf {
	highp vec4 light_pos;
	highp vec4 light_diff;
	highp vec4 light_spec;
	highp vec4 light_amb;
	highp float spec_power ;
} ubuf;

layout(location = 0) out vec4 fragColor;

void main(void)
{
  highp vec3 L = ubuf.light_pos.xyz - fP.xyz;
  highp vec3 V = -fP.xyz;

  highp vec3 N = normalize(fN);
  L = normalize(L);
  V = normalize(V);

  highp vec3 R = reflect(-L, N);
  highp vec4 diffuse = max(dot(N,L), 0.0) * ubuf.light_diff * fColor;
  highp vec4 specular = pow(max(dot(R,V), 0.0), ubuf.spec_power) * ubuf.light_spec;

  fragColor = ubuf.light_amb*fColor + diffuse;
}