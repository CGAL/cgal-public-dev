// Copyright (c) 2018  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>

#ifndef CGAL_BASIC_SHADERS_H
#define CGAL_BASIC_SHADERS_H

#include <CGAL/license/GraphicsView.h>

namespace CGAL
{

//------------------------------------------------------------------------------
const char VERTEX_SOURCE_COLOR[]=R"DELIM(
#version 150
in highp vec4 P;
in highp vec3 N;
in highp vec3 Color;

uniform highp mat4 u_Mvp;
uniform highp mat4 u_Mv;
uniform highp float u_PointSize;

out highp vec4 vs_fP; // view space position
out highp vec4 ls_fP; // local space position 
out highp vec3 fN;
out highp vec4 fColor;

void main(void)
{
  ls_fP = P;
  vs_fP = u_Mv * P;

  fN = mat3(u_Mv)* N;
  fColor = vec4(Color, 1.0);

  gl_Position = u_Mvp * P;
  gl_PointSize = u_PointSize;
}
)DELIM";

const char FRAGMENT_SOURCE_COLOR[]=R"DELIM(
#version 150
in highp vec4 vs_fP;
in highp vec4 ls_fP;
in highp vec3 fN;
in highp vec4 fColor;

uniform highp vec4  u_LightPos;
uniform highp vec4  u_LightDiff;
uniform highp vec4  u_LightSpec;
uniform highp vec4  u_LightAmb;
uniform highp float u_SpecPower;

uniform highp vec4  u_ClipPlane;
uniform highp vec4  u_PointPlane;
uniform highp float u_RenderingMode;
uniform highp float u_RenderingTransparency;

out highp vec4 out_color;

void main(void)
{
  highp vec3 L = u_LightPos.xyz - vs_fP.xyz;
  highp vec3 V = -vs_fP.xyz;

  highp vec3 N = normalize(fN);
  L = normalize(L);
  V = normalize(V);

  highp vec3 R = reflect(-L, N);
  highp vec4 diffuse = vec4(max(dot(N,L), 0.0) * u_LightDiff.rgb * fColor.rgb, 1.0);
  highp vec4 ambient = vec4(u_LightAmb.rgb * fColor.rgb, 1.0);
  highp vec4 specular = pow(max(dot(R,V), 0.0), u_SpecPower) * u_LightSpec;

  // onPlane == 1: inside clipping plane, should be solid;
  // onPlane == -1: outside clipping plane, should be transparent;
  // onPlane == 0: on clipping plane, whatever;
  float onPlane = sign(dot((ls_fP.xyz-u_PointPlane.xyz), u_ClipPlane.xyz));

  // rendering_mode == -1: draw all solid;
  // rendering_mode == 0: draw solid only;
  // rendering_mode == 1: draw transparent only;
  if (u_RenderingMode == (onPlane+1)/2) {
    // discard other than the corresponding half when rendering
    discard;
  }

  // draw corresponding part
  out_color = u_RenderingMode < 1 ? (diffuse + ambient) :
                      vec4(diffuse.rgb + ambient.rgb, u_RenderingTransparency);
}
)DELIM";

const char VERTEX_SOURCE_P_L[]=R"DELIM(
#version 150
in highp vec4 P;
in highp vec3 Color;

uniform highp mat4  u_Mvp;
uniform highp float u_PointSize;
uniform highp float u_UseGeometryShader;

out highp vec4 fColor; // use if not u_UseGeometryShader
out highp vec4 gColor; // else 
out highp vec4 ls_fP; 

void main(void)
{
  if(u_UseGeometryShader > 0.0)
  {
    gColor      = vec4(Color, 1.0);
    gl_Position = P;
  }
  else 
  {
    fColor = vec4(Color, 1.0);
    ls_fP  = P;

    gl_PointSize = u_PointSize;
    gl_Position  = u_Mvp * P;
  }
}
)DELIM";

// This two values are hardware specific (using a fragment shader with signed distance function instead of gs could be a better choice for cylinder edge and sphere vertex feature)
// GL_MAX_GEOMETRY_OUTPUT_VERTICES = 128
// GL_MAX_GEOMETRY_OUTPUT_COMPONENTS = 1024

// To have a defined behavior for the gs shader, we must stay under these bounds. 
// We can compute the max_vertices value like so : max_vertices = min(GL_MAX_GEOMETRY_OUTPUT_VERTICES, GL_MAX_GEOMETRY_OUTPUT_COMPONENTS / components per vertex) 
// Components are values passed to the fragment shader per vertex (below, we have fColor, ls_fP and gl_Position, 3 * vec4 = 12 components)
// To draw a UV sphere, we provide an even resolution value. In the following shader, the maximum resolution is 8 because it needs 72 vertices to be rendered
// and the max_vertices value is about : 1024 / 12 = 85 vertices.
// A resolution of 10 would necessit 110 (max 128) vertices, and 1320 (max 1024) components which go beyond the bounds.   
const char GEOMETRY_SOURCE_SPHERE[]=R"DELIM(
#version 150 
layout(points) in; 
layout(triangle_strip, max_vertices = 72) out; // max_vertices = (resolution+1) * 2 * latResolution 

#define PI 3.14159265358979323846

in highp vec4 gColor[];

uniform highp mat4 u_Mvp;
uniform highp float u_Radius;

out highp vec4 fColor;
out highp vec4 ls_fP;

void drawSphere(in vec4 center, in float radius, in float resolution)
{
  float latResolution = resolution*0.5;
  float stepTheta = PI/latResolution;
  float stepPhi = 2*PI/resolution;
  for(int i=0; i<latResolution; ++i)
  {
    float theta1 = stepTheta*i;
    float theta2 = stepTheta*(i+1);
    for(int j=0; j<=resolution; ++j)
    {
      float phi = stepPhi*j;
      float x1 = center.x + radius * sin(theta1) * cos(phi);
      float y1 = center.y + radius * sin(theta1) * sin(phi);
      float z1 = center.z + radius * cos(theta1);
      ls_fP = vec4(x1, y1, z1, 1.0);
      gl_Position = u_Mvp * ls_fP;
      EmitVertex();

      float x2 = center.x + radius * sin(theta2) * cos(phi);
      float y2 = center.y + radius * sin(theta2) * sin(phi);
      float z2 = center.z + radius * cos(theta2); 
      ls_fP = vec4(x2, y2, z2, 1.0);
      gl_Position = u_Mvp * ls_fP;
      EmitVertex();
    }
    EndPrimitive();
  }
}

void main(void)
{
  fColor = gColor[0];

  int resolution = 8;
  vec4 center = gl_in[0].gl_Position; 

  drawSphere(center, u_Radius, resolution);
}
)DELIM";

const char GEOMETRY_SOURCE_CYLINDER[]=R"DELIM(
#version 150 
layout(lines) in; 
layout(triangle_strip, max_vertices = 22) out; 

#define PI 3.14159265358979323846

in highp vec4 gColor[];

uniform highp mat4 u_Mvp;
uniform highp float u_Radius;

out highp vec4 fColor;
out highp vec4 ls_fP;

void drawCylinder(in vec3 u, in vec3 v, in vec4 bot, in vec4 top, in float radius, in float resolution)
{
  float step = 2*PI/resolution;
  for(int i=0; i<=resolution; ++i)
  {
    float theta = step*i;
    float cosf = radius*cos(theta);
    float sinf = radius*sin(theta);
    vec3 xAxis = cosf*u.xyz;
    vec3 yAxis = sinf*v.xyz;
    ls_fP = vec4(top.xyz+xAxis.xyz+yAxis.xyz, 1.0);
    gl_Position = u_Mvp * ls_fP;
    EmitVertex();
    ls_fP = vec4(bot.xyz+xAxis.xyz+yAxis.xyz, 1.0);
    gl_Position = u_Mvp * ls_fP;
    EmitVertex();
  }
  EndPrimitive(); 
}

void main(void)
{
  fColor = gColor[0];

  vec4 a = gl_in[0].gl_Position;
  vec4 b = gl_in[1].gl_Position;

  vec3 w = normalize(vec3(-b.z, b.x, b.y));
  
  vec3 n = normalize(vec3(b.x-a.x, b.y-a.y, b.z-a.z)); // compute normal 

  // Axis vectors 
  vec3 u = normalize(cross(n, w));
  vec3 v = normalize(cross(n, u));

  int resolution = 10;

  drawCylinder(u, v, a, b, u_Radius, resolution);
}
)DELIM";

const char FRAGMENT_SOURCE_P_L[]=R"DELIM(
#version 150
in highp vec4 fColor;
in highp vec4 ls_fP;

uniform highp vec4  u_ClipPlane;
uniform highp vec4  u_PointPlane;
uniform highp float u_RenderingMode;

out highp vec4 out_color;

void main(void)
{
  // onPlane == 1: inside clipping plane, should be solid;
  // onPlane == -1: outside clipping plane, should be transparent;
  // onPlane == 0: on clipping plane, whatever;
  float onPlane = sign(dot((ls_fP.xyz-u_PointPlane.xyz), u_ClipPlane.xyz));

  // rendering_mode == -1: draw both inside and outside;
  // rendering_mode == 0: draw inside only;
  // rendering_mode == 1: draw outside only;
  if (u_RenderingMode == (onPlane+1)/2) {
    // discard other than the corresponding half when rendering
    discard;
  }

  out_color = fColor;
}
)DELIM";

const char VERTEX_SOURCE_CLIPPING_PLANE[]=R"DELIM(
#version 150
in highp vec4 P;

uniform highp mat4 u_Vp;
uniform highp mat4 u_M;

void main(void)
{
  gl_Position = u_Vp * u_M * P;
}
)DELIM";

const char FRAGMENT_SOURCE_CLIPPING_PLANE[]=R"DELIM(
#version 150

out highp vec4 out_color;

void main(void)
{
  out_color = vec4(0.0, 0.0, 0.0, 1.0);
}
)DELIM";

const char VERTEX_SOURCE_LINE[]=R"DELIM(
#version 150
in highp vec3 P;
in highp vec3 Color;

out VS_OUT {
  highp vec4 color;
} vs_out; // vertex shader output

void main(void)
{
  vs_out.color = vec4(Color, 1.0);

  gl_Position = vec4(P, 1.0);
}
)DELIM";

const char GEOMETRY_SOURCE_ARROW[]=R"DELIM(
#version 150
layout(lines) in; 
layout(triangle_strip, max_vertices = 82) out; // max_vertices = resolution * 2 + 2 (cylinder) + resolution * 3 (disc) + resolution * 3 (cone)

#define PI 3.14159265358979323846

in VS_OUT {
  highp vec4 color; 
} gs_in[]; // geometry shader input

uniform mat4  u_Mvp;
uniform float u_SceneRadius;

out highp vec4 fColor; 

void drawTriangle(in vec4 v1, in vec4 v2, in vec4 v3)
{
  gl_Position = u_Mvp*v1;
  EmitVertex();
  gl_Position = u_Mvp*v2;
  EmitVertex();
  gl_Position = u_Mvp*v3;
  EmitVertex();
  EndPrimitive();
}

void drawTriangleFan(in vec3 u, in vec3 v, in vec4 center, in vec4 edge0, in float radius, in int resolution)
{
  float step = 2*PI/resolution;
  for(int i=0; i<resolution; ++i)
  {
    float theta = step*i;
    float cosf = radius*cos(theta);
    float sinf = radius*sin(theta);
    vec3 xAxis = cosf*u.xyz;
    vec3 yAxis = sinf*v.xyz;
    vec4 edge1 = vec4(edge0.xyz+xAxis.xyz+yAxis.xyz, 1.0);
    theta = step*(i+1);
    cosf = radius*cos(theta);
    sinf = radius*sin(theta);
    xAxis = cosf*u.xyz;
    yAxis = sinf*v.xyz;
    vec4 edge2 = vec4(edge0.xyz+xAxis.xyz+yAxis.xyz, 1.0);
    drawTriangle(center, edge1, edge2);
  } 
}

void drawDisc(in vec3 u, in vec3 v, in vec4 center, in float radius, in int resolution)
{
  drawTriangleFan(u, v, center, center, radius, resolution);
}

void drawCone(in vec3 u, in vec3 v, in vec3 n, in vec4 center, in float radius, in float height, in int resolution)
{
  drawTriangleFan(u, v, center, vec4(center.xyz-height*n.xyz, 1.0), radius, resolution);
}

void drawCylinder(in vec3 u, in vec3 v, in vec4 bot, in vec4 top, in float radius, in float resolution)
{
  float step = 2*PI/resolution;
  for(int i=0; i<=resolution; ++i)
  {
    float theta = step*i;
    float cosf = radius*cos(theta);
    float sinf = radius*sin(theta);
    vec3 xAxis = cosf*u.xyz;
    vec3 yAxis = sinf*v.xyz;
    gl_Position = u_Mvp * vec4(top.xyz+xAxis.xyz+yAxis.xyz, 1.0);
    EmitVertex();
    gl_Position = u_Mvp * vec4(bot.xyz+xAxis.xyz+yAxis.xyz, 1.0);
    EmitVertex();
  }
  EndPrimitive(); 
}

void main(void)
{
  fColor = gs_in[0].color;

  vec4 a = gl_in[0].gl_Position;
  vec4 b = gl_in[1].gl_Position;

  vec3 w = normalize(vec3(-b.z, b.x, b.y));
  
  vec3 n = normalize(vec3(b.x-a.x, b.y-a.y, b.z-a.z)); // compute normal 

  // Axis vectors 
  vec3 u = normalize(cross(n, w));
  vec3 v = normalize(cross(n, u));

  float radius = 0.013 * u_SceneRadius;
  float height = 0.035 * u_SceneRadius;
  int resolution = 10;

  vec4 c = vec4(b.xyz-height*n.xyz, 1.0);
  drawDisc(u, v, c, radius, resolution);
  drawCone(u, v, n, b, radius, height, resolution);
  drawCylinder(u, v, a, c, radius*0.5, resolution);
}
)DELIM";

const char GEOMETRY_SOURCE_LINE[]=R"DELIM(
#version 150
layout(lines) in; 
layout(line_strip, max_vertices = 2) out; 

in VS_OUT {
  highp vec4 color;
} gs_in[]; // geometry shader input

uniform highp mat4 u_Mvp;

out highp vec4 fColor; 

void main(void)
{
  fColor = gs_in[0].color;

  gl_Position = u_Mvp * gl_in[0].gl_Position;
  EmitVertex(); 
  gl_Position = u_Mvp * gl_in[1].gl_Position; 
  EmitVertex(); 

  EndPrimitive(); 
}
)DELIM";

const char FRAGMENT_SOURCE_LINE[]=R"DELIM(
#version 150

in highp vec4 fColor;

out highp vec4 out_color;

void main(void)
{
  out_color = fColor;
}
)DELIM";

const char VERTEX_SOURCE_NORMAL[]=R"DELIM(
#version 150
in highp vec4 P;
in highp vec3 N;
 
uniform highp mat4  u_Mv;
uniform highp vec4  u_Color;
uniform highp float u_UseMonoColor;

out VS_OUT {
  highp vec4  color;
  highp vec3  normal;
} vs_out; // vertex shader output

void main(void)
{
  if (u_UseMonoColor > 0.0)
  {
    vs_out.color = u_Color; 
  }
  else 
  {
    vs_out.color = vec4(abs(normalize(N)), 1); 
  }

  mat3 normalMatrix = mat3(transpose(inverse(u_Mv)));
  vs_out.normal = normalize(vec3(vec4(normalMatrix * N, 0.0)));
  
  gl_Position = u_Mv * P; 
}
)DELIM";

const char GEOMETRY_SOURCE_NORMAL[]=R"DELIM(
#version 150 
layout (triangles) in;
layout (line_strip, max_vertices = 6) out;

in VS_OUT {
  highp vec4  color;
  highp vec3  normal;
} gs_in[]; // geometry shader input

uniform highp mat4  u_Projection;
uniform highp float u_Factor;
uniform highp float u_SceneRadius;

out GS_OUT {
  highp vec4 color;
} gs_out; // geometry shader output

void GenerateLine(int index)
{
  gs_out.color = gs_in[index].color; 
  
  gl_Position = u_Projection * gl_in[index].gl_Position;
  EmitVertex();

  gl_Position = u_Projection * (gl_in[index].gl_Position + vec4(gs_in[index].normal, 0.0) * u_SceneRadius * u_Factor);
  EmitVertex();

  EndPrimitive();
}

void main()
{
  GenerateLine(0); // first vertex normal
  GenerateLine(1); // second vertex normal
  GenerateLine(2); // third vertex normal
}
)DELIM";

const char FRAGMENT_SOURCE_NORMAL[]=R"DELIM(
#version 150

in GS_OUT {
  highp vec4 color;
} fs_in; // fragment shader input

out highp vec4 out_color;

void main()
{
  out_color = fs_in.color;
}
)DELIM";

//------------------------------------------------------------------------------
//  compatibility shaders

const char VERTEX_SOURCE_COLOR_COMP[]=R"DELIM(
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
)DELIM";

const char FRAGMENT_SOURCE_COLOR_COMP[]=R"DELIM(
varying highp vec4 fP;
varying highp vec3 fN;
varying highp vec4 fColor;

uniform highp vec4 u_LightPos;
uniform highp vec4 u_LightDiff;
uniform highp vec4 u_LightSpec;
uniform highp vec4 u_LightAmb;
uniform highp float u_SpecPower ;

void main(void)
{
  highp vec3 L = u_LightPos.xyz - fP.xyz;
  highp vec3 V = -fP.xyz;

  highp vec3 N = normalize(fN);
  L = normalize(L);
  V = normalize(V);

  highp vec3 R = reflect(-L, N);
  highp vec4 diffuse = max(dot(N,L), 0.0) * u_LightDiff * fColor;
  highp vec4 specular = pow(max(dot(R,V), 0.0), u_SpecPower) * u_LightSpec;

  gl_FragColor = u_LightAmb*fColor + diffuse;
}
)DELIM";

const char VERTEX_SOURCE_P_L_COMP[]=R"DELIM(
varying highp vec4 vertex;
varying highp vec3 color;

uniform highp mat4 mvp_matrix;
uniform highp float point_size;

varying highp vec4 fColor;

void main(void)
{
  gl_PointSize = point_size;
  fColor = vec4(color, 1.0);
  gl_Position = mvp_matrix * vertex;
}
)DELIM";

const char FRAGMENT_SOURCE_P_L_COMP[]=R"DELIM(
varying highp vec4 fColor;
void main(void)
{
  gl_FragColor = fColor;
}
)DELIM";

/* const char vertex_source_clipping_plane_comp[]=R"DELIM(
attribute highp vec4 vertex;

uniform highp mat4 vp_matrix;
uniform highp mat4 m_matrix;

void main(void)
{
  gl_Position = vp_matrix * m_matrix * vertex;
}
)DELIM";

const char fragment_source_clipping_plane_comp[]=R"DELIM(
out highp vec4 out_color;
void main(void)
{
  out_color = vec4(0.0, 0.0, 0.0, 1.0);
}
)DELIM";
*/

}

#endif // CGAL_BASIC_SHADERS_H
