#version 440

layout(location = 0) in vec4 position;
layout(location = 1) in vec3 color;

layout(location = 0) out vec3 v_color;

layout(std140, binding = 0) uniform buf {
    mat4 mvp;
} ubuf;

out gl_PerVertex { 
	vec4 gl_Position;  
	float gl_PointSize;
};

void main()
{
	gl_PointSize = 10.0f;
    v_color = vec3(.345f,0.0f,0.0f);
    gl_Position = ubuf.mvp * position;
}
