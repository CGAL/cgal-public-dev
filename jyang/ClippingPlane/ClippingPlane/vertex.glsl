attribute highp vec4 a_position;
attribute highp vec2 a_texcoord;
attribute highp vec3 a_normal;

uniform highp mat4 u_mMatrix;
uniform highp mat4 u_pMatrix;
uniform highp mat4 u_vMatrix;

varying highp vec4 v_position;
varying highp vec2 v_texcoord;
varying highp vec3 v_normal;

varying highp vec3 va_clipPlane;
varying highp vec4 va_position;

void main(void) {

	vec4 clipPlane = vec4(1.0, 0.0, 0.0, 0.0);

	v_position = u_vMatrix * u_mMatrix * a_position;
	v_texcoord = a_texcoord;
	v_normal = normalize(vec3(u_vMatrix * u_mMatrix * vec4(a_normal, 0.0)));

	va_clipPlane = clipPlane;
	va_position = a_position;

	gl_Position = u_pMatrix * u_vMatrix * u_mMatrix * a_position;
}