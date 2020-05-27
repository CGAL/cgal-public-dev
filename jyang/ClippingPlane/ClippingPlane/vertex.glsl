attribute highp vec4 a_position;
attribute highp vec2 a_texcoord;
attribute highp vec3 a_normal;

uniform highp mat4 u_mMatrix;
uniform highp mat4 u_pMatrix;
uniform highp mat4 u_vMatrix;

varying highp vec4 v_position;
varying highp vec2 v_texcoord;
varying highp vec3 v_normal;

void main(void) {

	v_position = u_vMatrix * u_mMatrix * a_position;
	v_texcoord = a_texcoord;
	v_normal = normalize(vec3(u_vMatrix * u_mMatrix * vec4(a_normal, 0.0)));

	gl_Position = u_pMatrix * u_vMatrix * u_mMatrix * a_position;
}