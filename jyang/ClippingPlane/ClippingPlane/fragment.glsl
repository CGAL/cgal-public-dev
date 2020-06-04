varying highp vec4 v_position;
varying highp vec2 v_texcoord;
varying highp vec3 v_normal;

varying highp vec3 va_clipPlane;
varying highp vec4 va_position;

uniform sampler2D u_texture;

void main(void) {
	// init vector
	vec4 eyePosition = vec4(0.0, 0.0, 0.0, 1.0);
	vec4 lightPosition = vec4(0.0, 0.0, 0.0, 1.0);

	// init parameter
	float dist = length(v_position.xyz - eyePosition.xyz);
	float specularFactor = 10.0;
	float ambientFactor = 0.2;
	float lightPower = 5.0;

	// get vector
	vec3 eyeVec = normalize(v_position.xyz - eyePosition.xyz);
	vec3 lightVec = normalize(v_position.xyz - lightPosition.xyz);
	vec3 reflectVec = normalize(reflect(lightVec, v_normal));

	// calculate color
	vec4 diffuse = vec4(texture2D(u_texture, v_texcoord).rgb * lightPower * max(0.0, dot(v_normal, -lightVec)) / (1.0 + 0.25 * pow(dist, 2.0)), texture2D(u_texture, v_texcoord).a);
	vec4 ambient = vec4(texture2D(u_texture, v_texcoord).rgb * ambientFactor, texture2D(u_texture, v_texcoord).a) ;
	vec4 specular = vec4(1.0, 1.0, 1.0, 1.0) * lightPower * pow(max(0.0, dot(reflectVec, -eyeVec)), specularFactor) / (1.0 + 0.25 * pow(dist, 2.0));

	// update inside/outside fade factor
	//float fadeFactor = (dot(va_position.xyz, va_clipPlane) / (1e-5 + abs(dot(va_position.xyz, va_clipPlane))) + 1.0) * 0.4 + 0.2; // [0.2 or 1.0] -> [outside or inside]
	//diffuse = diffuse * vec4(1.0, 1.0, 1.0, fadeFactor);
	//ambient = ambient * vec4(1.0, 1.0, 1.0, fadeFactor);
	//specular = specular * vec4(1.0, 1.0, 1.0, fadeFactor);

	// 
	float onPlane = sign(dot(va_position.xyz, va_clipPlane));
	diffuse = abs(onPlane) * diffuse * vec4(1.0, 1.0, 1.0, 0.6 + 0.4*onPlane) + (1 - abs(onPlane)) * vec4(1.0, 1.0, 1.0, 1.0);
	ambient = abs(onPlane) * ambient * vec4(1.0, 1.0, 1.0, 0.6 + 0.4*onPlane) + (1 - abs(onPlane)) * vec4(1.0, 1.0, 1.0, 1.0);
	specular = abs(onPlane) * specular * vec4(1.0, 1.0, 1.0, 0.6 + 0.4*onPlane) + (1 - abs(onPlane)) * vec4(1.0, 1.0, 1.0, 1.0);

	// return
	gl_FragColor = diffuse + ambient + specular;
}