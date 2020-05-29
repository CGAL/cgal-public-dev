varying highp vec4 v_position;
varying highp vec2 v_texcoord;
varying highp vec3 v_normal;

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

	// init plane
	vec3 plane = vec3(1.0, 0.0, 0.0);

	// get vector
	vec3 eyeVec = normalize(v_position.xyz - eyePosition.xyz);
	vec3 lightVec = normalize(v_position.xyz - lightPosition.xyz);
	vec3 reflectVec = normalize(reflect(lightVec, v_normal));

	// calculate color
	vec4 diffuse = vec4(texture2D(u_texture, v_texcoord).rgb * lightPower * max(0.0, dot(v_normal, -lightVec)) / (1.0 + 0.25 * pow(dist, 2.0)), 1.0);
	vec4 ambient = vec4(texture2D(u_texture, v_texcoord).rgb * ambientFactor, 1.0) ;
	vec4 specular = vec4(1.0, 1.0, 1.0, 1.0) * lightPower * pow(max(0.0, dot(reflectVec, -eyeVec)), specularFactor) / (1.0 + 0.25 * pow(dist, 2.0));

	// update inside/outside fade factor
	float fadeFactor = (dot(v_position.xyz, plane) / abs(dot(v_position.xyz, plane)) + 1.0) * 0.35 + 0.3; // [0.3, 1.0] -> [outside,inside]
	diffuse = diffuse * vec4(1.0, 1.0, 1.0, fadeFactor);
	ambient = ambient * vec4(1.0, 1.0, 1.0, fadeFactor);
	specular = specular * vec4(1.0, 1.0, 1.0, fadeFactor);

	// return
	gl_FragColor = diffuse + ambient + specular;
}